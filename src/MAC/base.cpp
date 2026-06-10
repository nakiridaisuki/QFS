#include "MAC/base.h"
#include <algorithm>
#include <math.h>
#include <vector>

// PUBLIC
MACSimulatorBase::MACSimulatorBase(int width, int height)
    : BaseSimulator(width, height) {

    u.resize((nx + 1) * ny, 0.0);
    v.resize(nx * (ny + 1), 0.0);
    u_old.resize((nx + 1) * ny, 0.0);
    v_old.resize(nx * (ny + 1), 0.0);
    weight_u.resize((nx + 1) * ny, 0.0);
    weight_v.resize(nx * (ny + 1), 0.0);

    cell_type.resize(nx * ny, 0);
    current_count.resize(nx * ny, 0);
    fluid_map.resize(nx * ny, 0);
    triplets.reserve(nx * ny * 5);

    G     = 150.f;
    Sigma = 0.0f;

    phi[0].resize(nx * ny, 0.0f);
    phi[1].resize(nx * ny, 0.0f);
    nx_n.resize(nx * ny, 0.0f);
    ny_n.resize(nx * ny, 0.0f);
    kappa.resize(nx * ny, 0.0f);

    valid_u.resize((nx + 1) * ny, 0);
    valid_v.resize(nx * (ny + 1), 0);

    // solver.setMaxIterations(80);
    solver.setTolerance(0.1);
    // Eigen::setNbThreads(6);
}

std::vector<Line> MACSimulatorBase::getLines() const {
    std::vector<Line> lines;

    for (int i = 1; i < nx; i++) {
        lines.push_back({(float)i, 0, (float)i, (float)ny});
    }
    for (int i = 1; i < ny; i++) {
        lines.push_back({0, (float)i, (float)nx, (float)i});
    }
    return lines;
}

void MACSimulatorBase::setBoundaries(
    std::vector<float> &ufield, std::vector<float> &vfield
) {

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        ufield[IX_u(0, j)]  = 0.0;
        ufield[IX_u(nx, j)] = 0.0;
    }
#pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        vfield[IX_v(i, 0)]  = 0.0;
        vfield[IX_v(i, ny)] = 0.0;
    }
}

void MACSimulatorBase::project() {
    int N = nx * ny;

    std::fill(fluid_map.begin(), fluid_map.end(), -1);
    int fluid_count = 0;
    for (int i = 0; i < N; i++) {
        if (cell_type[i] == 1) {
            fluid_map[i] = fluid_count++;
        }
    }

    if (fluid_count == 0)
        return;

    Eigen::VectorXf div(fluid_count);
    div.setZero();
    triplets.clear();

    int max_threads = omp_get_max_threads();
    static std::vector<std::vector<Eigen::Triplet<float>>> thread_triplets(
        max_threads
    );

#pragma omp parallel
    {
        int tid       = omp_get_thread_num();
        auto &local_t = thread_triplets[tid];
        local_t.clear();

#pragma omp for
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = IX(i, j);

                if (cell_type[idx] == 0) {
                    continue;
                }

                int row  = fluid_map[idx];
                float d  = u[IX_u(i + 1, j)] -
                           u[IX_u(i, j)] +
                           v[IX_v(i, j + 1)] -
                           v[IX_v(i, j)];
                div[row] = -d;

                int count = 0;
                if (i > 0) {
                    if (cell_type[IX(i - 1, j)] != 0) {
                        local_t.push_back(
                            {row, fluid_map[IX(i - 1, j)], -1.0f}
                        );
                    }
                    count++;
                }
                if (i < nx - 1) {
                    if (cell_type[IX(i + 1, j)] != 0) {
                        local_t.push_back(
                            {row, fluid_map[IX(i + 1, j)], -1.0f}
                        );
                    }
                    count++;
                }
                if (j > 0) {
                    if (cell_type[IX(i, j - 1)] != 0) {
                        local_t.push_back(
                            {row, fluid_map[IX(i, j - 1)], -1.0f}
                        );
                    }
                    count++;
                }
                if (j < ny - 1) {
                    if (cell_type[IX(i, j + 1)] != 0) {
                        local_t.push_back(
                            {row, fluid_map[IX(i, j + 1)], -1.0f}
                        );
                    }
                    count++;
                }

                local_t.push_back({row, row, (float)count});
            }
        }
    }
    for (const auto &local_t : thread_triplets)
        triplets.insert(triplets.end(), local_t.begin(), local_t.end());

    // Eigen solvers
    Eigen::SparseMatrix<float> A(fluid_count, fluid_count);
    A.setFromTriplets(triplets.begin(), triplets.end());

    solver.compute(A);
    Eigen::VectorXf pressure = solver.solve(div);

    PCGItertimes = solver.iterations();

    auto get_pressure = [&](int i, int j) {
        int idx   = IX(i, j);
        int f_idx = fluid_map[idx];
        if (f_idx == -1)
            return 0.0f;
        return pressure[f_idx];
    };

#pragma omp parallel
    {
#pragma omp for nowait
        for (int j = 0; j < ny; j++) {
            for (int i = 1; i < nx; i++) {
                float p1 = get_pressure(i, j);
                float p2 = get_pressure(i - 1, j);
                u[IX_u(i, j)] -= (p1 - p2);
            }
        }

#pragma omp for
        for (int j = 1; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                float p1 = get_pressure(i, j);
                float p2 = get_pressure(i, j - 1);
                v[IX_v(i, j)] -= (p1 - p2);
            }
        }
    }
}

void MACSimulatorBase::velExtrapolation() {
    const int ext_layers = 3;
    const int AIR        = 99;

#pragma omp parallel
    {
#pragma omp for nowait
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i <= nx; i++) {
                bool fluid_left  = (i > 0) && (cell_type[IX(i - 1, j)] == 1);
                bool fluid_right = (i < nx) && (cell_type[IX(i, j)] == 1);
                if (fluid_left || fluid_right) {
                    valid_u[IX_u(i, j)] = 1;
                } else
                    valid_u[IX_u(i, j)] = AIR;
            }
        }

#pragma omp for
        for (int j = 0; j <= ny; j++) {
            for (int i = 0; i < nx; i++) {
                bool fluid_bottom = (j > 0) && (cell_type[IX(i, j - 1)] == 1);
                bool fluid_top    = (j < ny) && (cell_type[IX(i, j)] == 1);
                if (fluid_bottom || fluid_top) {
                    valid_v[IX_v(i, j)] = 1;
                } else
                    valid_v[IX_v(i, j)] = AIR;
            }
        }

        for (int iter = 1; iter <= ext_layers; iter++) {
#pragma omp for nowait
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i <= nx; i++) {
                    if (valid_u[IX_u(i, j)] == AIR) {
                        float sum = 0;
                        int count = 0;
                        if (i > 0 && valid_u[IX_u(i - 1, j)] <= iter) {
                            sum += u[IX_u(i - 1, j)];
                            count++;
                        }
                        if (i < nx && valid_u[IX_u(i + 1, j)] <= iter) {
                            sum += u[IX_u(i + 1, j)];
                            count++;
                        }
                        if (j > 0 && valid_u[IX_u(i, j - 1)] <= iter) {
                            sum += u[IX_u(i, j - 1)];
                            count++;
                        }
                        if (j < ny - 1 && valid_u[IX_u(i, j + 1)] <= iter) {
                            sum += u[IX_u(i, j + 1)];
                            count++;
                        }

                        if (count > 0) {
                            u[IX_u(i, j)]       = sum / (float)count;
                            valid_u[IX_u(i, j)] = iter + 1;
                        }
                    }
                }
            }

#pragma omp for
            for (int j = 0; j <= ny; j++) {
                for (int i = 0; i < nx; i++) {
                    if (valid_v[IX_v(i, j)] == AIR) {
                        float sum = 0;
                        int count = 0;
                        if (i > 0 && valid_v[IX_v(i - 1, j)] <= iter) {
                            sum += v[IX_v(i - 1, j)];
                            count++;
                        }
                        if (i < nx - 1 && valid_v[IX_v(i + 1, j)] <= iter) {
                            sum += v[IX_v(i + 1, j)];
                            count++;
                        }
                        if (j > 0 && valid_v[IX_v(i, j - 1)] <= iter) {
                            sum += v[IX_v(i, j - 1)];
                            count++;
                        }
                        if (j < ny && valid_v[IX_v(i, j + 1)] <= iter) {
                            sum += v[IX_v(i, j + 1)];
                            count++;
                        }

                        if (count > 0) {
                            v[IX_v(i, j)]       = sum / (float)count;
                            valid_v[IX_v(i, j)] = iter + 1;
                        }
                    }
                }
            }
        }
    }
}

void MACSimulatorBase::applyGravity(float dt) {
#pragma omp parallel for
    for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx; i++) {
            if ((j < ny && cell_type[IX(i, j)] == 1) ||
                (j > 0 && cell_type[IX(i, j - 1)] == 1)) {
                v[IX_v(i, j)] += G * dt;
            }
        }
    }
}

float MACSimulatorBase::bilerp(
    const std::vector<float> &field, int w, int h, float x, float y
) const {
    x = std::max(0.0f, std::min((float)w - 1.001f, x));
    y = std::max(0.0f, std::min((float)h - 1.001f, y));

    int i    = (int)x;
    int j    = (int)y;
    float fx = x - i;
    float fy = y - j;

    float c00 = field[i + j * w];
    float c10 = field[(i + 1) + j * w];
    float c01 = field[i + (j + 1) * w];
    float c11 = field[(i + 1) + (j + 1) * w];

    return (c00 * (1 - fx) + c10 * fx) * (1 - fy) +
           (c01 * (1 - fx) + c11 * fx) * fy;
}

void MACSimulatorBase::bidistri(
    std::vector<float> &field, int w, int h, float x, float y, float value
) {
    x = std::max(0.0f, std::min((float)w - 1.001f, x));
    y = std::max(0.0f, std::min((float)h - 1.001f, y));

    int i    = (int)x;
    int j    = (int)y;
    float fx = x - i;
    float fy = y - j;

    float c0 = value * (1 - fx);
    float c1 = value * fx;

    float c00 = c0 * (1 - fy);
    float c01 = c0 * fy;
    float c10 = c1 * (1 - fy);
    float c11 = c1 * fy;

    field[i + j * w] += c00;
    field[(i + 1) + j * w] += c10;
    field[i + (j + 1) * w] += c01;
    field[(i + 1) + (j + 1) * w] += c11;
}
