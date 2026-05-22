#include "MAC.h"
#include <algorithm>
#include <cstdlib>
#include <math.h>

// PUBLIC
MACSimulator::MACSimulator(int width, int height, int fps = 60)
    : nx(width), ny(height), fps(fps) {
    u.resize((nx + 1) * ny, 0.0);
    v.resize(nx * (ny + 1), 0.0);
    u_old.resize((nx + 1) * ny, 0.0);
    v_old.resize(nx * (ny + 1), 0.0);
    weight_u.resize((nx + 1) * ny, 0.0);
    weight_v.resize(nx * (ny + 1), 0.0);

    cell_type.resize(nx * ny, 0);
    particles_count.resize(nx * ny, 0);
    current_count.resize(nx * ny, 0);
    fluid_map.resize(nx * ny, 0);
    triplets.reserve(nx * ny * 5);
    p.resize(nx * ny, 0.0);
    G = 150.f;

    solver.setMaxIterations(80);
    solver.setTolerance(1e-3);
    Eigen::setNbThreads(6);
}

void MACSimulator::update(float dt) {
    static int fps_cnt = 0;
    fps_cnt++;

    particleToGrid();

    u_old = u;
    v_old = v;

    setBoundaries(u_old, v_old);

    markFluidCells();
    applyGravity(dt);

    setBoundaries(u, v);
    project();
    setBoundaries(u, v);

    gridToParticle();

    resampleParticles();

    advectParticles(dt);
}

void MACSimulator::addWater(float x, float y, float radius) {
    float r2 = radius * radius;
    for (float px = x - radius; px < x + radius; px += 0.5) {
        for (float py = y - radius; py < y + radius; py += 0.5) {
            if (((px - x) * (px - x) + (py - y) * (py - y)) < r2)
                particles.push_back({px, py, 0.f, 0.f});
        }
    }
}
void MACSimulator::delWater(float x, float y, float radius) {
    float r2 = radius * radius;

    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](auto &p) {
            return ((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y)) < r2;
        });

    particles.erase(new_end, particles.end());
}

void MACSimulator::setBoundaries(
    std::vector<float> &ufield, std::vector<float> &vfield
) {
    for (int j = 0; j < ny; j++) {
        if (ufield[IX_u(0, j)] < 0.f)
            ufield[IX_u(0, j)] = 0.0;
        if (ufield[IX_u(nx, j)] > 0.f)
            ufield[IX_u(nx, j)] = 0.0;
    }
    for (int i = 0; i < nx; i++) {
        if (vfield[IX_v(i, 0)] < 0.f)
            vfield[IX_v(i, 0)] = 0.0;
        if (vfield[IX_v(i, ny)] > 0.f)
            vfield[IX_v(i, ny)] = 0.0;
    }
}

void MACSimulator::particleToGrid() {

    std::fill(u.begin(), u.end(), 0.f);
    std::fill(v.begin(), v.end(), 0.f);
    std::fill(weight_u.begin(), weight_u.end(), 0.f);
    std::fill(weight_v.begin(), weight_v.end(), 0.f);

    for (const auto &p : particles) {
        bidistri(u, nx + 1, ny, p.x, p.y - 0.5, p.u);
        bidistri(weight_u, nx + 1, ny, p.x, p.y - 0.5, 1);

        bidistri(v, nx, ny + 1, p.x - 0.5, p.y, p.v);
        bidistri(weight_v, nx, ny + 1, p.x - 0.5, p.y, 1);
    }

    for (int i = 0; i < u.size(); i++) {
        if (weight_u[i] > 1e-6)
            u[i] /= weight_u[i];
    }
    for (int i = 0; i < v.size(); i++) {
        if (weight_v[i] > 1e-6)
            v[i] /= weight_v[i];
    }
}

void MACSimulator::gridToParticle() {
    float flipRatio = 0.95f;

#pragma omp parallel for
    for (auto &p : particles) {
        float u_new_grid = bilerp(u, nx + 1, ny, p.x, p.y - 0.5);
        float v_new_grid = bilerp(v, nx, ny + 1, p.x - 0.5, p.y);

        float u_old_grid = bilerp(u_old, nx + 1, ny, p.x, p.y - 0.5);
        float v_old_grid = bilerp(v_old, nx, ny + 1, p.x - 0.5, p.y);

        float u_pic = u_new_grid;
        float v_pic = v_new_grid;

        float u_flip = p.u + (u_new_grid - u_old_grid);
        float v_flip = p.v + (v_new_grid - v_old_grid);

        p.u = flipRatio * u_flip + (1.f - flipRatio) * u_pic;
        p.v = flipRatio * v_flip + (1.f - flipRatio) * v_pic;
    }
}

void MACSimulator::advectParticles(float dt) {
#pragma omp parallel for
    for (auto &p : particles) {
        p.x += p.u * dt;
        p.y += p.v * dt;

        p.x = std::max(0.001f, std::min(float(nx) - 0.001f, p.x));
        p.y = std::max(0.001f, std::min(float(ny) - 0.001f, p.y));
    }
}

void MACSimulator::project() {
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

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = IX(i, j);

            if (cell_type[idx] == 0) {
                continue;
            }

            int row = fluid_map[idx];
            float d = u[IX_u(i + 1, j)] -
                      u[IX_u(i, j)] +
                      v[IX_v(i, j + 1)] -
                      v[IX_v(i, j)];
            div[row] = -d;

            int count = 0;
            if (i > 0) {
                if (cell_type[IX(i - 1, j)] != 0) {
                    triplets.push_back({row, fluid_map[IX(i - 1, j)], -1.0f});
                }
                count++;
            }
            if (i < nx - 1) {
                if (cell_type[IX(i + 1, j)] != 0) {
                    triplets.push_back({row, fluid_map[IX(i + 1, j)], -1.0f});
                }
                count++;
            }
            if (j > 0) {
                if (cell_type[IX(i, j - 1)] != 0) {
                    triplets.push_back({row, fluid_map[IX(i, j - 1)], -1.0f});
                }
                count++;
            }
            if (j < ny - 1) {
                if (cell_type[IX(i, j + 1)] != 0) {
                    triplets.push_back({row, fluid_map[IX(i, j + 1)], -1.0f});
                }
                count++;
            }

            triplets.push_back({row, row, (float)count});
        }
    }
    // Eigen solvers
    Eigen::SparseMatrix<float> A(fluid_count, fluid_count);
    A.setFromTriplets(triplets.begin(), triplets.end());

    solver.compute(A);
    Eigen::VectorXf pressure = solver.solve(div);

    std::fill(p.begin(), p.end(), 0.0);
    for (int i = 0; i < N; i++)
        if (fluid_map[i] != -1)
            p[i] = pressure[fluid_map[i]];

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 1; i < nx; i++) {
            u[IX_u(i, j)] -= (p[IX(i, j)] - p[IX(i - 1, j)]);
        }
    }

#pragma omp parallel for
    for (int j = 1; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            v[IX_v(i, j)] -= (p[IX(i, j)] - p[IX(i, j - 1)]);
        }
    }
}

void MACSimulator::applyGravity(float dt) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            if (cell_type[IX(i, j)] == 1 ||
                (j > 0 && cell_type[IX(i, j - 1)] == 1)) {
                v[IX_v(i, j)] += G * dt;
            }
        }
    }
}

void MACSimulator::markFluidCells() {
    std::fill(cell_type.begin(), cell_type.end(), 0);

    for (const auto &p : particles) {
        int i = (int)p.x;
        int j = (int)p.y;
        if (i >= 0 && i < nx && j >= 0 && j < ny) {
            cell_type[IX(i, j)] = 1; // 1 = 液體
        }
    }
}

void MACSimulator::resampleParticles() {
    const int target_ppc = 4;
    const int max_ppc = 8;
    const int min_ppc = 3;

    std::fill(particles_count.begin(), particles_count.end(), 0);
    for (const auto &p : particles) {
        int i = std::clamp((int)p.x, 0, nx - 1);
        int j = std::clamp((int)p.y, 0, ny - 1);
        particles_count[IX(i, j)]++;
    }

    std::fill(current_count.begin(), current_count.end(), 0);
    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](const auto &p) {
            int i = std::clamp((int)p.x, 0, nx - 1);
            int j = std::clamp((int)p.y, 0, ny - 1);
            int idx = IX(i, j);

            current_count[idx]++;
            if (current_count[idx] > max_ppc)
                return true;

            return false;
        });
    particles.erase(new_end, particles.end());

    for (int j = 1; j < ny - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            int idx = IX(i, j);

            bool is_inter =
                (cell_type[idx] &&
                 cell_type[IX(i + 1, j)] &&
                 cell_type[IX(i - 1, j)] &&
                 cell_type[IX(i, j + 1)] &&
                 cell_type[IX(i, j - 1)]);

            if (is_inter && particles_count[idx] < min_ppc) {
                int needed = target_ppc - particles_count[idx];

                for (int k = 0; k < needed; k++) {
                    float rx = i + ((float)rand() / RAND_MAX);
                    float ry = j + ((float)rand() / RAND_MAX);

                    float new_u = bilerp(u, nx + 1, ny, rx, ry - 0.5f);
                    float new_v = bilerp(v, nx, ny + 1, rx - 0.5f, ry);

                    particles.push_back({rx, ry, new_u, new_v});
                }
            }
        }
    }
}

float MACSimulator::bilerp(
    const std::vector<float> &field, int w, int h, float x, float y
) const {
    x = std::max(0.0f, std::min((float)w - 1.001f, x));
    y = std::max(0.0f, std::min((float)h - 1.001f, y));

    int i = (int)x;
    int j = (int)y;
    float fx = x - i;
    float fy = y - j;

    float c00 = field[i + j * w];
    float c10 = field[(i + 1) + j * w];
    float c01 = field[i + (j + 1) * w];
    float c11 = field[(i + 1) + (j + 1) * w];

    return (c00 * (1 - fx) + c10 * fx) * (1 - fy) +
           (c01 * (1 - fx) + c11 * fx) * fy;
}

void MACSimulator::bidistri(
    std::vector<float> &field, int w, int h, float x, float y, float value
) {
    x = std::max(0.0f, std::min((float)w - 1.001f, x));
    y = std::max(0.0f, std::min((float)h - 1.001f, y));

    int i = (int)x;
    int j = (int)y;
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
