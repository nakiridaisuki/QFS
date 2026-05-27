#include "MAC.h"
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <vector>

// PUBLIC
MACSimulator::MACSimulator(int width, int height) : nx(width), ny(height) {
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
    Sigma = 0.0f;

    solver.setMaxIterations(80);
    solver.setTolerance(1e-3);
    Eigen::setNbThreads(6);
}

void MACSimulator::update(float dt) {
    max_u = 0.0f;
    max_v = 0.0f;

    particleToGrid();

    u_old = u;
    v_old = v;

    setBoundaries(u_old, v_old);

    applyGravity(dt);
    applySurfaceTension(dt);

    setBoundaries(u, v);
    project();
    setBoundaries(u, v);

    velExtrapolation();

    gridToParticle();

    advectParticles(dt);
    markFluidCells();
    resampleParticles();
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
    markFluidCells();
}

std::vector<Line> MACSimulator::getLines() const {
    std::vector<Line> lines;

    for (int i = 1; i < nx; i++) {
        lines.push_back({(float)i, 0, (float)i, (float)ny});
    }
    for (int i = 1; i < ny; i++) {
        lines.push_back({0, (float)i, (float)nx, (float)i});
    }
    return lines;
}

void MACSimulator::setBoundaries(
    std::vector<float> &ufield, std::vector<float> &vfield
) {
    for (int j = 0; j < ny; j++) {
        // if (ufield[IX_u(0, j)] < 0.f)
        ufield[IX_u(0, j)] = 0.0;
        // if (ufield[IX_u(nx, j)] > 0.f)
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

void MACSimulator::velExtrapolation() {
    const int ext_layers = 3;
    const int AIR = 99;

    static std::vector<int> valid_u, valid_v;
    valid_u.assign((nx + 1) * ny, AIR);
    valid_v.assign(nx * (ny + 1), AIR);

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i <= nx; i++) {
            bool fluid_left = (i > 0) && (cell_type[IX(i - 1, j)] == 1);
            bool fluid_right = (i < nx) && (cell_type[IX(i, j)] == 1);
            if (fluid_left || fluid_right) {
                valid_u[IX_u(i, j)] = 1;
            }
        }
    }

#pragma omp parallel for
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i < nx; i++) {
            bool fluid_bottom = (j > 0) && (cell_type[IX(i, j - 1)] == 1);
            bool fluid_top = (j < ny) && (cell_type[IX(i, j)] == 1);
            if (fluid_bottom || fluid_top) {
                valid_v[IX_v(i, j)] = 1;
            }
        }
    }

    for (int iter = 1; iter <= ext_layers; iter++) {
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
                        u[IX_u(i, j)] = sum / (float)count;
                        valid_u[IX_u(i, j)] = iter + 1;
                    }
                }
            }
        }

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
                        v[IX_v(i, j)] = sum / (float)count;
                        valid_v[IX_v(i, j)] = iter + 1;
                    }
                }
            }
        }
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
    float local_max_u = 0.0f;
    float local_max_v = 0.0f;

#pragma omp parallel for reduction(max : local_max_u)                          \
    reduction(max : local_max_v)
    for (auto &p : particles) {
        p.x += p.u * dt;
        p.y += p.v * dt;

        p.x = std::max(0.001f, std::min(float(nx) - 0.001f, p.x));
        p.y = std::max(0.001f, std::min(float(ny) - 0.001f, p.y));

        local_max_u = std::max(std::abs(p.u), max_u);
        local_max_v = std::max(std::abs(p.v), max_v);
    }

    max_u = local_max_u;
    max_v = local_max_v;
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
    for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx; i++) {
            if ((j < ny && cell_type[IX(i, j)] == 1) ||
                (j > 0 && cell_type[IX(i, j - 1)] == 1)) {
                v[IX_v(i, j)] += G * dt;
            }
        }
    }
}

void MACSimulator::applySurfaceTension(float dt) {
    if (Sigma <= 0.0f)
        return;

    static std::vector<float> phi[2];
    phi[0].assign(nx * ny, 0.0f);
    phi[1].assign(nx * ny, 0.0f);

#pragma omp parallel for
    for (int i = 0; i < nx * ny; i++) {
        phi[0][i] = (float)cell_type[i];
    }

    for (int iter = 0; iter < 2; iter++) {
#pragma omp parallel for
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                phi[1 - iter][IX(i, j)] = (phi[iter][IX(i, j)] * 4.0f +
                                           phi[iter][IX(i - 1, j)] +
                                           phi[iter][IX(i + 1, j)] +
                                           phi[iter][IX(i, j - 1)] +
                                           phi[iter][IX(i, j + 1)]) /
                                          8.0f;
            }
        }
    }

    static std::vector<float> nx_n, ny_n;
    nx_n.assign(nx * ny, 0.0f);
    ny_n.assign(nx * ny, 0.0f);

#pragma omp parallel for
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            // 利用中心差分計算梯度
            float dx = (phi[0][IX(i + 1, j)] - phi[0][IX(i - 1, j)]) * 0.5f;
            float dy = (phi[0][IX(i, j + 1)] - phi[0][IX(i, j - 1)]) * 0.5f;
            float len =
                std::sqrt(dx * dx + dy * dy + 1e-8f); // 加上 1e-8 避免除以零
            nx_n[IX(i, j)] = dx / len;
            ny_n[IX(i, j)] = dy / len;
        }
    }

    static std::vector<float> kappa;
    kappa.assign(nx * ny, 0.0f);

#pragma omp parallel for
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            float div_x = (nx_n[IX(i + 1, j)] - nx_n[IX(i - 1, j)]) * 0.5f;
            float div_y = (ny_n[IX(i, j + 1)] - ny_n[IX(i, j - 1)]) * 0.5f;
            kappa[IX(i, j)] = -(div_x + div_y);
        }
    }

#pragma omp parallel for
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 2; i < nx - 1; ++i) { // 避開最外圈邊界
            // u 面的曲率取左右格子的平均
            float k = (kappa[IX(i, j)] + kappa[IX(i - 1, j)]) * 0.5f;
            // u 面的顏色梯度 (直接用相鄰兩格相減)
            float grad_phi = phi[0][IX(i, j)] - phi[0][IX(i - 1, j)];

            // F = sigma * kappa * grad(phi)
            u[IX_u(i, j)] += dt * Sigma * k * grad_phi;
        }
    }

#pragma omp parallel for
    for (int j = 2; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            // v 面的曲率取上下格子的平均
            float k = (kappa[IX(i, j)] + kappa[IX(i, j - 1)]) * 0.5f;
            float grad_phi = phi[0][IX(i, j)] - phi[0][IX(i, j - 1)];

            v[IX_v(i, j)] += dt * Sigma * k * grad_phi;
        }
    }
}

void MACSimulator::markFluidCells() {
    std::fill(cell_type.begin(), cell_type.end(), 0);

#pragma omp parallel for
    for (const auto &p : particles) {
        int i = (int)p.x;
        int j = (int)p.y;
        if (i >= 0 && i < nx && j >= 0 && j < ny) {
#pragma omp atomic write
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

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = IX(i, j);

            bool is_inter = 1;
            if (i < nx - 1)
                is_inter &= cell_type[IX(i + 1, j)];
            if (i > 0)
                is_inter &= cell_type[IX(i - 1, j)];
            if (j < ny - 1)
                is_inter &= cell_type[IX(i, j + 1)];
            if (j > 0)
                is_inter &= cell_type[IX(i, j - 1)];

            if (is_inter && particles_count[idx] < min_ppc) {
                int needed = target_ppc - particles_count[idx];

                for (int k = 0; k < needed; k++) {
                    float rx = i + ((float)rand() / RAND_MAX);
                    float ry = j + ((float)rand() / RAND_MAX);

                    float new_u = bilerp(u, nx + 1, ny, rx, ry - 0.5f);
                    float new_v = bilerp(v, nx, ny + 1, rx - 0.5f, ry);

                    particles.push_back({rx, ry, new_u, new_v});
                }
                cell_type[idx] = 1;
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
