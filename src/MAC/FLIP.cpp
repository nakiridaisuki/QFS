#include "MAC/FLIP.h"
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <vector>

// PUBLIC
MACFLIP::MACFLIP(int width, int height) : MACSimulatorBase(width, height) {
    simulator_type = SimType::MAC_FLIP;
    particles_count.resize(nx * ny, 0);
}

void MACFLIP::update(float dt) {
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

void MACFLIP::addWater(float x, float y, float radius) {
    float r2 = radius * radius;
    for (float px = x - radius; px < x + radius; px += 0.5) {
        for (float py = y - radius; py < y + radius; py += 0.5) {
            if (((px - x) * (px - x) + (py - y) * (py - y)) < r2)
                particles.push_back({px, py, 0.f, 0.f});
        }
    }
}
void MACFLIP::delWater(float x, float y, float radius) {
    float r2 = radius * radius;

    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](auto &p) {
            return ((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y)) < r2;
        });

    particles.erase(new_end, particles.end());
    markFluidCells();
}

void MACFLIP::particleToGrid() {

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

void MACFLIP::gridToParticle() {
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

void MACFLIP::advectParticles(float dt) {
    float local_max_u = 0.0f;
    float local_max_v = 0.0f;

#pragma omp parallel for reduction(max : local_max_u)                          \
    reduction(max : local_max_v)
    for (auto &p : particles) {
        p.x += p.u * dt;
        p.y += p.v * dt;

        p.x = std::max(0.001f, std::min(float(nx) - 0.001f, p.x));
        p.y = std::max(0.001f, std::min(float(ny) - 0.001f, p.y));
    }
}

void MACFLIP::applySurfaceTension(float dt) {
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
            float div_x     = (nx_n[IX(i + 1, j)] - nx_n[IX(i - 1, j)]) * 0.5f;
            float div_y     = (ny_n[IX(i, j + 1)] - ny_n[IX(i, j - 1)]) * 0.5f;
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
            float k        = (kappa[IX(i, j)] + kappa[IX(i, j - 1)]) * 0.5f;
            float grad_phi = phi[0][IX(i, j)] - phi[0][IX(i, j - 1)];

            v[IX_v(i, j)] += dt * Sigma * k * grad_phi;
        }
    }
}

void MACFLIP::markFluidCells() {
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

void MACFLIP::resampleParticles() {
    const int target_ppc = 4;
    const int max_ppc    = 8;
    const int min_ppc    = 3;

    std::fill(particles_count.begin(), particles_count.end(), 0);
    for (const auto &p : particles) {
        int i = std::clamp((int)p.x, 0, nx - 1);
        int j = std::clamp((int)p.y, 0, ny - 1);
        particles_count[IX(i, j)]++;
    }

    std::fill(current_count.begin(), current_count.end(), 0);
    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](const auto &p) {
            int i   = std::clamp((int)p.x, 0, nx - 1);
            int j   = std::clamp((int)p.y, 0, ny - 1);
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
