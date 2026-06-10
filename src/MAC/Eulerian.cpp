#include "MAC/Eulerian.h"
#include <algorithm>
#include <math.h>
#include <vector>

// PUBLIC
MACEulerian::MACEulerian(int width, int height)
    : MACSimulatorBase(width, height) {
    simulator_type = SimType::MAC_Eulerian;
    density.resize(nx * ny, 0);
    density_old.resize(nx * ny, 0);
}

void MACEulerian::update(float dt) {
    u_old       = u;
    v_old       = v;
    density_old = density;
    setBoundaries(u_old, v_old);

    advectDatas(dt);
    markFluidCells();

    applyGravity(dt);
    applySurfaceTension(dt);

    setBoundaries(u, v);
    project();
    setBoundaries(u, v);

    velExtrapolation();
}

void MACEulerian::addWater(float x, float y, float radius) {
    float r2 = radius * radius;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            float cx = i + 0.5f;
            float cy = j + 0.5f;
            if (distance2(cx, cy, x, y) < r2) {
                density[IX(i, j)]   = 1.f;
                cell_type[IX(i, j)] = 1;
            }
        }
    }
}
void MACEulerian::delWater(float x, float y, float radius) {
    float r2 = radius * radius;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            float cx = i + 0.5f;
            float cy = j + 0.5f;
            if (distance2(cx, cy, x, y) < r2) {
                density[IX(i, j)]   = 0.f;
                cell_type[IX(i, j)] = 0;
            }
        }
    }
}

void MACEulerian::applySurfaceTension(float dt) {
    if (Sigma <= 0.0f)
        return;

#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < nx * ny; i++) {
            phi[0][i] = density[i];
            phi[1][i] = 0.f;
            nx_n[i]   = 0.f;
            ny_n[i]   = 0.f;
            kappa[i]  = 0.f;
        }

        for (int iter = 0; iter < 2; iter++) {
#pragma omp for
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

#pragma omp for
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                // 利用中心差分計算梯度
                float dx = (phi[0][IX(i + 1, j)] - phi[0][IX(i - 1, j)]) * 0.5f;
                float dy = (phi[0][IX(i, j + 1)] - phi[0][IX(i, j - 1)]) * 0.5f;
                float len = std::sqrt(
                    dx * dx + dy * dy + 1e-8f
                ); // 加上 1e-8 避免除以零
                nx_n[IX(i, j)] = dx / len;
                ny_n[IX(i, j)] = dy / len;
            }
        }

#pragma omp for
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                float div_x = (nx_n[IX(i + 1, j)] - nx_n[IX(i - 1, j)]) * 0.5f;
                float div_y = (ny_n[IX(i, j + 1)] - ny_n[IX(i, j - 1)]) * 0.5f;
                kappa[IX(i, j)] = -(div_x + div_y);
            }
        }

#pragma omp for
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

#pragma omp for
        for (int j = 2; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                // v 面的曲率取上下格子的平均
                float k        = (kappa[IX(i, j)] + kappa[IX(i, j - 1)]) * 0.5f;
                float grad_phi = phi[0][IX(i, j)] - phi[0][IX(i, j - 1)];

                v[IX_v(i, j)] += dt * Sigma * k * grad_phi;
            }
        }
    }
}

void MACEulerian::advectDatas(float dt) {
#pragma omp parallel
    {
#pragma omp for nowait
        // Advect U velocity
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i <= nx; i++) {
                float x = (float)i;
                float y = (float)j + 0.5f;

                float u_vel = u_old[IX_u(i, j)];
                float v_vel = bilerp(v_old, nx, ny + 1, x - 0.5f, y);

                float x_prev = x - u_vel * dt;
                float y_prev = y - v_vel * dt;

                u[IX_u(i, j)] =
                    bilerp(u_old, nx + 1, ny, x_prev, y_prev - 0.5f);
            }
        }

#pragma omp for nowait
        // Advect V velocity
        for (int j = 0; j <= ny; j++) {
            for (int i = 0; i < nx; i++) {
                float x = (float)i + 0.5f;
                float y = (float)j;

                float u_vel = bilerp(u_old, nx + 1, ny, x, y - 0.5f);
                float v_vel = v_old[IX_v(i, j)];

                float x_prev = x - u_vel * dt;
                float y_prev = y - v_vel * dt;

                v[IX_v(i, j)] =
                    bilerp(v_old, nx, ny + 1, x_prev - 0.5f, y_prev);
            }
        }

#pragma omp for nowait
        // Advect density
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                float x = (float)i + 0.5f;
                float y = (float)j + 0.5f;

                float u_vel = bilerp(u_old, nx + 1, ny, x, y - 0.5f);
                float v_vel = bilerp(v_old, nx, ny + 1, x - 0.5f, y);

                float x_prev = x - u_vel * dt;
                float y_prev = y - v_vel * dt;

                density[IX(i, j)] =
                    bilerp(density_old, nx, ny, x_prev - 0.5f, y_prev - 0.5f);
            }
        }
    }
}

void MACEulerian::markFluidCells() {
    std::fill(cell_type.begin(), cell_type.end(), 0);

#pragma omp parallel for
    for (int i = 0; i < nx * ny; i++) {
        if (density[i] > 0.1f)
            cell_type[i] = 1;
        else
            cell_type[i] = 0;
    }
}
