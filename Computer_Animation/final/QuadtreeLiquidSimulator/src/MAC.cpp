#include "MAC.h"
#include <vector>

// PUBLIC
MACSimulator::MACSimulator(int width, int height) : nx(width), ny(height) {
    u.resize((nx + 1) * ny, 0.0);
    v.resize(nx * (ny + 1), 0.0);
    p.resize(nx * ny, 0.0);
    dye.resize(nx * ny, 0.0);

    initEigen();
}

void MACSimulator::update(float dt) {
    advect(dt);
    setBoundaries();
    project();
    setBoundaries();
}

void MACSimulator::addForce(
    float x, float y, float dx, float dy, float radius
) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double dist2 = (i - x) * (i - x) + (j - y) * (j - y);
            if (dist2 < radius * radius) {
                double falloff = exp(-dist2 / (radius * 0.5));
                u[IX_u(i, j)] += dx * falloff * 0.1;
                v[IX_v(i, j)] += dy * falloff * 0.1;
                dye[IX(i, j)] += 5.0 * falloff; // 加一點染料
            }
        }
    }
}

// PRIVATE
void MACSimulator::initEigen() {
    int N = nx * ny;
    A.resize(N, N);
    std::vector<Eigen::Triplet<double>> triplets;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = IX(i, j);

            int count = 0;
            if (i > 0) {
                triplets.push_back({idx, IX(i - 1, j), -1.0});
                count++;
            }
            if (i < nx - 1) {
                triplets.push_back({idx, IX(i + 1, j), -1.0});
                count++;
            }
            if (j > 0) {
                triplets.push_back({idx, IX(i, j - 1), -1.0});
                count++;
            }
            if (j < ny - 1) {
                triplets.push_back({idx, IX(i, j + 1), -1.0});
                count++;
            }

            double diag = (double)count;
            if (idx == 0) {
                diag += 1.0;
            }
            triplets.push_back({idx, idx, diag});
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    solver.compute(A);
}

void MACSimulator::setBoundaries() {
    for (int j = 0; j < ny; j++) {
        u[IX_u(0, j)] = 0.0;
        u[IX_u(nx, j)] = 0.0;
    }
    for (int i = 0; i < nx; i++) {
        v[IX_v(i, 0)] = 0.0;
        v[IX_v(i, ny)] = 0.0;
    }
}

void MACSimulator::advect(float dt) {
    std::vector<double> next_u(u.size()), next_v(v.size()),
        next_dye(dye.size());

    // 1. Advect Dye (位於網格中心)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double u_avg = (u[IX_u(i, j)] + u[IX_u(i + 1, j)]) * 0.5;
            double v_avg = (v[IX_v(i, j)] + v[IX_v(i, j + 1)]) * 0.5;
            double x_prev = i - u_avg * dt;
            double y_prev = j - v_avg * dt;
            next_dye[IX(i, j)] = bilerp(dye, nx, ny, x_prev, y_prev) * 0.995;
        }
    }

    // 2. Advect U (位於網格的垂直面上，X 座標有 -0.5 的偏移)
    for (int j = 0; j < ny; j++) {
        for (int i = 1; i < nx; i++) {
            double u_vel = u[IX_u(i, j)]; // u 自身的速度
            // 取周圍 4 個 v 面的平均來求出 u 面上的 v 速度
            double v_vel = (v[IX_v(i - 1, j)] +
                            v[IX_v(i, j)] +
                            v[IX_v(i - 1, j + 1)] +
                            v[IX_v(i, j + 1)]) *
                           0.25;

            double x_prev = (i - 0.5) - u_vel * dt;
            double y_prev = j - v_vel * dt;
            next_u[IX_u(i, j)] = bilerp(u, nx + 1, ny, x_prev + 0.5, y_prev);
        }
    }

    // 3. Advect V (位於網格的水平面上，Y 座標有 -0.5 的偏移)
    for (int j = 1; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            // 取周圍 4 個 u 面的平均來求出 v 面上的 u 速度
            double u_vel = (u[IX_u(i, j - 1)] +
                            u[IX_u(i + 1, j - 1)] +
                            u[IX_u(i, j)] +
                            u[IX_u(i + 1, j)]) *
                           0.25;
            double v_vel = v[IX_v(i, j)]; // v 自身的速度

            double x_prev = i - u_vel * dt;
            double y_prev = (j - 0.5) - v_vel * dt;
            next_v[IX_v(i, j)] = bilerp(v, nx, ny + 1, x_prev, y_prev + 0.5);
        }
    }
    u = std::move(next_u);
    v = std::move(next_v);
    dye = std::move(next_dye);
}

void MACSimulator::project() {
    int N = nx * ny;
    Eigen::VectorXd div(N);

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double d = u[IX_u(i + 1, j)] -
                       u[IX_u(i, j)] +
                       v[IX_v(i, j + 1)] -
                       v[IX_v(i, j)];
            div[IX(i, j)] = -d;
        }
    }

    Eigen::VectorXd pressure = solver.solve(div);
    for (int i = 0; i < N; i++)
        p[i] = pressure[i];

    for (int j = 0; j < ny; j++) {
        for (int i = 1; i < nx; i++) {
            u[IX_u(i, j)] -= (p[IX(i, j)] - p[IX(i - 1, j)]);
        }
    }
    for (int j = 1; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            v[IX_v(i, j)] -= (p[IX(i, j)] - p[IX(i, j - 1)]);
        }
    }
}

double MACSimulator::bilerp(
    const std::vector<double> &field, int w, int h, double x, double y
) const {
    x = std::max(0.0, std::min((double)w - 1.001, x));
    y = std::max(0.0, std::min((double)h - 1.001, y));

    int i = (int)x;
    int j = (int)y;
    double fx = x - i;
    double fy = y - j;

    double c00 = field[i + j * w];
    double c10 = field[(i + 1) + j * w];
    double c01 = field[i + (j + 1) * w];
    double c11 = field[(i + 1) + (j + 1) * w];

    return (c00 * (1 - fx) + c10 * fx) * (1 - fy) +
           (c01 * (1 - fx) + c11 * fx) * fy;
}
