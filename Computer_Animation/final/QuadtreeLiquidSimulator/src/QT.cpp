#include "QT.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>

// PUBLIC
QTSimulator::QTSimulator(int width, int height) : nx(width), ny(height) {
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

    root = new QuadtreeNode{(float)nx / 2, (float)ny / 2, (float)nx, 0};
    initQuadtree(root, 3);
}
void QTSimulator::reset() {
    particles.clear();
    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(u_old.begin(), u_old.end(), 0.0);
    std::fill(v_old.begin(), v_old.end(), 0.0);
    std::fill(p.begin(), p.end(), 0.0);
    std::fill(cell_type.begin(), cell_type.end(), 0);

    recursiveFree(root);
    root = new QuadtreeNode{(float)nx / 2, (float)ny / 2, (float)nx, 0};
    initQuadtree(root, 3);
}

void QTSimulator::update(float dt) {
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

    buildNewTree(dt);

    advectParticles(dt);
    markFluidCells();
    resampleParticles();
}

void QTSimulator::addWater(float x, float y, float radius) {
    float r2 = radius * radius;
    for (float px = x - radius; px < x + radius; px += 0.5) {
        for (float py = y - radius; py < y + radius; py += 0.5) {
            if (((px - x) * (px - x) + (py - y) * (py - y)) < r2)
                particles.push_back({px, py, 0.f, 0.f});
        }
    }

    recursiveUpdatePhi(root, x, y, radius);
}
void QTSimulator::delWater(float x, float y, float radius) {
    float r2 = radius * radius;

    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](auto &p) {
            return ((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y)) < r2;
        });

    particles.erase(new_end, particles.end());
    markFluidCells();
}

void QTSimulator::recursiveGetLines(
    QuadtreeNode *node, std::vector<Line> &lines
) const {
    int x1 = node->x - node->size / 2.f;
    int y1 = node->y - node->size / 2.f;
    int x2 = node->x + node->size / 2.f;
    int y2 = node->y + node->size / 2.f;

    lines.push_back({x1, y1, x2, y1});
    lines.push_back({x1, y1, x1, y2});
    lines.push_back({x1, y2, x2, y2});
    lines.push_back({x2, y1, x2, y2});

    if (node->is_leaf)
        return;

    for (int i = 0; i < 4; i++) {
        recursiveGetLines(node->children[i], lines);
    }
}
std::vector<Line> QTSimulator::getLines() const {
    std::vector<Line> lines;

    recursiveGetLines(root, lines);

    return lines;
}

void QTSimulator::setBoundaries(
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

void QTSimulator::particleToGrid() {

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

void QTSimulator::velExtrapolation() {
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

void QTSimulator::gridToParticle() {
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

void QTSimulator::advectParticles(float dt) {
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

void QTSimulator::project() {
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

void QTSimulator::applyGravity(float dt) {
    for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx; i++) {
            if ((j < ny && cell_type[IX(i, j)] == 1) ||
                (j > 0 && cell_type[IX(i, j - 1)] == 1)) {
                v[IX_v(i, j)] += G * dt;
            }
        }
    }
}

void QTSimulator::applySurfaceTension(float dt) {
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

void QTSimulator::markFluidCells() {
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

void QTSimulator::resampleParticles() {
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

float QTSimulator::bilerp(
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

void QTSimulator::bidistri(
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

void QTSimulator::initQuadtree(QuadtreeNode *node, int max_depth) {
    if (node->depth >= max_depth)
        return;

    subdivideNode(node);
    for (int i = 0; i < 4; i++) {
        initQuadtree(node->children[i], max_depth);
    }
}

void QTSimulator::subdivideNode(QuadtreeNode *node) {
    node->is_leaf = false;

    float child_size = node->size / 2.f;
    for (int i = 0; i < 4; i++) {
        float x = node->x + (i & 1 ? 1.f : -1.f) * child_size / 2.f;
        float y = node->y + (i & 2 ? 1.f : -1.f) * child_size / 2.f;
        node->children[i] = new QuadtreeNode{x, y, child_size, node->depth + 1};
    }
}

void QTSimulator::recursiveUpdatePhi(
    QuadtreeNode *node, float cx, float cy, float radius
) {
    float new_phi = circleSDF(cx, cy, radius, node->x, node->y);
    node->phi = std::min(node->phi, new_phi);

    if (!node->is_leaf) {
        for (int i = 0; i < 4; i++)
            recursiveUpdatePhi(node->children[i], cx, cy, radius);
        return;
    }

    if (node->size > 1.f && std::abs(new_phi) < node->size) {
        subdivideNode(node);
        for (int i = 0; i < 4; i++)
            recursiveUpdatePhi(node->children[i], cx, cy, radius);
    }
}

void QTSimulator::recursiveFree(QuadtreeNode *node) {
    if (node->is_leaf) {
        delete node;
        return;
    }

    for (int i = 0; i < 4; i++) {
        recursiveFree(node->children[i]);
    }
    delete node;
}

QuadtreeNode *QTSimulator::getNodeAt(QuadtreeNode *node, float x, float y) {
    if (node->is_leaf)
        return node;

    int child_idx = 0;
    if (x > node->x)
        child_idx |= 1;
    if (y > node->y)
        child_idx |= 2;

    return getNodeAt(node->children[child_idx], x, y);
}

void QTSimulator::getNodesIn(
    QuadtreeNode *node,
    float x,
    float y,
    float radius,
    std::vector<QuadtreeNode *> &nodes
) {
    if (!node)
        return;

    float half_size = node->size / 2.f;
    if (node->x + half_size < x - radius ||
        node->x - half_size > x + radius ||
        node->y + half_size < y - radius ||
        node->y - half_size > y + radius) {
        return;
    }

    if (node->is_leaf) {
        nodes.push_back(node);
    } else {
        for (int i = 0; i < 4; i++) {
            getNodesIn(node->children[i], x, y, radius, nodes);
        }
    }
}

void QTSimulator::advectQuadtreePhi(QuadtreeNode *node, float dt) {
    if (!node->is_leaf) {
        for (int i = 0; i < 4; i++)
            advectQuadtreePhi(node->children[i], dt);
        return;
    }

    float u_val = bilerp(u, nx + 1, ny, node->x, node->y - 0.5f);
    float v_val = bilerp(v, nx, ny + 1, node->x - 0.5f, node->y);

    float past_x = std::clamp(node->x - u_val * dt, 0.0f, (float)nx);
    float past_y = std::clamp(node->y - v_val * dt, 0.0f, (float)ny);

    node->phi = MLSinterpolate(past_x, past_y).phi;
}

void QTSimulator::computeSizingFunction(QuadtreeNode *node) {
    if (!node->is_leaf) {
        for (int i = 0; i < 4; i++)
            computeSizingFunction(node->children[i]);
        return;
    }

    if (std::abs(node->phi) < node->size * 1.5f) {
        float gamma_phi = 4.f;
        float gamma_u = 3.f;

        float geom_term = gamma_phi / (std::abs(node->phi) + 0.1f);

        int i = std::clamp((int)node->x, 1, nx - 2);
        int j = std::clamp((int)node->y, 1, ny - 2);

        float du_dx = std::abs(u[IX_u(i + 1, j)] - u[IX_u(i, j)]);
        float dv_dy = std::abs(v[IX_v(i, j + 1)] - v[IX_v(i, j)]);
        float vel_term = gamma_u * std::sqrt(du_dx * du_dx + dv_dy * dv_dy);

        node->S = geom_term + vel_term;
    } else {
        node->S = 0.0f;
    }
}

void QTSimulator::propagateSizingFunction() {
    std::vector<QuadtreeNode *> leaves;
    collectLeafNodes(root, leaves);

    std::vector<float> new_S(leaves.size(), 0.f);

    for (int iter = 0; iter < 5; iter++) {
        for (int i = 0; i < leaves.size(); i++) {
            auto node = leaves[i];
            float max_S = node->S;

            float step = node->size;
            max_S =
                std::max(max_S, getNodeAt(root, node->x + step, node->y)->S);
            max_S =
                std::max(max_S, getNodeAt(root, node->x - step, node->y)->S);
            max_S =
                std::max(max_S, getNodeAt(root, node->x, node->y + step)->S);
            max_S =
                std::max(max_S, getNodeAt(root, node->x, node->y - step)->S);

            new_S[i] = max_S;
        }

        for (int i = 0; i < leaves.size(); i++)
            leaves[i]->S = new_S[i];
    }
}

void QTSimulator::collectLeafNodes(
    QuadtreeNode *node, std::vector<QuadtreeNode *> &leaves
) {
    if (node->is_leaf) {
        leaves.push_back(node);
        return;
    }
    for (int i = 0; i < 4; i++) {
        collectLeafNodes(node->children[i], leaves);
    }
}

InterpolatedData QTSimulator::MLSinterpolate(float x, float y) {

    InterpolatedData result = {0.0f, 1000.f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    QuadtreeNode *target = getNodeAt(root, x, y);
    if (!target)
        return result;

    float search_radius = target->size * 1.5f;
    std::vector<QuadtreeNode *> neighbors;
    getNodesIn(root, x, y, search_radius, neighbors);

    Eigen::Matrix3f A = Eigen::Matrix3f::Zero();
    Eigen::MatrixXf b(3, 7);
    b.setZero();

    for (auto n : neighbors) {
        float dx = std::abs(n->x - x);
        float dy = std::abs(n->y - y);

        float h = n->size;

        float wx = std::max(1.f - (dx / h), 0.0f);
        float wy = std::max(1.f - (dy / h), 0.0f);
        float weight = wx * wy;

        if (weight <= 0.0f)
            continue;

        Eigen::Vector3f z_i(n->x, n->y, 1.f);

        A += weight * (z_i * z_i.transpose());
        b.col(0) += weight * z_i * n->S;
        b.col(1) += weight * z_i * n->phi;
        b.col(2) += weight * z_i * n->pressure;
        b.col(3) += weight * z_i * n->ul;
        b.col(4) += weight * z_i * n->ur;
        b.col(5) += weight * z_i * n->vl;
        b.col(6) += weight * z_i * n->vr;
    }

    A += Eigen::Matrix3f::Identity() * 1e-6;

    Eigen::MatrixXf c = A.ldlt().solve(b);
    Eigen::Vector3f query(x, y, 1.f);
    Eigen::VectorXf datas = c.transpose() * query;

    result.S = datas[0];
    result.phi = datas[1];
    result.pressure = datas[2];
    result.ul = datas[3];
    result.ur = datas[4];
    result.vl = datas[5];
    result.vr = datas[6];

    return result;
}

void QTSimulator::recursiveBuildTree(QuadtreeNode *node) {
    if (node->size < 1.5f)
        return;

    auto data = MLSinterpolate(node->x, node->y);
    float phi_p = data.phi;
    float S_p = data.S;

    float max_vel = std::max(max_u, max_v);
    float safe_band = node->size * 1.5f + max_vel * 0.016f;

    if (std::abs(phi_p) < safe_band && S_p > (1.f / node->size)) {
        subdivideNode(node);

        for (int i = 0; i < 4; i++) {
            recursiveBuildTree(node->children[i]);
        }
    }
}
void QTSimulator::buildNewTree(float dt) {
    auto new_root =
        new QuadtreeNode{(float)nx / 2, (float)ny / 2, (float)nx, 0};

    computeSizingFunction(root);
    propagateSizingFunction();

    recursiveBuildTree(new_root);

    advectQuadtreePhi(new_root, dt);

    recursiveFree(root);
    root = new_root;
}
