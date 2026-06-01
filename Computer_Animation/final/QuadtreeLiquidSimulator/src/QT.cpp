#include "QT.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

// PUBLIC
QTSimulator::QTSimulator(int width, int height) : nx(width), ny(height) {
    cell_type.resize(nx * ny, 0);
    particle_idx.resize(nx * ny, std::vector<int>(8));
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
    std::fill(cell_type.begin(), cell_type.end(), 0);

    recursiveFree(root);
    root = new QuadtreeNode{(float)nx / 2, (float)ny / 2, (float)nx, 0};
    initQuadtree(root, 3);
    cacheLeaves(root);
    cacheNeighbors(root);
    findAllEdges();
    QTu = QTu_new;
    QTv = QTv_new;
}

void QTSimulator::update(float dt) {
    max_u = 0.0f;
    max_v = 0.0f;

    // Generate new Quad Tree
    computeSizingFunction(dt);
    commitQuadtreeS();
    propagateSizingFunction();

    auto new_root =
        new QuadtreeNode{(float)nx / 2, (float)ny / 2, (float)nx, 0};
    recursiveBuildTree(new_root, dt);

    smoothing(new_root);
    cacheLeaves(new_root);
    cacheNeighbors(new_root);

    findAllEdges();

    advectQuadtreeDatas(dt);
    recursiveFree(root);
    root = new_root;

    QTu = QTu_new;
    QTv = QTv_new;

    particleToGrid();
    // particleSurfaceToGrid();

    QTapplyGravity(dt);
    for (auto &face : QTu)
        face.val_old = face.val;
    for (auto &face : QTv)
        face.val_old = face.val;
    QTproject();

    gridToParticle();
    advectParticles(dt);
    resampleParticles();

    // // Redistance
    updateParticleIdx();
    redistancing();
}

void QTSimulator::particleSurfaceToGrid() {
    for (auto leaf : cached_leaves) {
        if (leaf->size < 1.5f) {
            auto nearest_p = nearestParticle(leaf->x, leaf->y, 2.f);
            if (nearest_p.x != 1e6) {
                float dist =
                    std::sqrt(
                        (nearest_p.x - leaf->x) * (nearest_p.x - leaf->x) +
                        ((nearest_p.y - leaf->y) * (nearest_p.y - leaf->y))
                    ) -
                    particle_radius;
                leaf->phi = std::min(leaf->phi, dist);
            }
        }
    }
}

void QTSimulator::particleToGrid() {
    float search_radius = 1.5f;
    float weight_eps = 1e-6;

    auto updateFace = [&](QuadtreeEdge &face, bool is_u) {
        if (face.length > 1.5f)
            return;

        float sum_w = 0;
        float sum_val = 0;

        int ix = std::clamp((int)face.x, 0, nx - 1);
        int iy = std::clamp((int)face.y, 0, ny - 1);

        for (int dj = -2; dj <= 2; dj++) {
            for (int di = -2; di <= 2; di++) {
                int ni = ix + di, nj = iy + dj;

                if (ni < 0 || ni >= nx || nj < 0 || nj >= ny)
                    continue;

                for (int idx : particle_idx[IX(ni, nj)]) {
                    auto &p = particles[idx];

                    float dx = std::abs(p.x - face.x);
                    float dy = std::abs(p.y - face.y);

                    float wx = std::max(0.f, 1.f - dx / face.length);
                    float wy = std::max(0.f, 1.f - dy / face.length);
                    float w = wx * wy;

                    if (w > 0) {
                        sum_w += w;
                        if (is_u)
                            sum_val += w * p.u;
                        else
                            sum_val += w * p.v;
                    }
                }
            }
        }
        if (sum_w > weight_eps)
            face.val = sum_val / sum_w;
    };

    for (int i = 0; i < QTu.size(); i++)
        updateFace(QTu[i], true);

    for (int i = 0; i < QTv.size(); i++)
        updateFace(QTv[i], false);
}

void QTSimulator::updateParticleIdx() {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            particle_idx[IX(i, j)].clear();
        }
    }
    for (int i = 0; i < particles.size(); i++) {
        auto &p = particles[i];
        int ix = p.x;
        int iy = p.y;

        ix = std::clamp(ix, 0, nx - 1);
        iy = std::clamp(iy, 0, ny - 1);

        particle_idx[IX(ix, iy)].push_back(i);
    }
}

void QTSimulator::addWater(float x, float y, float radius) {
    recursiveUpdatePhi(root, x, y, radius);
    smoothing(root);
    cacheLeaves(root);
    cacheNeighbors(root);
}
void QTSimulator::delWater(float x, float y, float radius) {
    // float r2 = radius * radius;
    //
    // auto new_end =
    //     std::remove_if(particles.begin(), particles.end(), [&](auto &p) {
    //         return ((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y)) < r2;
    //     });
    //
    // particles.erase(new_end, particles.end());
    // markFluidCells();
}

void QTSimulator::recursiveGetLines(
    QuadtreeNode *node, std::vector<Line> &lines
) const {
    float x1 = node->x - node->size / 2.f;
    float y1 = node->y - node->size / 2.f;
    float x2 = node->x + node->size / 2.f;
    float y2 = node->y + node->size / 2.f;

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

float compute_Wk(QuadtreeEdge &face) {
    float total = 0, liquid = 0;
    for (int i = 0; i < face.adj_cells.size(); i++) {
        total += face.grad_coeff[i] * face.adj_cells[i]->phi;
        if (face.adj_cells[i]->phi <= 0.f)
            liquid += face.grad_coeff[i] * face.adj_cells[i]->phi;
    }
    if (std::abs(liquid) < 1e-6f)
        return 1.f;
    return total / liquid;
}
void QTSimulator::QTproject() {

    int fluid_count = 0;
    for (auto leaf : cached_leaves) {
        if (leaf->phi <= 0.0)
            leaf->fluid_id = fluid_count++;
        else
            leaf->fluid_id = -1;
    }

    if (fluid_count == 0)
        return;

    Eigen::VectorXf div(fluid_count);
    div.setZero();
    std::vector<Eigen::Triplet<float>> QTtriplets;

    auto addFace = [&](QuadtreeEdge &face) {
        float u_star = face.val;
        float VA = face.length * (1.f - face.solid_fraction);

        bool has_fluid = false;
        for (auto cell : face.adj_cells)
            has_fluid |= (cell->fluid_id != -1);
        if (!has_fluid)
            return;

        float W_k = compute_Wk(face);
        float W_prime = std::max(W_k, 0.01f);
        float face_weight = VA * W_prime;

        for (int i = 0; i < face.adj_cells.size(); i++) {
            auto cell_i = face.adj_cells[i];
            if (cell_i->fluid_id == -1)
                continue;
            int row = cell_i->fluid_id;

            div[row] += VA * u_star * face.grad_coeff[i];

            for (int j = 0; j < face.adj_cells.size(); j++) {
                auto cell_j = face.adj_cells[j];
                if (cell_j->fluid_id == -1)
                    continue;
                int col = cell_j->fluid_id;

                float val =
                    face_weight * face.grad_coeff[i] * face.grad_coeff[j];
                QTtriplets.push_back({row, col, val});
            }
        }
    };
    for (auto &face : QTu)
        addFace(face);
    for (auto &face : QTv)
        addFace(face);

    // Eigen solvers
    Eigen::SparseMatrix<float> A(fluid_count, fluid_count);
    A.setFromTriplets(QTtriplets.begin(), QTtriplets.end());
    solver.compute(A);
    Eigen::VectorXf pressure = solver.solve(div);

    auto updateFace = [&](QuadtreeEdge &face) {
        bool has_fluid = false;
        for (auto cell : face.adj_cells)
            has_fluid |= (cell->fluid_id != -1);
        if (!has_fluid)
            return;

        float grad_p = 0;

        for (int i = 0; i < face.adj_cells.size(); i++) {
            auto cell = face.adj_cells[i];
            float p_val =
                (cell->fluid_id != -1 ? pressure[cell->fluid_id] : 0.0f);
            grad_p += face.grad_coeff[i] * p_val;
        }

        float W_prime = std::max(compute_Wk(face), 0.01f);
        face.val -= W_prime * grad_p;
    };

    for (auto &face : QTu)
        updateFace(face);
    for (auto &face : QTv)
        updateFace(face);
}

void QTSimulator::advectParticles(float dt) {
    for (auto &p : particles) {
        auto d1 = MLSinterpolate(p.x, p.y);
        float mid_x = p.x + d1.u * dt * 0.5f;
        float mid_y = p.y + d1.v * dt * 0.5f;

        auto d2 = MLSinterpolate(mid_x, mid_y);
        p.x += d2.u * dt;
        p.y += d2.v * dt;

        float eps = 1e-3;
        p.x = std::clamp(p.x, eps, (float)nx - eps);
        p.y = std::clamp(p.y, eps, (float)ny - eps);
    }
}

void QTSimulator::QTapplyGravity(float dt) {
    for (auto &face : QTv) {
        bool is_water = false;
        for (auto node : face.adj_cells)
            is_water |= (node->phi < 0.0f);

        if (is_water)
            face.val += G * dt;
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

// Quad Tree
void QTSimulator::initQuadtree(QuadtreeNode *node, int max_depth) {
    if (node->depth >= max_depth)
        return;

    subdivideNode(node);
    for (int i = 0; i < 4; i++) {
        initQuadtree(node->children[i], max_depth);
    }
}

// NeighborData QTSimulator::getNeighborData(QuadtreeNode *node, int direction)
// {
//     // 0: left
//     // 1: up
//     // 2: right
//     // 3: down
//
//     NeighborData data;
//     float h_C = node->size;
//
//     std::vector<std::pair<int, QuadtreeNode *>> neighbors;
//     getNeighbors(node, neighbors);
//     std::vector<QuadtreeNode *> tmp;
//
//     for (auto [dir, node] : neighbors)
//         if (dir == direction)
//             tmp.push_back(node);
//
//     if (tmp.empty()) {
//         data.datas = InterpolatedData(node);
//         data.distance = h_C;
//         return data;
//     }
//
//     if (tmp.size() == 1) {
//         // 情況 1 & 2：只有一個鄰居 (等大或較大)
//         QuadtreeNode *nb = tmp[0];
//         float h_nb = nb->size;
//
//         data.datas = InterpolatedData(nb);
//         data.distance = (h_C + h_nb) * 0.5f; // 中心點距離公式: (h_C + h_nb)
//         / 2
//     } else if (tmp.size() == 2) {
//         // 情況 3：T-junction，有兩個較小的鄰居
//         auto f0 = InterpolatedData(tmp[0]);
//         auto f1 = InterpolatedData(tmp[1]);
//
//         data.datas = (f0 + f1) * 0.5f;             // 數值取平均
//         data.distance = (h_C + h_C * 0.5f) * 0.5f; // 距離為 0.75 * h_C
//     }
//
//     return data;
// }

void QTSimulator::computeSizingFunction(float dt) {

    for (auto leaf : cached_leaves) {
        if (std::abs(leaf->phi) < leaf->size * 1.5f) {
            float gamma_phi = 4.f;
            float gamma_u = 3.f;

            float T1 = 0.9, T2 = 0.01;
            float R_t = std::pow(T1, dt / T2);

            // // 0: left
            // // 1: up
            // // 2: right
            // // 3: down
            // NeighborData L = getNeighborData(leaf, 0);
            // NeighborData U = getNeighborData(leaf, 1);
            // NeighborData R = getNeighborData(leaf, 2);
            // NeighborData D = getNeighborData(leaf, 3);
            //
            // float f_C = leaf->phi;
            // float d2f_dx2 =
            //     2.0f * (L.datas.phi / (L.distance * (L.distance +
            //     R.distance)) -
            //             f_C / (L.distance * R.distance) +
            //             R.datas.phi / (R.distance * (L.distance +
            //             R.distance)));
            // float d2f_dy2 =
            //     2.0f * (D.datas.phi / (D.distance * (D.distance +
            //     U.distance)) -
            //             f_C / (D.distance * U.distance) +
            //             U.datas.phi / (U.distance * (D.distance +
            //             U.distance)));
            //
            // float laplacian = std::abs(d2f_dx2 + d2f_dy2);
            // float geom_term = gamma_phi * laplacian;
            float geom_term = gamma_phi / (std::abs(leaf->phi) + 0.1f);

            float ul = getVelocity(QTu, leaf->ul_id);
            float ur = getVelocity(QTu, leaf->ur_id);
            float vl = getVelocity(QTv, leaf->vl_id);
            float vr = getVelocity(QTv, leaf->vr_id);

            float du_dx = (ul - ur) / leaf->size;
            float dv_dy = (vl - vr) / leaf->size;

            float vel_term = gamma_u * std::sqrt(du_dx * du_dx + dv_dy * dv_dy);

            float S = geom_term + vel_term;

            auto data = advect(leaf->x, leaf->y, dt);
            float S_star = data.S;

            leaf->S_new = std::max(R_t * S_star, S);
        } else {
            leaf->S_new = 0.0f;
        }
    }
}

void QTSimulator::commitQuadtreeS() {
    for (auto leaf : cached_leaves) {
        leaf->S = leaf->S_new;
    }
}

void QTSimulator::commitQuadtreePhi() {
    for (auto leaf : cached_leaves) {
        leaf->phi = leaf->phi_new;
    }
}

void QTSimulator::propagateSizingFunction() {
    int iterations = 5;
    for (int iter = 0; iter < iterations; iter++) {
        for (int i = 0; i < cached_leaves.size(); i++) {
            auto node = cached_leaves[i];

            float self_area = node->size * node->size;
            node->S_new = node->S * self_area;

            float total_area = self_area;
            for (auto n : node->cached_neighbors) {

                float area = n->size * n->size;
                node->S_new += std::max(n->S, node->S) * area;
                total_area += area;
            }

            node->S_new /= total_area;
        }

        for (int i = 0; i < cached_leaves.size(); i++)
            cached_leaves[i]->S = cached_leaves[i]->S_new;
    }
}
void QTSimulator::recursiveBuildTree(QuadtreeNode *node, float dt) {
    if (node->size < 1.5f)
        return;

    auto advected_data = advect(node->x, node->y, dt);
    auto data = MLSinterpolate(node->x, node->y);

    float exp_phi = advected_data.phi;
    float exp_S = data.S;

    if (std::abs(exp_phi) < node->size && exp_S > (1.f / node->size)) {
        subdivideNode(node);

        for (int i = 0; i < 4; i++) {
            recursiveBuildTree(node->children[i], dt);
        }
    }
}

InterpolatedData QTSimulator::advect(float x, float y, float dt) {
    auto current_data = MLSinterpolate(x, y);
    // float u_val = bilerp(u, nx + 1, ny, x, y);
    // float v_val = bilerp(v, nx, ny + 1, x, y);
    float u_val = current_data.u;
    float v_val = current_data.v;

    float past_x = std::clamp(x - u_val * dt, 0.0f, (float)nx);
    float past_y = std::clamp(y - v_val * dt, 0.0f, (float)ny);

    return MLSinterpolate(past_x, past_y);
}

void QTSimulator::advectQuadtreeDatas(float dt) {

    for (auto leaf : cached_leaves) {
        auto advected_leaf = advect(leaf->x, leaf->y, dt);
        leaf->phi = advected_leaf.phi;
        leaf->S = advected_leaf.S;
        leaf->pressure = advected_leaf.pressure;
    }

    for (auto &face : QTu_new) {
        auto advected_face = advect(face.x, face.y, dt);
        face.val = advected_face.u;
    }
    for (auto &face : QTv_new) {
        auto advected_face = advect(face.x, face.y, dt);
        face.val = advected_face.v;
    }
}
void QTSimulator::gridToParticle() {
    float flipRatio = 0.95f;

#pragma omp parallel for
    for (auto &p : particles) {
        auto data = MLSinterpolate(p.x, p.y);
        float u_new_grid = data.u;
        float v_new_grid = data.v;

        float u_old_grid = data.u_old;
        float v_old_grid = data.v_old;

        float u_pic = u_new_grid;
        float v_pic = v_new_grid;

        float u_flip = p.u + (u_new_grid - u_old_grid);
        float v_flip = p.v + (v_new_grid - v_old_grid);

        p.u = flipRatio * u_flip + (1.f - flipRatio) * u_pic;
        p.v = flipRatio * v_flip + (1.f - flipRatio) * v_pic;
    }
}
void QTSimulator::resampleParticles() {
    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](const auto &p) {
            QuadtreeNode *node = getNodeAt(root, p.x, p.y);

            if (!node)
                return true;

            if (node->size > 1.5f)
                return true;

            if (node->phi < -3.f * node->size)
                return true;

            node->has_particle = true;
            return false;
        });
    particles.erase(new_end, particles.end());

    for (auto leaf : cached_leaves) {
        if (leaf->phi < 0 &&
            leaf->size < 1.5f &&
            leaf->phi > -leaf->size &&
            !leaf->has_particle) {

            float px =
                leaf->x + (((rand() % 100) / 100.0f) - 0.5f) * leaf->size;
            float py =
                leaf->y + (((rand() % 100) / 100.0f) - 0.5f) * leaf->size;
            auto data = MLSinterpolate(px, py);
            particles.push_back({px, py, data.u, data.v});
        }
    }
}

void QTSimulator::redistancing() {
    std::priority_queue<
        QuadtreeNode *,
        std::vector<QuadtreeNode *>,
        QTNSurfaceFirst>
        pq;

    for (auto leaf : cached_leaves) {
        // if (std::abs(leaf->phi) < 0.5) {
        //     leaf->phi_new = leaf->phi;
        if (std::abs(leaf->phi) < leaf->size * 1.5f) {

            if (leaf->size < 1.5f) {
                auto nearest_p = nearestParticle(leaf->x, leaf->y, 2.f);

                if (nearest_p.x == 1e6 || nearest_p.y == 1e6) {
                    leaf->phi_new = std::numeric_limits<float>::infinity();
                    continue;
                }

                float dist = std::sqrt(
                    (leaf->x - nearest_p.x) * (leaf->x - nearest_p.x) +
                    (leaf->y - nearest_p.y) * (leaf->y - nearest_p.y)
                );
                dist -= particle_radius;

                leaf->phi = std::min(leaf->phi, dist);
            }
            leaf->phi_new = leaf->phi;
            pq.push(leaf);
        } else {
            leaf->phi_new = std::numeric_limits<float>::infinity();
        }
    }

    if (pq.empty()) {
        return;
    }

    while (!pq.empty()) {
        auto curr_node = pq.top();
        pq.pop();

        if (curr_node->known)
            continue;
        curr_node->known = true;

        for (auto neighbor : curr_node->cached_neighbors) {
            float dist = std::sqrt(
                pow(curr_node->x - neighbor->x, 2) +
                pow(curr_node->y - neighbor->y, 2)
            );

            dist = neighbor->phi < 0 ? -dist : dist;

            float proposed_phi_new = curr_node->phi_new + dist;

            if (std::abs(proposed_phi_new) < std::abs(neighbor->phi_new)) {
                neighbor->phi_new = proposed_phi_new;
                pq.push(neighbor);
            }
        }
    }
    commitQuadtreePhi();
}

Particle QTSimulator::nearestParticle(float x, float y, float radius) {
    int ix = std::clamp((int)x, 0, nx - 1);
    int iy = std::clamp((int)y, 0, ny - 1);

    Particle result{1e6, 1e6};
    float min_dis = 1e6;
    for (int j = -radius; j <= radius; j++) {
        for (int i = -radius; i <= radius; i++) {
            int ni = ix + i, nj = iy + j;

            if (ni < 0 || ni >= nx || nj < 0 || nj >= ny)
                continue;

            for (int idx : particle_idx[IX(ni, nj)]) {
                auto &p = particles[idx];

                float tmp_dis = (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y);
                if (tmp_dis < min_dis) {
                    min_dis = tmp_dis;
                    result = p;
                }
            }
        }
    }
    return result;
}

void QTSimulator::getNeighbors(
    QuadtreeNode *new_root,
    QuadtreeNode *node,
    std::vector<std::pair<int, QuadtreeNode *>> &neighbors
) {
    // 0: left
    // 1: up
    // 2: right
    // 3: down

    if (!node)
        return;

    float hs = node->size / 2.f;
    float qs = node->size / 4.f;

    float eps = 0.1;

    std::vector<std::pair<float, float>> sample_points = {
        // left
        {node->x - hs - eps, node->y + qs},
        {node->x - hs - eps, node->y - qs},

        // up
        {node->x + qs, node->y - hs - eps},
        {node->x - qs, node->y - hs - eps},

        // right
        {node->x + hs + eps, node->y + qs},
        {node->x + hs + eps, node->y - qs},

        // down
        {node->x + qs, node->y + hs + eps},
        {node->x - qs, node->y + hs + eps},
    };

    for (int i = 0; i < 8; i++) {
        auto [qx, qy] = sample_points[i];
        if (qx < 0.0f || qx >= (float)nx || qy < 0.0f || qy >= (float)ny) {
            continue;
        }

        QuadtreeNode *neighbor = getNodeAt(new_root, qx, qy);

        if (neighbor && neighbor != node) {
            neighbors.push_back({i / 2, neighbor});
        }
    }

    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(
        std::unique(neighbors.begin(), neighbors.end()), neighbors.end()
    );
}

void QTSimulator::subdivideNode(QuadtreeNode *node) {
    node->is_leaf = false;

    float child_size = node->size / 2.f;
    for (int i = 0; i < 4; i++) {
        float x = node->x + (i & 1 ? 1.f : -1.f) * child_size / 2.f;
        float y = node->y + (i & 2 ? 1.f : -1.f) * child_size / 2.f;
        node->children[i] = new QuadtreeNode{x, y, child_size, node->depth + 1};
        node->children[i]->ul_id = node->ul_id;
        node->children[i]->ur_id = node->ur_id;
        node->children[i]->vl_id = node->vl_id;
        node->children[i]->vr_id = node->vr_id;
        node->cached_neighbors.resize(8);
    }
}
void QTSimulator::smoothing(QuadtreeNode *new_root) {

    std::vector<QuadtreeNode *> leaves;
    collectLeafNodes(new_root, leaves);

    std::priority_queue<
        QuadtreeNode *,
        std::vector<QuadtreeNode *>,
        QTNLeafNodeFirst>
        pq(QTNLeafNodeFirst(), std::move(leaves));

    while (!pq.empty()) {
        auto node = pq.top();
        pq.pop();

        if (!node->is_leaf)
            continue;

        float step = node->size / 2.f + 0.1f;
        int direction[5] = {1, 0, -1, 0, 1};
        for (int i = 0; i < 4; i++) {
            float dx = direction[i] * step;
            float dy = direction[i + 1] * step;

            float x = node->x + dx;
            float y = node->y + dy;

            if (x < 0 || x > nx || y < 0 || y > ny)
                continue;

            auto neighbor = getNodeAt(new_root, x, y);

            if (neighbor && node->depth - neighbor->depth > 1) {
                subdivideNode(neighbor);
                for (int c = 0; c < 4; c++) {
                    neighbor->children[c]->phi = neighbor->phi;
                    pq.push(neighbor->children[c]);
                }
            }
        }
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
    if (node->size < 1.5f && node->phi < 0) {
        particles.push_back({node->x, node->y});
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

QuadtreeNode *
QTSimulator::getNodeAt(QuadtreeNode *node, float x, float y) const {
    if (x < 0 || x > nx || y < 0 || y > ny)
        return nullptr;
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

void QTSimulator::collectLeafNodes(
    QuadtreeNode *node, std::vector<QuadtreeNode *> &leaves
) const {
    if (node->is_leaf) {
        leaves.push_back(node);
        return;
    }
    for (int i = 0; i < 4; i++) {
        collectLeafNodes(node->children[i], leaves);
    }
}

// ######################### MLS #########################
struct SamplePoint {
    float x, y;    // 採樣點的物理座標（細胞中心或邊中心）
    float h;       // 該點對應的網格尺寸
    float val;     // 該點的數值
    float val_old; // 該點的數值
};

// 1D 變數的 MLS 求解輔助函式
std::pair<float, float>
solveMLS(float x, float y, const std::vector<SamplePoint> &samples) {
    if (samples.empty())
        return {0, 0};

    Eigen::Matrix3f A = Eigen::Matrix3f::Zero();
    Eigen::Matrix<float, 3, 2> b;
    b.setZero();

    for (const auto &p : samples) {
        float dx = std::abs(p.x - x);
        float dy = std::abs(p.y - y);
        float h = p.h;

        float eps = 1e-2f;
        float wx = std::max(1.f - (dx / h), eps);
        float wy = std::max(1.f - (dy / h), eps);
        float weight = wx * wy;

        float lx = p.x - x;
        float ly = p.y - y;
        Eigen::Vector3f z_i(lx, ly, 1.f);

        A += weight * (z_i * z_i.transpose());
        b.col(0) += weight * z_i * p.val;
        b.col(1) += weight * z_i * p.val_old;
    }

    // 脊迴歸懲罰項，確保 A 在共線或極端情況下仍可逆
    A += Eigen::Matrix3f::Identity() * 1e-6f;

    Eigen::MatrixXf c = A.ldlt().solve(b);
    return {
        c(2, 0), c(2, 1)
    }; // 局部座標系下，常數項即為 (0,0) 擬合值，對應索引 2
}

InterpolatedData QTSimulator::MLSinterpolate(float x, float y) {
    InterpolatedData result = {0.0f, 1000.f, 0.0f, 0.0f, 0.0f};

    QuadtreeNode *target = getNodeAt(root, x, y);
    if (!target)
        return result;

    float search_radius = target->size * 3.f;
    std::vector<QuadtreeNode *> neighbors;
    getNodesIn(root, x, y, search_radius, neighbors);

    float min_x = 0.0f;
    float max_x = nx;
    float min_y = 0.0f;
    float max_y = ny;

    // ==========================================================
    // 1. 細胞中心變數插值 (S, phi, pressure)
    // ==========================================================
    struct CellNode {
        float x, y, h;
        float S, phi, pressure;
    };
    std::vector<CellNode> cell_points;

    for (auto n : neighbors) {
        cell_points.push_back({n->x, n->y, n->size, n->S, n->phi, n->pressure});

        // 檢查細胞邊界鏡像 (Mirroring)
        float dist_left = n->x - min_x;
        float dist_right = max_x - n->x;
        float dist_top = n->y - min_y;
        float dist_bottom = max_y - n->y;
        float threshold = n->size * 1.5f;

        bool mirror_x = false;
        float mx = n->x;
        if (dist_left < threshold) {
            mirror_x = true;
            mx = min_x - dist_left;
        } else if (dist_right < threshold) {
            mirror_x = true;
            mx = max_x + dist_right;
        }

        bool mirror_y = false;
        float my = n->y;
        if (dist_top < threshold) {
            mirror_y = true;
            my = min_y - dist_top;
        } else if (dist_bottom < threshold) {
            mirror_y = true;
            my = max_y + dist_bottom;
        }

        if (mirror_x)
            cell_points.push_back(
                {mx, n->y, n->size, n->S, n->phi, n->pressure}
            );
        if (mirror_y)
            cell_points.push_back(
                {n->x, my, n->size, n->S, n->phi, n->pressure}
            );
        if (mirror_x && mirror_y)
            cell_points.push_back({mx, my, n->size, n->S, n->phi, n->pressure});
    }

    // 求解細胞中心變數
    Eigen::Matrix3f A_cell = Eigen::Matrix3f::Zero();
    Eigen::Matrix<float, 3, 3> b_cell;
    b_cell.setZero();

    for (const auto &p_node : cell_points) {
        float dx = std::abs(p_node.x - x);
        float dy = std::abs(p_node.y - y);
        float h = p_node.h;

        float eps = 1e-2f;
        float wx = std::max(1.f - (dx / h), eps);
        float wy = std::max(1.f - (dy / h), eps);
        float weight = wx * wy;

        float lx = p_node.x - x;
        float ly = p_node.y - y;
        Eigen::Vector3f z_i(lx, ly, 1.f);

        A_cell += weight * (z_i * z_i.transpose());
        b_cell.col(0) += weight * z_i * p_node.S;
        b_cell.col(1) += weight * z_i * p_node.phi;
        b_cell.col(2) += weight * z_i * p_node.pressure;
    }
    A_cell += Eigen::Matrix3f::Identity() * 1e-6f;
    Eigen::MatrixXf c_cell = A_cell.ldlt().solve(b_cell);

    result.S = c_cell(2, 0);
    result.phi = c_cell(2, 1);
    result.pressure = c_cell(2, 2);

    // ==========================================================
    // 2. u 速度插值 (垂直邊上的速度)
    // ==========================================================
    std::vector<SamplePoint> u_samples;
    std::unordered_set<int> visited_u_faces;

    for (auto n : neighbors) {
        // 假設你的 Node 中存有左右垂直邊的面索引 (face indices)
        // 論文提到 T-junction 會自動對齊到大 parent face 取得單一速度值
        int left_face = n->ul_id;
        int right_face = n->ur_id;

        auto add_u_face = [&](int face_idx) {
            if (face_idx == -1 || visited_u_faces.count(face_idx))
                return;
            visited_u_faces.insert(face_idx);

            auto &face = QTu[face_idx];
            float val = face.val;
            float val_old = face.val_old;
            float face_x = face.x;
            float face_y = face.y;
            float face_h = face.length;
            u_samples.push_back({face_x, face_y, face_h, val, val_old});

            // 垂直邊邊界鏡像 (以物理邊界條件處理速度)
            float dist_left = face_x - min_x;
            float dist_right = max_x - face_x;
            float dist_top = face_y - min_y;
            float dist_bottom = max_y - face_y;
            float threshold = face_h * 1.5f;

            bool mirror_x = false;
            float mx = face_x;
            float m_val_x = val;
            float m_val_old_x = val_old;
            if (dist_left < threshold) {
                mirror_x = true;
                mx = min_x - dist_left;
                m_val_x = -val;         // 固體邊界法向速度反向 (消去穿透流)
                m_val_old_x = -val_old; // 固體邊界法向速度反向 (消去穿透流)
            } else if (dist_right < threshold) {
                mirror_x = true;
                mx = max_x + dist_right;
                m_val_x = -val;
                m_val_old_x = -val_old;
            }

            bool mirror_y = false;
            float my = face_y;
            float m_val_y = val;
            float m_val_old_y = val_old;
            if (dist_top < threshold) {
                mirror_y = true;
                my = min_y - dist_top;
                m_val_y = val;         // 滑動邊界 (Slip BC) 切向速度同向
                m_val_old_y = val_old; // 滑動邊界 (Slip BC) 切向速度同向
            } else if (dist_bottom < threshold) {
                mirror_y = true;
                my = max_y + dist_bottom;
                m_val_y = val;
                m_val_old_y = val_old;
            }

            if (mirror_x)
                u_samples.push_back({mx, face_y, face_h, m_val_x, m_val_old_x});
            if (mirror_y)
                u_samples.push_back({face_x, my, face_h, m_val_y, m_val_old_y});
            if (mirror_x && mirror_y)
                u_samples.push_back({mx, my, face_h, -m_val_x, -m_val_old_x});
        };

        // 垂直邊的 X 座標在細胞中心左右，Y 座標與細胞同高
        float h = n->size;
        add_u_face(left_face);
        add_u_face(right_face);
    }
    auto tmp_result_u = solveMLS(x, y, u_samples);
    result.u = tmp_result_u.first;
    result.u_old = tmp_result_u.second;

    // ==========================================================
    // 3. v 速度插值 (水平邊上的速度)
    // ==========================================================
    std::vector<SamplePoint> v_samples;
    std::unordered_set<int> visited_v_faces;

    for (auto n : neighbors) {
        int top_face = n->vl_id;
        int bottom_face = n->vr_id;

        auto add_v_face = [&](int face_idx) {
            if (face_idx == -1 || visited_v_faces.count(face_idx))
                return;
            visited_v_faces.insert(face_idx);

            auto &face = QTv[face_idx];
            float val = face.val;         // 全域速度陣列
            float val_old = face.val_old; // 全域速度陣列
            float face_x = face.x;
            float face_y = face.y;
            float face_h = face.length;
            v_samples.push_back({face_x, face_y, face_h, val, val_old});

            // 水平邊邊界鏡像
            float dist_left = face_x - min_x;
            float dist_right = max_x - face_x;
            float dist_top = face_y - min_y;
            float dist_bottom = max_y - face_y;
            float threshold = face_h * 1.5f;

            bool mirror_x = false;
            float mx = face_x;
            float m_val_x = val;
            float m_val_old_x = val_old;
            if (dist_left < threshold) {
                mirror_x = true;
                mx = min_x - dist_left;
                m_val_x = val;         // 滑動邊界 (Slip BC) 切向速度同向
                m_val_old_x = val_old; // 滑動邊界 (Slip BC) 切向速度同向
            } else if (dist_right < threshold) {
                mirror_x = true;
                mx = max_x + dist_right;
                m_val_x = val;
                m_val_old_x = val_old;
            }

            bool mirror_y = false;
            float my = face_y;
            float m_val_y = val;
            float m_val_old_y = val_old;
            if (dist_top < threshold) {
                mirror_y = true;
                my = min_y - dist_top;
                m_val_y = -val;         // 固體邊界法向速度反向 (消去穿透流)
                m_val_old_y = -val_old; // 固體邊界法向速度反向 (消去穿透流)
            } else if (dist_bottom < threshold) {
                mirror_y = true;
                my = max_y + dist_bottom;
                m_val_y = -val;
                m_val_old_y = -val_old;
            }

            if (mirror_x)
                v_samples.push_back({mx, face_y, face_h, m_val_x, m_val_old_x});
            if (mirror_y)
                v_samples.push_back({face_x, my, face_h, m_val_y, m_val_old_y});
            if (mirror_x && mirror_y)
                v_samples.push_back({mx, my, face_h, -m_val_y, -m_val_old_y});
        };

        // 水平邊的 Y 座標在細胞中心上下，X 座標與細胞同寬
        float h = n->size;
        add_v_face(bottom_face);
        add_v_face(top_face);
    }
    auto tmp_result_v = solveMLS(x, y, v_samples);
    result.v = tmp_result_v.first;
    result.v_old = tmp_result_v.second;

    return result;
}

void QTSimulator::findAllEdges() {
    QTu_new.clear();
    QTv_new.clear();
    for (auto leaf : cached_leaves) {

        float half_size = leaf->size / 2.f;
        float quad_size = leaf->size / 4.f;
        float step = half_size + 0.1f;

        int node_idx = 0;
        auto &nodes = leaf->cached_neighbors;
        for (int dir = 0; dir < 4; dir++) {
            int neighbot_cnt = leaf->neighbor_cnt[dir];

            float x = leaf->x, y = leaf->y, solid_frac = 0;
            int face_sgn;
            std::vector<QuadtreeNode *> adj_cells = {leaf};
            std::vector<float> coefs;

            // up or down, need update y
            if (dir & 1) {
                y += (dir & 2 ? half_size : -half_size);
            }
            // left or right, need update x
            else {
                x += (dir & 2 ? half_size : -half_size);
            }

            if (dir == 0)
                leaf->ul_id = QTu_new.size();
            if (dir == 1)
                leaf->vl_id = QTv_new.size();
            if (dir == 2)
                leaf->ur_id = QTu_new.size();
            if (dir == 3)
                leaf->vr_id = QTv_new.size();

            if (dir <= 1) // left and up
                face_sgn = -1;
            else // right and down
                face_sgn = 1;

            if (neighbot_cnt == 0) {
                coefs.assign({0});
                solid_frac = 1;
            } else if (neighbot_cnt == 2) {
                adj_cells.push_back(nodes[node_idx]);
                adj_cells.push_back(nodes[node_idx + 1]);

                float coef = face_sgn * 1.f / (1.5f * nodes[node_idx]->size);
                coefs.assign({-1.f * coef, 0.5f * coef, 0.5f * coef});
            } else if (neighbot_cnt == 1) {
                if (nodes[node_idx]->depth == leaf->depth) {
                    if (dir <= 1)
                        goto OVERLAY_NODE;

                    adj_cells.push_back(nodes[node_idx]);
                    float coef = face_sgn * 1.f / leaf->size;
                    coefs.assign({-coef, coef});
                    goto UPDATE;
                }

            OVERLAY_NODE:
                switch (dir) {
                case 0:
                    leaf->ul_id = nodes[0]->ur_id;
                    break;
                case 1:
                    leaf->vl_id = nodes[0]->vr_id;
                    break;
                case 2:
                    leaf->ur_id = nodes[0]->ul_id;
                    break;
                case 3:
                    leaf->vr_id = nodes[0]->vl_id;
                    break;
                }
                node_idx += neighbot_cnt;
                continue;
            } else {
                std::cout << "Some ERROR appear" << std::endl;
                continue;
            }

        UPDATE:
            // up or down, push into v
            if (dir & 1)
                QTv_new.push_back({
                    x,
                    y,
                    leaf->size,
                    0,
                    solid_frac,
                    adj_cells,
                    coefs,
                });
            // left or right, push into u
            else
                QTu_new.push_back({
                    x,
                    y,
                    leaf->size,
                    0,
                    solid_frac,
                    adj_cells,
                    coefs,
                });

            node_idx += neighbot_cnt;
        }
    }
}

float QTSimulator::getVelocity(std::vector<QuadtreeEdge> &field, int id) {
    if (id == -1)
        return 0.0f;
    return field[id].val;
}

void QTSimulator::cacheNeighbors(QuadtreeNode *new_root) {
    for (auto leaf : cached_leaves) {

        std::vector<std::pair<int, QuadtreeNode *>> neighbors;
        getNeighbors(new_root, leaf, neighbors);

        std::vector<QuadtreeNode *> direc_neigh[4];
        for (auto [dir, node] : neighbors) {
            leaf->cached_neighbors.push_back(node);
            leaf->neighbor_cnt[dir]++;
        }
    }
}

void QTSimulator::cacheLeaves(QuadtreeNode *new_root) {
    cached_leaves.clear();
    collectLeafNodes(new_root, cached_leaves);

    sort(
        cached_leaves.begin(),
        cached_leaves.end(),
        [](QuadtreeNode *a, QuadtreeNode *b) {
            if (a->depth != b->depth)
                return a->depth < b->depth;
            if (a->x != b->x)
                return a->x < b->x;
            return a->y < b->y;
        }
    );
}
