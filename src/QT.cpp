#include "QT.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <math.h>
#include <omp.h>
#include <queue>
#include <random>
#include <utility>
#include <vector>

// PUBLIC
QTSimulator::QTSimulator(int width, int height) : nx(width), ny(height) {
    phash_head.resize(nx * ny, -1);
    phash_next.clear();
    leaf_table[0].clear();
    leaf_table[1].clear();
    G     = 150.f;
    Sigma = 0.0f;

    // solver.setMaxIterations(80);
    // solver.setTolerance(1e-3);
    // Eigen::setNbThreads(6);

    root_list = 0;
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, root_list);
    initQuadtree(3, root_list);
}
void QTSimulator::reset() {
    particles.clear();
    phash_head.assign(nx * ny, -1);
    phash_next.clear();
    node_pool[0].clear();
    node_pool[1].clear();
    leaf_table[0].clear();
    leaf_table[1].clear();

    root_list = 0;
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, root_list);
    initQuadtree(3, root_list);
}

void QTSimulator::addWater(float x, float y, float radius) {
    recursiveUpdatePhi(x, y, radius, false, root_list);
    leaf_table[root_list].clear();
    smoothing(root_list);
    cacheLeaves(root_list);
    cacheNeighbors(root_list);
    updateParticleIdx();
}
void QTSimulator::delWater(float x, float y, float radius) {
    float r2 = radius * radius;
    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](auto &p) {
            return ((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y)) < r2;
        });
    particles.erase(new_end, particles.end());
    recursiveUpdatePhi(x, y, radius, true, root_list);
    leaf_table[root_list].clear();
    smoothing(root_list);
    cacheLeaves(root_list);
    cacheNeighbors(root_list);
    updateParticleIdx();
}

void QTSimulator::recursiveGetLines(
    std::vector<Line> &lines, int list_idx, int node_idx
) const {
    auto &node = getNode(list_idx, node_idx);
    if (node.is_leaf) {
        float x1 = node.x - node.size / 2.f;
        float y1 = node.y - node.size / 2.f;
        float x2 = node.x + node.size / 2.f;
        float y2 = node.y + node.size / 2.f;

        lines.push_back({x1, y1, x2, y1});
        lines.push_back({x1, y1, x1, y2});
        lines.push_back({x1, y2, x2, y2});
        lines.push_back({x2, y1, x2, y2});
        return;
    }

    for (int i = 0; i < 4; i++) {
        recursiveGetLines(lines, list_idx, node.children_idx[i]);
    }
}
std::vector<Line> QTSimulator::getLines() const {
    std::vector<Line> lines;

    recursiveGetLines(lines, root_list);

    return lines;
}

void QTSimulator::update(float dt) {
    max_u = 0.0f;
    max_v = 0.0f;

    // Generate new Quad Tree
    computeSizingFunction(dt);
    propagateSizingFunction();

    int new_root_list = 1 - root_list;
    node_pool[new_root_list].clear();
    leaf_table[new_root_list].clear();
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, new_root_list);
    recursiveBuildTree(dt, new_root_list);

    smoothing(new_root_list);
    cacheLeaves(new_root_list);
    cacheNeighbors(new_root_list);

    findAllEdges(new_root_list);

    advectQuadtreeDatas(dt, new_root_list);
    root_list = new_root_list;

    QTu = QTu_new;
    QTv = QTv_new;

    particleToGrid();

    for (auto &face : QTu)
        face.val_old = face.val;
    for (auto &face : QTv)
        face.val_old = face.val;

    QTapplyGravity(dt);

    setBoundaries();
    QTproject();
    setBoundaries();

    velExtrapolation();

    gridToParticle();
    advectParticles(dt);

    redistancing();
    resampleParticles();

    // Redistance
    updateParticleIdx();
    // reconstructSurface();
    // redistancing();
}

void QTSimulator::setBoundaries() {
    /*
     * For all edges.
     * Set the boundary to 0
     */
    for (int i = 0; i < QTu.size(); i++) {
        if (QTu[i].solid_fraction >= 1.f) {
            QTu[i].val     = 0.0f;
            QTu[i].val_old = 0.0f;
        }
    }

    for (int i = 0; i < QTv.size(); i++) {
        if (QTv[i].solid_fraction >= 1.f) {
            QTv[i].val     = 0.0f;
            QTv[i].val_old = 0.0f;
        }
    }
}

void QTSimulator::computeSizingFunction(float dt) {
    /*
     * For all leaves.
     * Compute the size function value
     */

    float T1 = 0.9, T2 = 0.01;
    float R_t = std::pow(T1, dt / T2);

#pragma omp parallel for
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        if (std::abs(leaf.phi) < leaf.size * 1.5f) {
            float gamma_phi = 4.f;
            float gamma_u   = 3.f;

            // // 0: left
            // // 1: up
            // // 2: right
            // // 3: down
            // NeighborData L = getNeighborData(leaf, 0);
            // NeighborData U = getNeighborData(leaf, 1);
            // NeighborData R = getNeighborData(leaf, 2);
            // NeighborData D = getNeighborData(leaf, 3);
            //
            // float f_C = leaf.phi;
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
            float geom_term = gamma_phi / (std::abs(leaf.phi) + 0.1f);

            float ul = getVelocity(QTu, leaf.ul_id);
            float ur = getVelocity(QTu, leaf.ur_id);
            float vl = getVelocity(QTv, leaf.vl_id);
            float vr = getVelocity(QTv, leaf.vr_id);

            float du_dx = (ul - ur) / leaf.size;
            float dv_dy = (vl - vr) / leaf.size;

            float vel_term = gamma_u * std::sqrt(du_dx * du_dx + dv_dy * dv_dy);

            float S = geom_term + vel_term;

            auto data    = advect(leaf.x, leaf.y, dt, OPT_S);
            float S_star = data.S;

            leaf.S_new = std::max(R_t * S_star, S);
        } else {
            leaf.S_new = 0.0f;
        }
    }

#pragma omp parallel for
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        leaf.S     = leaf.S_new;
    }
}

void QTSimulator::propagateSizingFunction() {
    /*
     * For all leaves.
     * propagate the size function value to around nodes
     */

    int iterations = 5;
#pragma omp parallel
    for (int iter = 0; iter < iterations; iter++) {
#pragma omp for
        for (int leaf_idx : cached_leaves_idx) {
            auto &leaf = getNode(root_list, leaf_idx);

            float self_area = leaf.size * leaf.size;
            leaf.S_new      = leaf.S * self_area;

            float total_area = self_area;
            for (int i = 0; i < leaf.cached_neighbors_cnt; i++) {
                int neighbor_idx = leaf.cached_neighbors_idx[i];
                auto &neighbor   = getNode(root_list, neighbor_idx);

                float area = neighbor.size * neighbor.size;
                leaf.S_new += std::max(neighbor.S, leaf.S) * area;
                total_area += area;
            }

            leaf.S_new /= total_area;
        }

#pragma omp for
        for (auto leaf_idx : cached_leaves_idx) {
            auto &leaf = getNode(root_list, leaf_idx);
            leaf.S     = leaf.S_new;
        }
    }
}

void QTSimulator::recursiveBuildTree(float dt, int list_idx, int node_idx) {
    /*
     * Recursive
     * build the tree until size to 1.f
     */

    auto &node = getNode(list_idx, node_idx);
    if (node.size < 1.5f)
        return;

    auto advected_data = advect(node.x, node.y, dt, OPT_PHI);
    auto data          = MLSinterpolate(node.x, node.y, OPT_S);

    float exp_phi = advected_data.phi;
    float exp_S   = data.S;

    // if (std::abs(exp_phi) < node.size && exp_S > (1.f / node.size)) {
    if (std::abs(exp_phi) < node.size) {
        subdivideNode(list_idx, node_idx);

        for (int i = 0; i < 4; i++) {
            recursiveBuildTree(
                dt, list_idx, getNode(list_idx, node_idx).children_idx[i]
            );
        }
    }
}

void QTSimulator::findAllEdges(int list_idx) {
    /*
     * For all leaves.
     * find all leaves' edge
     */

    QTu_new.clear();
    QTv_new.clear();
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(list_idx, leaf_idx);

        float half_size = leaf.size / 2.f;
        float quad_size = leaf.size / 4.f;
        float step      = half_size + 0.1f;

        int neighbor_idx = 0;
        auto &nodes_idx  = leaf.cached_neighbors_idx;
        for (int dir = 0; dir < 4; dir++) {
            int neighbot_cnt = leaf.neighbor_cnt[dir];

            float x = leaf.x, y = leaf.y, solid_frac = 0;
            int face_sgn;
            std::vector<int> adj_cells_idx = {leaf_idx};
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
                leaf.ul_id = QTu_new.size();
            if (dir == 1)
                leaf.vl_id = QTv_new.size();
            if (dir == 2)
                leaf.ur_id = QTu_new.size();
            if (dir == 3)
                leaf.vr_id = QTv_new.size();

            if (dir <= 1) // left and up
                face_sgn = -1;
            else // right and down
                face_sgn = 1;

            if (neighbot_cnt == 0) {
                coefs.assign({0});
                solid_frac = 1;
            } else if (neighbot_cnt == 2) {
                auto &n_0 = getNode(list_idx, nodes_idx[neighbor_idx]);
                adj_cells_idx.push_back(nodes_idx[neighbor_idx]);
                adj_cells_idx.push_back(nodes_idx[neighbor_idx + 1]);

                float coef = face_sgn * 1.f / (1.5f * n_0.size);
                coefs.assign({-1.f * coef, 0.5f * coef, 0.5f * coef});
            } else if (neighbot_cnt == 1) {
                auto &n_0 = getNode(list_idx, nodes_idx[neighbor_idx]);
                if (n_0.depth == leaf.depth) {
                    if (dir <= 1)
                        goto OVERLAY_NODE;

                    adj_cells_idx.push_back(nodes_idx[neighbor_idx]);
                    float coef = face_sgn * 1.f / leaf.size;
                    coefs.assign({-coef, coef});
                    goto UPDATE;
                }

            OVERLAY_NODE:
                switch (dir) {
                case 0:
                    leaf.ul_id = n_0.ur_id;
                    break;
                case 1:
                    leaf.vl_id = n_0.vr_id;
                    break;
                case 2:
                    leaf.ur_id = n_0.ul_id;
                    break;
                case 3:
                    leaf.vr_id = n_0.vl_id;
                    break;
                }
                neighbor_idx += neighbot_cnt;
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
                    leaf.size,
                    0,
                    solid_frac,
                    adj_cells_idx,
                    coefs,
                });
            // left or right, push into u
            else
                QTu_new.push_back({
                    x,
                    y,
                    leaf.size,
                    0,
                    solid_frac,
                    adj_cells_idx,
                    coefs,
                });

            neighbor_idx += neighbot_cnt;
        }
    }
}

void QTSimulator::advectQuadtreeDatas(float dt, int list_idx) {
    /*
     * For all leaves.
     * For all edges.
     * advect the Phi, Size function and velocity data
     */

#pragma omp parallel for
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf         = getNode(list_idx, leaf_idx);
        auto advected_leaf = advect(leaf.x, leaf.y, dt, OPT_CELL_ALL);

        leaf.phi = advected_leaf.phi;
        leaf.S   = advected_leaf.S;
    }

#pragma omp parallel for
    for (auto &face : QTu_new) {
        auto advected_face = advect(face.x, face.y, dt, OPT_U_ALL);
        face.val           = advected_face.u;
    }
#pragma omp parallel for
    for (auto &face : QTv_new) {
        auto advected_face = advect(face.x, face.y, dt, OPT_V_ALL);
        face.val           = advected_face.v;
    }
}

void QTSimulator::particleToGrid() {
    /*
     * For all edges.
     * get velocity data from particles
     */

    float search_radius = 1.5f;
    float weight_eps    = 1e-6;

    auto updateFace = [&](QuadtreeEdge &face, bool is_u) {
        if (face.length > 1.5f)
            return;

        float sum_w   = 0;
        float sum_val = 0;

        int ix = std::clamp((int)face.x, 0, nx - 1);
        int iy = std::clamp((int)face.y, 0, ny - 1);

        for (int dj = -2; dj <= 2; dj++) {
            for (int di = -2; di <= 2; di++) {
                int ni = ix + di, nj = iy + dj;

                if (ni < 0 || ni >= nx || nj < 0 || nj >= ny)
                    continue;

                int cell_idx = IX(ni, nj);
                for (int idx = phash_head[cell_idx]; idx != -1;
                     idx     = phash_next[idx]) {
                    auto &p = particles[idx];

                    float dx = std::abs(p.x - face.x);
                    float dy = std::abs(p.y - face.y);

                    float wx = std::max(0.f, 1.f - dx / search_radius);
                    float wy = std::max(0.f, 1.f - dy / search_radius);
                    float w  = wx * wy;

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

#pragma omp parallel for
    for (int i = 0; i < QTu.size(); i++)
        updateFace(QTu[i], true);

#pragma omp parallel for
    for (int i = 0; i < QTv.size(); i++)
        updateFace(QTv[i], false);
}

void QTSimulator::QTapplyGravity(float dt) {
    /*
     * For all edges.
     * add gravity on horizontial edges
     */

#pragma omp parallel for
    for (auto &face : QTv) {
        bool is_water = false;
        for (int node_idx : face.adj_cells_idx)
            is_water |= (getNode(root_list, node_idx).phi < 0.0f);

        if (is_water)
            face.val += G * dt;
    }
}

void QTSimulator::QTproject() {
    /*
     * For all edges.
     * compute pressure and update velocity
     */

    auto compute_Wk = [&](QuadtreeEdge &face) {
        float total = 0, liquid = 0;
        for (int i = 0; i < face.adj_cells_idx.size(); i++) {
            auto &cell = getNode(root_list, face.adj_cells_idx[i]);
            total += face.grad_coeff[i] * cell.phi;
            if (cell.phi <= 0.f)
                liquid += face.grad_coeff[i] * cell.phi;
        }
        if (std::abs(liquid) < 1e-6f)
            return 1.f;
        return total / liquid;
    };

    int fluid_count = 0;
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        if (leaf.phi <= 0.0)
            leaf.fluid_id = fluid_count++;
        else
            leaf.fluid_id = -1;
    }

    if (fluid_count == 0)
        return;

    Eigen::VectorXf div(fluid_count);
    div.setZero();
    QTtriplets.clear();

    int max_threads = omp_get_max_threads();
    static std::vector<std::vector<Eigen::Triplet<float>>> thread_triplets(
        max_threads
    );

    auto addFace = [&](QuadtreeEdge &face,
                       std::vector<Eigen::Triplet<float>> &local_trip) {
        float u_star = face.val;
        float VA     = face.length * (1.f - face.solid_fraction);

        bool has_fluid = false;
        for (int cell_idx : face.adj_cells_idx)
            has_fluid |= (getNode(root_list, cell_idx).fluid_id != -1);
        if (!has_fluid)
            return;

        float W_k         = compute_Wk(face);
        float W_prime     = std::max(W_k, 0.01f);
        float face_weight = VA * W_prime;

        for (int i = 0; i < face.adj_cells_idx.size(); i++) {
            auto &cell_i = getNode(root_list, face.adj_cells_idx[i]);
            if (cell_i.fluid_id == -1)
                continue;
            int row = cell_i.fluid_id;

#pragma omp atomic
            div[row] += VA * u_star * face.grad_coeff[i];

            for (int j = 0; j < face.adj_cells_idx.size(); j++) {
                auto &cell_j = getNode(root_list, face.adj_cells_idx[j]);
                if (cell_j.fluid_id == -1)
                    continue;
                int col = cell_j.fluid_id;

                float val =
                    face_weight * face.grad_coeff[i] * face.grad_coeff[j];
                local_trip.push_back({row, col, val});
            }
        }
    };

#pragma omp parallel
    {
        int tid       = omp_get_thread_num();
        auto &local_t = thread_triplets[tid];
        local_t.clear();

#pragma omp for nowait
        for (auto &face : QTu)
            addFace(face, local_t);
#pragma omp for nowait
        for (auto &face : QTv)
            addFace(face, local_t);
    }

    for (const auto &local_t : thread_triplets)
        QTtriplets.insert(QTtriplets.end(), local_t.begin(), local_t.end());

    // Eigen solvers
    Eigen::SparseMatrix<float> A(fluid_count, fluid_count);
    A.setFromTriplets(QTtriplets.begin(), QTtriplets.end());
    solver.compute(A);
    Eigen::VectorXf pressure = solver.solve(div);

    auto updateFace = [&](QuadtreeEdge &face) {
        bool has_fluid = false;
        for (int cell_idx : face.adj_cells_idx)
            has_fluid |= (getNode(root_list, cell_idx).fluid_id != -1);
        if (!has_fluid)
            return;

        float grad_p = 0;

        for (int i = 0; i < face.adj_cells_idx.size(); i++) {
            auto &cell = getNode(root_list, face.adj_cells_idx[i]);
            float p_val =
                (cell.fluid_id != -1 ? pressure[cell.fluid_id] : 0.0f);
            grad_p += face.grad_coeff[i] * p_val;
        }

        float W_prime = std::max(compute_Wk(face), 0.01f);
        face.val -= W_prime * grad_p;
    };

#pragma omp parallel for
    for (auto &face : QTu)
        updateFace(face);
#pragma omp parallel for
    for (auto &face : QTv)
        updateFace(face);
}

void QTSimulator::velExtrapolation() {

    std::queue<std::pair<int, int>> Q;
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);

        if (leaf.phi < 0) {
            leaf.known = true;
            if (leaf.size < 1.5f)
                Q.push({leaf_idx, 0});
        }
    }

    if (Q.empty()) {
        return;
    }

    while (!Q.empty()) {
        auto [curr_node_idx, curr_dist] = Q.front();
        Q.pop();

        auto &curr_node = getNode(root_list, curr_node_idx);

        int idx = 0;
        for (int dir = 0; dir < 4; dir++) {
            for (int i = 0; i < curr_node.neighbor_cnt[dir]; i++) {
                int neighbor_idx = curr_node.cached_neighbors_idx[idx + i];
                auto &neighbor   = getNode(root_list, neighbor_idx);

                int dist = curr_dist + 1;
                if (neighbor.known || dist > 4)
                    continue;

                if (dir != 0)
                    QTu[neighbor.ur_id].val = QTu[curr_node.ur_id].val;
                if (dir != 1)
                    QTv[neighbor.vr_id].val = QTv[curr_node.vr_id].val;
                if (dir != 2)
                    QTu[neighbor.ul_id].val = QTu[curr_node.ul_id].val;
                if (dir != 3)
                    QTv[neighbor.vl_id].val = QTv[curr_node.vl_id].val;
                Q.push({neighbor_idx, dist});
                neighbor.known = true;
            }
            idx += curr_node.neighbor_cnt[dir];
        }
    }
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        leaf.known = false;
    }
}

void QTSimulator::gridToParticle() {
    /*
     * For all particles
     * get velocity from edges
     */

    float flipRatio = 0.95f;

#pragma omp parallel for
    for (auto &p : particles) {
        auto data        = MLSinterpolate(p.x, p.y, OPT_VEL_ALL);
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

void QTSimulator::advectParticles(float dt) {
    /*
     * For all particles
     * advect the velocity
     */

#pragma omp parallel for
    for (auto &p : particles) {
        auto d1     = MLSinterpolate(p.x, p.y, OPT_VEL_ALL);
        float mid_x = p.x + d1.u * dt * 0.5f;
        float mid_y = p.y + d1.v * dt * 0.5f;

        auto d2 = MLSinterpolate(mid_x, mid_y, OPT_VEL_ALL);
        p.x += d2.u * dt;
        p.y += d2.v * dt;

        float eps = 1e-3;
        p.x       = std::clamp(p.x, eps, (float)nx - eps);
        p.y       = std::clamp(p.y, eps, (float)ny - eps);
    }
}

void QTSimulator::resampleParticles() {
    /*
     * For all leaves.
     * resample particles
     */

    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](const auto &p) {
            int node_idx = getNodeIdxAt(p.x, p.y, root_list);

            if (node_idx == -1)
                return true;

            auto &node = getNode(root_list, node_idx);
            if (node.size > 1.5f)
                return true;

            if (node.phi < -3.f * node.size)
                return true;

            node.particle_cnt++;
            return false;
        });
    particles.erase(new_end, particles.end());

    int max_threads = omp_get_max_threads();
    static std::vector<std::vector<Particle>> thread_particles(max_threads);

#pragma omp parallel
    {
        int tid       = omp_get_thread_num();
        auto &local_p = thread_particles[tid];
        local_p.clear();

        thread_local std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<float> dis(-0.5f, 0.5f);

        const int target_ppc = 4;

#pragma omp for
        for (int leaf_idx : cached_leaves_idx) {
            auto &leaf = getNode(root_list, leaf_idx);
            if (leaf.phi < 0 &&
                leaf.size < 1.5f &&
                leaf.phi > -3.f * leaf.size &&
                leaf.particle_cnt < target_ppc) {

                int num_to_seed = target_ppc - leaf.particle_cnt;

                int seeded = 0;
                int tries  = 0;
                // 使用限制次數的嘗試，確保粒子能精準落在等值面內部
                while (seeded < num_to_seed && tries < 15) {
                    tries++;
                    float px = leaf.x + dis(gen) * leaf.size;
                    float py = leaf.y + dis(gen) * leaf.size;

                    // 關鍵修改：插值得到該候選點的精準 phi 值
                    auto phi_data = MLSinterpolate(px, py, OPT_PHI);

                    // 只有當該位置確實位於流體內部（保留一點安全邊界，例如 -0.1
                    // * size）才允許生成
                    if (phi_data.phi < -0.1f * leaf.size) {
                        auto vel_data = MLSinterpolate(px, py, OPT_VEL_ALL);
                        local_p.push_back({px, py, vel_data.u, vel_data.v});
                        seeded++;
                    }
                }
                // for (int k = 0; k < num_to_seed; k++) {
                //     float px      = leaf.x + dis(gen) * leaf.size;
                //     float py      = leaf.y + dis(gen) * leaf.size;
                //     auto vel_data = MLSinterpolate(px, py, OPT_VEL_ALL);
                //     local_p.push_back({px, py, vel_data.u, vel_data.v});
                // }
            }
        }
    }

    for (const auto &local_p : thread_particles) {
        particles.insert(particles.end(), local_p.begin(), local_p.end());
    }
}

void QTSimulator::updateParticleIdx() {
    /*
     * For all particles
     * update grid hash data
     */

    phash_head.assign(nx * ny, -1);
    phash_next.resize(particles.size());
    for (int i = 0; i < particles.size(); i++) {
        auto &p = particles[i];
        int ix  = std::clamp((int)p.x, 0, nx - 1);
        int iy  = std::clamp((int)p.y, 0, ny - 1);

        int cell_idx         = IX(ix, iy);
        phash_next[i]        = phash_head[cell_idx];
        phash_head[cell_idx] = i;
    }
}

void QTSimulator::reconstructSurface() {
#pragma omp parallel for
    for (int i = 0; i < cached_leaves_idx.size(); i++) {
        int leaf_idx = cached_leaves_idx[i];
        auto &leaf   = getNode(root_list, leaf_idx);

        // 僅在最細的網格層級（FLIP 粒子存在的區域）進行重建
        if (leaf.size > 1.5f || leaf.phi < 0)
            continue;

        float min_dist = nearestParticleDistance(leaf.x, leaf.y);

        if (min_dist != std::numeric_limits<float>::infinity()) {
            leaf.phi = min_dist;
            // float phi_p = min_dist - particle_radius;
            //
            // if (leaf.phi > 0 && phi_p < 0) {
            //     // 1. 水花細節：網格原本判定為空氣（phi >
            //     // 0），但周圍有飛散的粒子。 這裡我們強迫網格認可粒子，將 phi
            //     // 改為負值，使其在物理與視覺上被當作水。
            //     leaf.phi = phi_p;
            // } else if (leaf.phi < 0) {
            //     // 2. 體積保持：在液體內部，結合粒子與網格的
            //     // SDF（取聯集，避免耗散萎縮）。
            //     leaf.phi = std::min(leaf.phi, phi_p);
            // }
        }
    }
}

void QTSimulator::redistancing() {
    /*
     * For all leaves.
     * recompute the phi by particles
     */

    using PQElement = std::pair<float, int>; // {abs_phi, node_idx}

    std::vector<PQElement> init_elements;
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);

        float abs_leaf_phi = std::abs(leaf.phi);

        for (int i = 0; i < leaf.cached_neighbors_cnt; i++) {
            int neigh_idx = leaf.cached_neighbors_idx[i];
            auto &neigh   = getNode(root_list, neigh_idx);

            if ((leaf.phi < 0.f) == (neigh.phi < 0.f))
                continue;

            float demon = abs_leaf_phi + std::abs(neigh.phi);
            if (demon < 1e-6)
                continue;

            float theta = abs_leaf_phi / demon;
            float L     = (leaf.size + neigh.size) * 0.5f;

            leaf.phi_new = std::min(leaf.phi_new, theta * L);
            leaf.known   = true;

            neigh.phi_new = std::min(neigh.phi_new, (1.f - theta) * L);
            neigh.known   = true;
        }
    }
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);

        float sign = (leaf.phi < 0 ? -1 : 1);
        if (leaf.known) {
            init_elements.push_back({std::abs(leaf.phi_new), leaf_idx});
            leaf.phi_new *= sign;
        } else
            leaf.phi_new = std::numeric_limits<float>::infinity();

        leaf.known = false;
    }

    if (init_elements.empty()) {
        return;
    }

    FMMSolver(init_elements);
}

void QTSimulator::FMMSolver(std::vector<std::pair<float, int>> &init_datas) {
    auto get_axis = [](const QuadtreeNode &a, const QuadtreeNode &b) {
        float dx = std::abs(a.x - b.x);
        float dy = std::abs(a.y - b.y);
        if (dx > dy)
            return 0; // x 軸相鄰
        return 1;     // y 軸相鄰
    };

    auto solve_eikonal_2d = [](float u1, float h1, float u2, float h2) {
        float inv_h1_sq = 1.0f / (h1 * h1);
        float inv_h2_sq = 1.0f / (h2 * h2);

        float A = inv_h1_sq + inv_h2_sq;
        float B = -2.0f * (u1 * inv_h1_sq + u2 * inv_h2_sq);
        float C = (u1 * u1 * inv_h1_sq) + (u2 * u2 * inv_h2_sq) - 1.0f;

        float disc = B * B - 4.0f * A * C;
        if (disc < 0.0f)
            return std::min(u1 + h1, u2 + h2); // 數值異常時退化為 1D
        return (-B + std::sqrt(disc)) / (2.0f * A);
    };

    auto calculate_new_abs_phi = [&](const QuadtreeNode &neighbor) {
        struct ActiveNeighbor {
            float u = std::numeric_limits<float>::infinity();
            float h = 0.0f;
        };
        ActiveNeighbor axis_data[3]; // 0: x, 1: y, 2: z

        for (int i = 0; i < neighbor.cached_neighbors_cnt; i++) {
            int neigh_idx     = neighbor.cached_neighbors_idx[i];
            auto &neigh_neigh = getNode(root_list, neigh_idx);
            if (!neigh_neigh.known)
                continue; // 只使用已確定的鄰居

            float u_cand = std::abs(neigh_neigh.phi_new);
            float h_cand = 0.5f * (neighbor.size +
                                   neigh_neigh.size); // 考慮八叉樹不同網格大小
            int axis     = get_axis(neighbor, neigh_neigh);

            if (u_cand < axis_data[axis].u) {
                axis_data[axis].u = u_cand;
                axis_data[axis].h = h_cand;
            }
        }

        // 收集有有效已知鄰居的軸向
        struct AxisInfo {
            float u, h;
        };
        std::vector<AxisInfo> active_axes;
        for (int i = 0; i < 3; i++) {
            if (axis_data[i].u != std::numeric_limits<float>::infinity()) {
                active_axes.push_back({axis_data[i].u, axis_data[i].h});
            }
        }

        if (active_axes.empty())
            return std::numeric_limits<float>::infinity();

        // 依距離由小到大排序 (確保滿足因果關係 $u > u_i$)
        std::sort(
            active_axes.begin(),
            active_axes.end(),
            [](const AxisInfo &a, const AxisInfo &b) { return a.u < b.u; }
        );

        // 1D 嘗試
        float u = active_axes[0].u + active_axes[0].h;

        // 2D 嘗試
        if (active_axes.size() >= 2 && u > active_axes[1].u) {
            u = solve_eikonal_2d(
                active_axes[0].u,
                active_axes[0].h,
                active_axes[1].u,
                active_axes[1].h
            );
        }
        return u;
    };

    using PQElement = std::pair<float, int>; // {abs_phi, node_idx}
    std::priority_queue<
        PQElement,
        std::vector<PQElement>,
        std::greater<PQElement>>
        pq(std::greater<PQElement>(), std::move(init_datas));

    while (!pq.empty()) {
        auto [curr_abs_phi, curr_node_idx] = pq.top();
        pq.pop();

        auto &curr_node = getNode(root_list, curr_node_idx);
        if (curr_node.known)
            continue;
        curr_node.known = true;

        for (int i = 0; i < curr_node.cached_neighbors_cnt; i++) {
            int neighbor_idx = curr_node.cached_neighbors_idx[i];
            auto &neighbor   = getNode(root_list, neighbor_idx);

            if (neighbor.known)
                continue;

            float new_abs_phi = calculate_new_abs_phi(neighbor);

            if (new_abs_phi < std::abs(neighbor.phi_new)) {
                float sign       = (neighbor.phi < 0.f ? -1.f : 1.f);
                neighbor.phi_new = sign * new_abs_phi;
                pq.push({new_abs_phi, neighbor_idx});
            }
        }
    }
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        leaf.phi   = leaf.phi_new;
        leaf.known = false;
    }
}

// Quad Tree Build Functions
void QTSimulator::initQuadtree(int max_depth, int list_idx, int node_idx) {
    if (node_pool[list_idx][node_idx].depth >= max_depth)
        return;

    subdivideNode(list_idx, node_idx);
    for (int i = 0; i < 4; i++) {
        initQuadtree(
            max_depth, list_idx, node_pool[list_idx][node_idx].children_idx[i]
        );
    }
}

void QTSimulator::subdivideNode(int list_idx, int node_idx) {

    auto &pool = node_pool[list_idx];
    pool.reserve(pool.size() + 4);

    auto &node   = pool[node_idx];
    node.is_leaf = false;

    float child_size = node.size / 2.f;
    for (int i = 0; i < 4; i++) {
        float x = node.x + (i & 1 ? 1.f : -1.f) * child_size / 2.f;
        float y = node.y + (i & 2 ? 1.f : -1.f) * child_size / 2.f;
        node.children_idx[i] =
            allocate(x, y, child_size, node.depth + 1, list_idx);

        auto &child = getNode(list_idx, node.children_idx[i]);
        child.ul_id = node.ul_id;
        child.ur_id = node.ur_id;
        child.vl_id = node.vl_id;
        child.vr_id = node.vr_id;
        child.phi   = node.phi;
        child.S     = node.S;
    }
}

void QTSimulator::smoothing(int list_idx) {

    std::vector<int> leaves;
    collectLeafNodes(leaves, list_idx);

    auto surfaceFirst = [&](int a_idx, int b_idx) {
        return std::abs(getNode(list_idx, a_idx).depth) <
               std::abs(getNode(list_idx, b_idx).depth);
    };
    std::priority_queue<int, std::vector<int>, decltype(surfaceFirst)> pq(
        surfaceFirst, std::move(leaves)
    );

    while (!pq.empty()) {
        int node_idx = pq.top();
        pq.pop();

        float node_x, node_y, node_size;
        int node_depth;
        bool node_is_leaf;
        {
            auto &node   = getNode(list_idx, node_idx);
            node_x       = node.x;
            node_y       = node.y;
            node_size    = node.size;
            node_depth   = node.depth;
            node_is_leaf = node.is_leaf;
        }

        if (!node_is_leaf)
            continue;

        float step       = node_size / 2.f + 0.1f;
        int direction[5] = {1, 0, -1, 0, 1};
        for (int i = 0; i < 4; i++) {
            float dx = direction[i] * step;
            float dy = direction[i + 1] * step;

            float x = node_x + dx;
            float y = node_y + dy;

            if (x < 0 || x > nx || y < 0 || y > ny)
                continue;

            int neighbor_idx = getNodeIdxAt(x, y, list_idx);
            if (neighbor_idx == -1)
                continue;

            if (node_depth - getNode(list_idx, neighbor_idx).depth > 1) {
                subdivideNode(list_idx, neighbor_idx);
                for (int c = 0; c < 4; c++) {
                    auto &neighbor = getNode(list_idx, neighbor_idx);
                    auto &child = getNode(list_idx, neighbor.children_idx[c]);
                    child.phi   = neighbor.phi;
                    pq.push(neighbor.children_idx[c]);
                }
            }
        }
    }
}

void QTSimulator::recursiveUpdatePhi(
    float cx, float cy, float radius, bool is_delete, int list_idx, int node_idx
) {

    float node_x, node_y;
    float new_phi, node_phi, node_size;
    bool node_is_leaf;
    int node_particle_cnt;
    {
        auto &node = getNode(list_idx, node_idx);
        new_phi    = circleSDF(cx, cy, radius, node.x, node.y);

        if (is_delete && node.phi <= 0.f)
            node.phi = std::max(node.phi, -new_phi);
        else
            node.phi = std::min(node.phi, new_phi);

        node_x            = node.x;
        node_y            = node.y;
        node_size         = node.size;
        node_phi          = node.phi;
        node_is_leaf      = node.is_leaf;
        node_particle_cnt = node.particle_cnt;
    }

    if (!node_is_leaf) {
        for (int i = 0; i < 4; i++)
            recursiveUpdatePhi(
                cx,
                cy,
                radius,
                is_delete,
                list_idx,
                getNode(list_idx, node_idx).children_idx[i]
            );
        return;
    }

    if (node_size > 1.f && std::abs(node_phi) < node_size) {
        subdivideNode(list_idx, node_idx);
        for (int i = 0; i < 4; i++)
            recursiveUpdatePhi(
                cx,
                cy,
                radius,
                is_delete,
                list_idx,
                getNode(list_idx, node_idx).children_idx[i]
            );
    }
    if (node_size < 1.5f && node_phi < 0 && node_particle_cnt == 0) {
        particles.push_back({node_x, node_y});
    }
}

void QTSimulator::collectLeafNodes(
    std::vector<int> &leaves_idx, int list_idx
) const {
    for (int i = 0; i < node_pool[list_idx].size(); i++) {
        auto &node = getNode(list_idx, i);
        if (node.is_leaf)
            leaves_idx.push_back(i);
    }
}

void QTSimulator::cacheLeaves(int list_idx) {
    cached_leaves_idx.clear();
    collectLeafNodes(cached_leaves_idx, list_idx);

    sort(
        cached_leaves_idx.begin(),
        cached_leaves_idx.end(),
        [&](int a_idx, int b_idx) {
            auto &a = getNode(list_idx, a_idx);
            auto &b = getNode(list_idx, b_idx);
            if (a.depth != b.depth)
                return a.depth < b.depth;
            if (a.x != b.x)
                return a.x < b.x;
            return a.y < b.y;
        }
    );

    leaf_table[list_idx].assign(nx * ny, -1);
#pragma omp parallel for
    for (int leaf_idx : cached_leaves_idx) {
        const auto &leaf = getNode(list_idx, leaf_idx);

        int x_start = std::max(0, (int)(leaf.x - leaf.size / 2.f));
        int x_end   = std::min(nx, (int)(leaf.x + leaf.size / 2.f));
        int y_start = std::max(0, (int)(leaf.y - leaf.size / 2.f));
        int y_end   = std::min(ny, (int)(leaf.y + leaf.size / 2.f));

        for (int y = y_start; y < y_end; ++y) {
            for (int x = x_start; x < x_end; ++x) {
                leaf_table[list_idx][IX(x, y)] = leaf_idx;
            }
        }
    }
}

void QTSimulator::cacheNeighbors(int list_idx) {
#pragma omp parallel for
    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(list_idx, leaf_idx);

        std::vector<std::pair<int, int>> neighbors_idx;
        getNeighbors(neighbors_idx, list_idx, leaf_idx);

        for (auto [dir, node_idx] : neighbors_idx) {
            leaf.cached_neighbors_idx[leaf.cached_neighbors_cnt++] = node_idx;
            leaf.neighbor_cnt[dir]++;
        }
    }
}

void QTSimulator::getNeighbors(
    std::vector<std::pair<int, int>> &neighbors, int list_idx, int node_idx
) {
    // 0: left
    // 1: up
    // 2: right
    // 3: down

    auto &node = getNode(list_idx, node_idx);

    float hs = node.size / 2.f;
    float qs = node.size / 4.f;

    float eps = 0.1;

    std::vector<std::pair<float, float>> sample_points = {
        // left
        {node.x - hs - eps, node.y + qs},
        {node.x - hs - eps, node.y - qs},

        // up
        {node.x + qs, node.y - hs - eps},
        {node.x - qs, node.y - hs - eps},

        // right
        {node.x + hs + eps, node.y + qs},
        {node.x + hs + eps, node.y - qs},

        // down
        {node.x + qs, node.y + hs + eps},
        {node.x - qs, node.y + hs + eps},
    };

    for (int i = 0; i < 8; i++) {
        auto [qx, qy] = sample_points[i];
        if (qx < 0.0f || qx >= (float)nx || qy < 0.0f || qy >= (float)ny) {
            continue;
        }

        int neighbor_idx = getNodeIdxAt(qx, qy, list_idx);

        if (neighbor_idx != -1 && neighbor_idx != node_idx) {
            neighbors.push_back({i / 2, neighbor_idx});
        }
    }

    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(
        std::unique(neighbors.begin(), neighbors.end()), neighbors.end()
    );
}

int QTSimulator::getNodeIdxAt(
    float x, float y, int list_idx, int node_idx
) const {

    if (x < 0 || x > nx || y < 0 || y > ny)
        return -1;

    if (!leaf_table[list_idx].empty()) {
        int ix = std::clamp((int)x, 0, nx - 1);
        int iy = std::clamp((int)y, 0, ny - 1);
        return leaf_table[list_idx][IX(ix, iy)];
    }

    auto &node = getNode(list_idx, node_idx);
    if (node.is_leaf)
        return node_idx;

    int child_idx = 0;
    if (x > node.x)
        child_idx |= 1;
    if (y > node.y)
        child_idx |= 2;
    return getNodeIdxAt(x, y, list_idx, node.children_idx[child_idx]);
}

void QTSimulator::getNodesIdxIn(
    float x, float y, float radius_ratio, std::vector<int> &nodes_idx
) {
    nodes_idx.clear();

    int center_idx = getNodeIdxAt(x, y, root_list);
    if (center_idx == -1)
        return;

    auto &center          = getNode(root_list, center_idx);
    float search_radius_2 = center.size * radius_ratio;
    search_radius_2 *= search_radius_2;

    thread_local std::vector<uint8_t> visited;
    if (visited.size() < node_pool[root_list].size())
        visited.resize(node_pool[root_list].size(), false);

    visited[center_idx] = true;
    nodes_idx.push_back(center_idx);

    int read_head = 0;
    while (read_head < nodes_idx.size()) {
        int curr_idx = nodes_idx[read_head++];

        auto &curr          = getNode(root_list, curr_idx);
        auto &neighbors_idx = curr.cached_neighbors_idx;
        int neighbors_cnt   = curr.cached_neighbors_cnt;
        for (int i = 0; i < neighbors_cnt; i++) {
            int node_idx = neighbors_idx[i];
            auto &node   = getNode(root_list, node_idx);

            if (visited[node_idx])
                continue;

            float dist = distance2(node.x, node.y, x, y);
            if (dist <= search_radius_2) {
                visited[node_idx] = true;
                nodes_idx.push_back(node_idx);
            }
        }
    }

    for (int idx : nodes_idx)
        visited[idx] = false;
}

InterpolatedData
QTSimulator::advect(float x, float y, float dt, uint32_t opts) {
    auto current_data = MLSinterpolate(x, y, OPT_VEL_ALL);
    // float u_val = bilerp(u, nx + 1, ny, x, y);
    // float v_val = bilerp(v, nx, ny + 1, x, y);
    float u_val = current_data.u;
    float v_val = current_data.v;

    float past_x = std::clamp(x - u_val * dt, 0.0f, (float)nx);
    float past_y = std::clamp(y - v_val * dt, 0.0f, (float)ny);

    return MLSinterpolate(past_x, past_y, opts);
}

float QTSimulator::getVelocity(std::vector<QuadtreeEdge> &field, int id) {
    if (id == -1)
        return 0.0f;
    return field[id].val;
}

float QTSimulator::nearestParticleDistance(float x, float y) {
    int ix = std::clamp((int)x, 0, nx - 1);
    int iy = std::clamp((int)y, 0, ny - 1);

    float min_dis = std::numeric_limits<float>::infinity();
    for (int j = -1; j <= 1; j++) {
        for (int i = -1; i <= 1; i++) {
            int ni = ix + i, nj = iy + j;

            if (ni < 0 || ni >= nx || nj < 0 || nj >= ny)
                continue;

            int cell_idx = IX(ni, nj);
            for (int idx = phash_head[cell_idx]; idx != -1;
                 idx     = phash_next[idx]) {
                auto &p = particles[idx];

                float d = distance2(p.x, p.y, x, y);
                min_dis = std::min(min_dis, d);
            }
        }
    }
    return std::sqrt(min_dis);
}

float QTSimulator::distance2(float x1, float y1, float x2, float y2) {
    /*
     * Return the square of distance between two point
     */
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

// ######################### MLS #########################
InterpolatedData QTSimulator::MLSinterpolate(float x, float y, uint32_t opts) {
    InterpolatedData result = {0.0f, 1000.f, 0.0f, 0.0f, 0.0f};

    if (opts == OPT_NONE)
        return result;

    thread_local std::vector<int> nodes_to_process;
    getNodesIdxIn(x, y, 3.f, nodes_to_process);
    if (nodes_to_process.empty())
        return result;

    // ==========================================================
    // 1. 細胞中心變數插值 (S, phi)
    // ==========================================================
    if (opts & OPT_CELL_ALL) {
        thread_local std::vector<MLSSamplePoint> cell_points;
        cell_points.clear();

        for (int n_idx : nodes_to_process) {
            MLSMirrorNode(cell_points, n_idx);
        }
        auto tmp_result_cell = solveMLS(x, y, cell_points);
        result.S             = tmp_result_cell.first;
        result.phi           = tmp_result_cell.second;
    }

    // ==========================================================
    // 2. u 速度插值 (垂直邊上的速度)
    // ==========================================================
    if (opts & OPT_U_ALL) {
        thread_local std::vector<MLSSamplePoint> u_samples;
        thread_local std::vector<bool> visited_u_faces;
        u_samples.clear();
        visited_u_faces.assign(QTu.size(), 0);

        for (int n_idx : nodes_to_process) {
            auto &n = getNode(root_list, n_idx);
            // 假設你的 Node 中存有左右垂直邊的面索引 (face indices)
            // 論文提到 T-junction 會自動對齊到大 parent face 取得單一速度值
            int left_face  = n.ul_id;
            int right_face = n.ur_id;
            MLSMirrorEdge(u_samples, visited_u_faces, left_face, QTu, true);
            MLSMirrorEdge(u_samples, visited_u_faces, right_face, QTu, true);
        }
        auto tmp_result_u = solveMLS(x, y, u_samples);
        result.u          = tmp_result_u.first;
        result.u_old      = tmp_result_u.second;
    }

    // ==========================================================
    // 3. v 速度插值 (水平邊上的速度)
    // ==========================================================
    if (opts & OPT_V_ALL) {
        thread_local std::vector<MLSSamplePoint> v_samples;
        thread_local std::vector<bool> visited_v_faces;
        v_samples.clear();
        visited_v_faces.assign(QTv.size(), 0);

        for (int n_idx : nodes_to_process) {
            auto &n         = getNode(root_list, n_idx);
            int top_face    = n.vl_id;
            int bottom_face = n.vr_id;
            MLSMirrorEdge(v_samples, visited_v_faces, top_face, QTv, false);
            MLSMirrorEdge(v_samples, visited_v_faces, bottom_face, QTv, false);
        }
        auto tmp_result_v = solveMLS(x, y, v_samples);
        result.v          = tmp_result_v.first;
        result.v_old      = tmp_result_v.second;
    }

    return result;
}

// 2D 變數的 MLS 求解輔助函式
std::pair<float, float> QTSimulator::solveMLS(
    float x, float y, const std::vector<MLSSamplePoint> &samples
) {
    if (samples.empty())
        return {0, 0};

    Eigen::Matrix3f A = Eigen::Matrix3f::Zero();
    Eigen::Matrix<float, 3, 2> b;
    b.setZero();

    for (const auto &p : samples) {
        float dx = std::abs(p.x - x);
        float dy = std::abs(p.y - y);
        float h  = p.h;

        float eps    = 1e-2f;
        float wx     = std::max(1.f - (dx / h), eps);
        float wy     = std::max(1.f - (dy / h), eps);
        float weight = wx * wy;

        float lx = p.x - x;
        float ly = p.y - y;
        Eigen::Vector3f z_i(lx, ly, 1.f);

        A += weight * (z_i * z_i.transpose());
        b.col(0) += weight * z_i * p.val_1;
        b.col(1) += weight * z_i * p.val_2;
    }

    // 脊迴歸懲罰項，確保 A 在共線或極端情況下仍可逆
    A += Eigen::Matrix3f::Identity() * 1e-6f;

    Eigen::Matrix<float, 3, 2> c = A.ldlt().solve(b);
    return {
        c(2, 0), c(2, 1)
    }; // 局部座標系下，常數項即為 (0,0) 擬合值，對應索引 2
}

void QTSimulator::MLSMirrorNode(
    std::vector<MLSSamplePoint> &sample_points, int node_idx
) {

    float min_x = 0.0f;
    float max_x = nx;
    float min_y = 0.0f;
    float max_y = ny;

    auto &n = getNode(root_list, node_idx);
    sample_points.push_back({n.x, n.y, n.size, n.S, n.phi});

    // 檢查細胞邊界鏡像 (Mirroring)
    float dist_left   = n.x - min_x;
    float dist_right  = max_x - n.x;
    float dist_top    = n.y - min_y;
    float dist_bottom = max_y - n.y;
    float threshold   = n.size * 1.5f;

    bool mirror_x = false;
    float mx      = n.x;
    if (dist_left < threshold) {
        mirror_x = true;
        mx       = min_x - dist_left;
    } else if (dist_right < threshold) {
        mirror_x = true;
        mx       = max_x + dist_right;
    }

    bool mirror_y = false;
    float my      = n.y;
    if (dist_top < threshold) {
        mirror_y = true;
        my       = min_y - dist_top;
    } else if (dist_bottom < threshold) {
        mirror_y = true;
        my       = max_y + dist_bottom;
    }

    if (mirror_x)
        sample_points.push_back({mx, n.y, n.size, n.S, n.phi});
    if (mirror_y)
        sample_points.push_back({n.x, my, n.size, n.S, n.phi});
    if (mirror_x && mirror_y)
        sample_points.push_back({mx, my, n.size, n.S, n.phi});
}

void QTSimulator::MLSMirrorEdge(
    std::vector<MLSSamplePoint> &sample_points,
    std::vector<bool> &visited_faces,
    int face_idx,
    std::vector<QuadtreeEdge> &field,
    bool is_u
) {
    if (face_idx == -1)
        return;
    if (visited_faces[face_idx])
        return;
    visited_faces[face_idx] = 1;

    float min_x = 0.0f;
    float max_x = nx;
    float min_y = 0.0f;
    float max_y = ny;

    auto &face    = field[face_idx];
    float val     = face.val;
    float val_old = face.val_old;
    float face_x  = face.x;
    float face_y  = face.y;
    float face_h  = face.length;
    int sign      = (is_u ? 1 : -1);
    sample_points.push_back({face_x, face_y, face_h, val, val_old});

    // 垂直邊邊界鏡像 (以物理邊界條件處理速度)
    float dist_left   = face_x - min_x;
    float dist_right  = max_x - face_x;
    float dist_top    = face_y - min_y;
    float dist_bottom = max_y - face_y;
    float threshold   = face_h * 1.5f;

    bool mirror_x     = false;
    float mx          = face_x;
    float m_val_x     = val;
    float m_val_old_x = val_old;
    if (dist_left < threshold) {
        mirror_x    = true;
        mx          = min_x - dist_left;
        m_val_x     = -sign * val;     // 固體邊界法向速度反向 (消去穿透流)
        m_val_old_x = -sign * val_old; // 固體邊界法向速度反向 (消去穿透流)
    } else if (dist_right < threshold) {
        mirror_x    = true;
        mx          = max_x + dist_right;
        m_val_x     = -sign * val;
        m_val_old_x = -sign * val_old;
    }

    bool mirror_y     = false;
    float my          = face_y;
    float m_val_y     = val;
    float m_val_old_y = val_old;
    if (dist_top < threshold) {
        mirror_y    = true;
        my          = min_y - dist_top;
        m_val_y     = sign * val;     // 滑動邊界 (Slip BC) 切向速度同向
        m_val_old_y = sign * val_old; // 滑動邊界 (Slip BC) 切向速度同向
    } else if (dist_bottom < threshold) {
        mirror_y    = true;
        my          = max_y + dist_bottom;
        m_val_y     = sign * val;
        m_val_old_y = sign * val_old;
    }

    if (mirror_x)
        sample_points.push_back({mx, face_y, face_h, m_val_x, m_val_old_x});
    if (mirror_y)
        sample_points.push_back({face_x, my, face_h, m_val_y, m_val_old_y});
    if (mirror_x && mirror_y)
        sample_points.push_back({mx, my, face_h, -val, -val_old});
}
