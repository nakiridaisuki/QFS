#include "QT.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <queue>
#include <utility>
#include <vector>

// PUBLIC
QTSimulator::QTSimulator(int width, int height) : nx(width), ny(height) {
    particle_idx.resize(nx * ny, std::vector<int>(8));
    G     = 150.f;
    Sigma = 0.0f;

    solver.setMaxIterations(80);
    solver.setTolerance(1e-3);
    Eigen::setNbThreads(6);

    root_list = 0;
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, root_list);
    initQuadtree(3, root_list);
}
void QTSimulator::reset() {
    particles.clear();
    node_pool[0].clear();
    node_pool[1].clear();

    root_list = 0;
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, root_list);
    initQuadtree(3, root_list);
}

void QTSimulator::addWater(float x, float y, float radius) {
    recursiveUpdatePhi(x, y, radius, false, root_list);
    smoothing(root_list);
    cacheLeaves(root_list);
    cacheNeighbors(root_list);
}
void QTSimulator::delWater(float x, float y, float radius) {
    float r2 = radius * radius;
    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [&](auto &p) {
            return ((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y)) < r2;
        });
    particles.erase(new_end, particles.end());
    recursiveUpdatePhi(x, y, radius, true, root_list);
    smoothing(root_list);
    cacheLeaves(root_list);
    cacheNeighbors(root_list);
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
    particleSurfaceToGrid();

    for (auto &face : QTu)
        face.val_old = face.val;
    for (auto &face : QTv)
        face.val_old = face.val;

    QTapplyGravity(dt);

    setBoundaries();
    QTproject();
    setBoundaries();

    gridToParticle();
    advectParticles(dt);
    resampleParticles();

    // // Redistance
    updateParticleIdx();
    redistancing();
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

    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        if (std::abs(leaf.phi) < leaf.size * 1.5f) {
            float gamma_phi = 4.f;
            float gamma_u   = 3.f;

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
    for (int iter = 0; iter < iterations; iter++) {
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

    if (std::abs(exp_phi) < node.size && exp_S > (1.f / node.size)) {
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

        int node_idx = 0;
        auto &nodes  = leaf.cached_neighbors_idx;
        for (int dir = 0; dir < 4; dir++) {
            int neighbot_cnt = leaf.neighbor_cnt[dir];

            float x = leaf.x, y = leaf.y, solid_frac = 0;
            int face_sgn;
            std::vector<int> adj_cells_idx = {leaf_idx};
            std::vector<float> coefs;
            auto &first_neighbor = getNode(list_idx, nodes[node_idx]);

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
                adj_cells_idx.push_back(nodes[node_idx]);
                adj_cells_idx.push_back(nodes[node_idx + 1]);

                float coef = face_sgn * 1.f / (1.5f * first_neighbor.size);
                coefs.assign({-1.f * coef, 0.5f * coef, 0.5f * coef});
            } else if (neighbot_cnt == 1) {
                if (first_neighbor.depth == leaf.depth) {
                    if (dir <= 1)
                        goto OVERLAY_NODE;

                    adj_cells_idx.push_back(nodes[node_idx]);
                    float coef = face_sgn * 1.f / leaf.size;
                    coefs.assign({-coef, coef});
                    goto UPDATE;
                }

            OVERLAY_NODE:
                switch (dir) {
                case 0:
                    leaf.ul_id = first_neighbor.ur_id;
                    break;
                case 1:
                    leaf.vl_id = first_neighbor.vr_id;
                    break;
                case 2:
                    leaf.ur_id = first_neighbor.ul_id;
                    break;
                case 3:
                    leaf.vr_id = first_neighbor.vl_id;
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

            node_idx += neighbot_cnt;
        }
    }
}

void QTSimulator::advectQuadtreeDatas(float dt, int list_idx) {
    /*
     * For all leaves.
     * For all edges.
     * advect the Phi, Size function and velocity data
     */

    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf         = getNode(list_idx, leaf_idx);
        auto advected_leaf = advect(leaf.x, leaf.y, dt, OPT_CELL_ALL);
        leaf.phi           = advected_leaf.phi;
        leaf.S             = advected_leaf.S;
    }

    for (auto &face : QTu_new) {
        auto advected_face = advect(face.x, face.y, dt, OPT_U_ALL);
        face.val           = advected_face.u;
    }
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

                for (int idx : particle_idx[IX(ni, nj)]) {
                    auto &p = particles[idx];

                    float dx = std::abs(p.x - face.x);
                    float dy = std::abs(p.y - face.y);

                    float wx = std::max(0.f, 1.f - dx / face.length);
                    float wy = std::max(0.f, 1.f - dy / face.length);
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

    for (int i = 0; i < QTu.size(); i++)
        updateFace(QTu[i], true);

    for (int i = 0; i < QTv.size(); i++)
        updateFace(QTv[i], false);
}

void QTSimulator::particleSurfaceToGrid() {
    /*
     * For all leaves.
     * reconstruct the surface by particles
     */

    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        if (leaf.size < 1.5f && leaf.phi < leaf.size) {
            auto nearest_p = nearestParticle(leaf.x, leaf.y, 2.f);
            if (nearest_p.x != 1e6) {
                float dist = std::sqrt(
                    (nearest_p.x - leaf.x) * (nearest_p.x - leaf.x) +
                    ((nearest_p.y - leaf.y) * (nearest_p.y - leaf.y))
                );

                // if (dist < 0.f)
                leaf.phi = std::min(leaf.phi, dist);
            }
        }
    }
}

void QTSimulator::QTapplyGravity(float dt) {
    /*
     * For all edges.
     * add gravity on horizontial edges
     */

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
    std::vector<Eigen::Triplet<float>> QTtriplets;

    auto addFace = [&](QuadtreeEdge &face) {
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

            div[row] += VA * u_star * face.grad_coeff[i];

            for (int j = 0; j < face.adj_cells_idx.size(); j++) {
                auto &cell_j = getNode(root_list, face.adj_cells_idx[j]);
                if (cell_j.fluid_id == -1)
                    continue;
                int col = cell_j.fluid_id;

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

    for (auto &face : QTu)
        updateFace(face);
    for (auto &face : QTv)
        updateFace(face);
}

void QTSimulator::gridToParticle() {
    /*
     * For all particles
     * get velocity from edges
     */

    float flipRatio = 0.95f;

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

            node.has_particle = true;
            return false;
        });
    particles.erase(new_end, particles.end());

    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        if (leaf.phi < 0 &&
            leaf.size < 1.5f &&
            leaf.phi > -leaf.size &&
            !leaf.has_particle) {

            float px  = leaf.x + (((rand() % 100) / 100.0f) - 0.5f) * leaf.size;
            float py  = leaf.y + (((rand() % 100) / 100.0f) - 0.5f) * leaf.size;
            auto data = MLSinterpolate(px, py, OPT_VEL_ALL);
            particles.push_back({px, py, data.u, data.v});
        }
    }
}

void QTSimulator::updateParticleIdx() {
    /*
     * For all particles
     * update grid hash data
     */

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            particle_idx[IX(i, j)].clear();
        }
    }
    for (int i = 0; i < particles.size(); i++) {
        auto &p = particles[i];
        int ix  = p.x;
        int iy  = p.y;

        ix = std::clamp(ix, 0, nx - 1);
        iy = std::clamp(iy, 0, ny - 1);

        particle_idx[IX(ix, iy)].push_back(i);
    }
}

void QTSimulator::redistancing() {
    /*
     * For all leaves.
     * recompute the phi by particles
     */

    auto surfaceFirst = [&](int a_idx, int b_idx) {
        return std::abs(getNode(root_list, a_idx).phi_new) >
               std::abs(getNode(root_list, b_idx).phi_new);
    };
    std::priority_queue<int, std::vector<int>, decltype(surfaceFirst)> pq(
        surfaceFirst
    );

    for (int leaf_idx : cached_leaves_idx) {
        auto &leaf = getNode(root_list, leaf_idx);
        // if (std::abs(leaf.phi) < 0.5) {
        //     leaf.phi_new = leaf.phi;
        if (std::abs(leaf.phi) < leaf.size * 1.5f) {

            if (leaf.size < 1.5f) {
                auto nearest_p = nearestParticle(leaf.x, leaf.y, 2.f);

                if (nearest_p.x == 1e6 || nearest_p.y == 1e6) {
                    leaf.phi_new = std::numeric_limits<float>::infinity();
                    continue;
                }

                float dist = std::sqrt(
                    (leaf.x - nearest_p.x) * (leaf.x - nearest_p.x) +
                    (leaf.y - nearest_p.y) * (leaf.y - nearest_p.y)
                );
                dist -= particle_radius;

                leaf.phi = std::min(leaf.phi, dist);
            }
            leaf.phi_new = leaf.phi;
            pq.push(leaf_idx);
        } else {
            leaf.phi_new = std::numeric_limits<float>::infinity();
        }
    }

    if (pq.empty()) {
        return;
    }

    while (!pq.empty()) {
        int curr_node_idx = pq.top();
        auto &curr_node   = getNode(root_list, curr_node_idx);
        pq.pop();

        if (curr_node.known)
            continue;
        curr_node.known = true;

        for (int i = 0; i < curr_node.cached_neighbors_cnt; i++) {
            int neighbor_idx = curr_node.cached_neighbors_idx[i];
            auto &neighbor   = getNode(root_list, neighbor_idx);
            float dist       = std::sqrt(
                pow(curr_node.x - neighbor.x, 2) +
                pow(curr_node.y - neighbor.y, 2)
            );

            dist = neighbor.phi < 0 ? -dist : dist;

            float proposed_phi_new = curr_node.phi_new + dist;

            if (std::abs(proposed_phi_new) < std::abs(neighbor.phi_new)) {
                neighbor.phi_new = proposed_phi_new;
                pq.push(neighbor_idx);
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
    {
        auto &node = getNode(list_idx, node_idx);
        new_phi    = circleSDF(cx, cy, radius, node.x, node.y);

        // if (is_delete && node.phi <= 0.f)
        //     node.phi = std::max(node.phi, -new_phi);
        // else if (!is_delete && node.phi > 0.f)
        node.phi = std::min(node.phi, new_phi);

        node_x       = node.x;
        node_y       = node.y;
        node_size    = node.size;
        node_phi     = node.phi;
        node_is_leaf = node.is_leaf;
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

    if (node_size >= 1.5f && std::abs(node_phi) < node_size) {
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
    if (node_size < 1.5f && node_phi < 0) {
        particles.push_back({node_x, node_y});
    }
}

void QTSimulator::collectLeafNodes(
    std::vector<int> &leaves_idx, int list_idx, int node_idx
) const {
    auto &node = getNode(list_idx, node_idx);
    if (node.is_leaf) {
        leaves_idx.push_back(node_idx);
        return;
    }
    for (int i = 0; i < 4; i++) {
        collectLeafNodes(leaves_idx, list_idx, node.children_idx[i]);
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
}

void QTSimulator::cacheNeighbors(int list_idx) {
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

    center.known = true;
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

            if (node.known)
                continue;

            float dist = distance2(node.x, node.y, x, y);
            if (dist <= search_radius_2) {
                node.known = true;
                nodes_idx.push_back(node_idx);
            }
        }
    }

    for (int idx : nodes_idx) {
        getNode(root_list, idx).known = false;
    }
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
                    result  = p;
                }
            }
        }
    }
    return result;
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

    static std::vector<int> nodes_to_process;
    getNodesIdxIn(x, y, 3.f, nodes_to_process);
    if (nodes_to_process.empty())
        return result;

    // ==========================================================
    // 1. 細胞中心變數插值 (S, phi)
    // ==========================================================
    if (opts & OPT_CELL_ALL) {
        static std::vector<MLSSamplePoint> cell_points;
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
        static std::vector<MLSSamplePoint> u_samples;
        static std::vector<int> visited_u_faces;
        u_samples.clear();
        visited_u_faces.clear();

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
        static std::vector<MLSSamplePoint> v_samples;
        static std::vector<int> visited_v_faces;
        v_samples.clear();
        visited_v_faces.clear();

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
    std::vector<int> &visited_faces,
    int face_idx,
    std::vector<QuadtreeEdge> &field,
    bool is_u
) {
    if (face_idx == -1)
        return;
    if (std::find(visited_faces.begin(), visited_faces.end(), face_idx) !=
        visited_faces.end())
        return;
    visited_faces.push_back(face_idx);

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
