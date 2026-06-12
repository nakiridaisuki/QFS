#include "QT/EXNBFLIP.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <math.h>
#include <omp.h>
#include <random>
#include <vector>

// PUBLIC
QTEXNBFLIP::QTEXNBFLIP(int width, int height) : QTSimulatorBase(width, height) {
    simulator_type = SimType::QT_EXNBFLIP;
    phash_head.resize(nx * ny, -1);
    phash_next.clear();

    root_list = 0;
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, root_list);
    initQuadtree(3, root_list);
}
void QTEXNBFLIP::reset() {
    particles.clear();
    phash_head.assign(nx * ny, -1);
    phash_next.clear();

    root_list = 0;
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, root_list);
    initQuadtree(3, root_list);
}

void QTEXNBFLIP::update(float dt) {

    buildNewTree(dt);

    particleToGrid();

    for (auto &face : QTu)
        face.val_old = face.val;
    for (auto &face : QTv)
        face.val_old = face.val;

    applyGravity(dt);
    applySurfaceTension(dt);

    setBoundaries();
    project();
    setBoundaries();

    velExtrapolation();

    gridToParticle();
    advectParticles(dt);

    redistancing();
    resampleParticles();

    updateParticleIdx();
    // reconstructSurface();
    // redistancing();

    // Reset interact state
    add_water = false;
    del_water = false;
}

void QTEXNBFLIP::recursiveBuildTree(float dt, int node_idx) {
    /*
     * Recursive
     * build the tree until size to 1.f
     */

    int new_list = 1 - root_list;
    auto &node   = getNode(new_list, node_idx);

    auto advected_data = advect(node.x, node.y, dt, OPT_PHI);
    auto data          = MLSinterpolate(node.x, node.y, OPT_S);

    float exp_phi = advected_data.phi;
    float exp_S   = data.S;

    if (add_water || del_water) {
        float water_phi = circleSDF(node.x, node.y);

        if (add_water && water_phi < exp_phi)
            node.is_new = true;
        if (del_water && -water_phi > exp_phi)
            node.is_new = true;

        if (add_water)
            exp_phi = std::min(exp_phi, water_phi);
        if (del_water)
            exp_phi = std::max(exp_phi, -water_phi);
        exp_S = std::max(exp_S, 1.f / node.size + 1);
    }

    if (node.size < 1.5f) {
        if (add_water && node.is_new && exp_phi <= 0.f) {
            particles.push_back({node.x + 0.25f, node.y + 0.25f});
            particles.push_back({node.x + 0.25f, node.y - 0.25f});
            particles.push_back({node.x - 0.25f, node.y + 0.25f});
            particles.push_back({node.x - 0.25f, node.y - 0.25f});
        }
        return;
    }

    if (std::abs(exp_phi) < node.size && exp_S > (1.f / node.size)) {
        subdivideNode(new_list, node_idx);

        for (int i = 0; i < 4; i++) {
            recursiveBuildTree(dt, getNode(new_list, node_idx).children_idx[i]);
        }
    }
}

void QTEXNBFLIP::particleToGrid() {
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

void QTEXNBFLIP::gridToParticle() {
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

void QTEXNBFLIP::advectParticles(float dt) {
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

void QTEXNBFLIP::resampleParticles() {
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

void QTEXNBFLIP::updateParticleIdx() {
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

void QTEXNBFLIP::reconstructSurface() {
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

float QTEXNBFLIP::nearestParticleDistance(float x, float y) {
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
