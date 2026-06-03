#pragma once

#include <iostream>

#include "Base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <limits>
#include <vector>

struct QuadtreeNode {
    // x, y: position
    // size: edge length of node
    float x, y, size;
    int depth;

    bool is_leaf = true;
    bool known = false;
    bool has_particle = false;
    int children_idx[4];

    // signed distance to surface, < 0 if inter water
    float phi = std::numeric_limits<float>::infinity(), phi_new;
    // size function
    float S = 0.0f, S_new = 0.0f;
    float pressure = 0.0f;
    int ul_id = -1, ur_id = -1, vl_id = -1, vr_id = -1;
    int cached_neighbors_idx[8];
    int cached_neighbors_cnt = 0;
    int neighbor_cnt[4] = {0};
    int fluid_id = -1;

    bool operator<(const QuadtreeNode &other) const {
        return depth < other.depth;
    }
};

struct QTNSurfaceFirst {
    bool operator()(QuadtreeNode *a, QuadtreeNode *b) {
        return std::abs(a->phi_new) > std::abs(b->phi_new);
    }
};
struct QTNSmallNodeFirst {
    bool operator()(QuadtreeNode *a, QuadtreeNode *b) {
        return a->depth < b->depth;
    }
};

struct QuadtreeEdge {
    float x, y;
    float length;
    float val = 0;
    float solid_fraction = 0;
    std::vector<int> adj_cells_idx;
    std::vector<float> grad_coeff;
    float val_old = 0;
};

struct InterpolatedData {
    float S, phi, pressure;
    float u = 0, v = 0;
    float u_old = 0, v_old = 0;
};

class QTSimulator : public BaseSimulator {
  private:
    // QuadtreeNode *root;
    int root_list;
    std::vector<QuadtreeNode> node_pool[2];
    int nx, ny;
    float G;     // Gravity const
    float Sigma; // surface tension
    float max_u, max_v;
    float particle_radius = 0.5f;

    // QT grid data
    // u for row velocity
    // v for column velocity
    // p for pressure
    std::vector<Particle> particles;
    std::vector<std::vector<int>> particle_idx;
    std::vector<int> cached_leaves_idx;
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<float>,
        Eigen::Lower | Eigen::Upper>
        solver;

    std::vector<QuadtreeEdge> QTu, QTv, QTu_new, QTv_new;

    void updateParticleIdx();

    // Quad Tree simulation functions
    void computeSizingFunction(float dt);
    void propagateSizingFunction();
    void QTapplyGravity(float dt);
    void QTproject();
    void advectParticles(float dt);
    void particleToGrid();
    void particleSurfaceToGrid();
    void gridToParticle();
    void resampleParticles();
    void redistancing();
    void findAllEdges(int list_idx);
    void advectQuadtreeDatas(float dt, int list_idx);
    void setBoundaries();

    // Quad Tree functions
    int allocate(float x, float y, float size, int depth, int list_idx) {
        int idx = node_pool[list_idx].size();
        node_pool[list_idx].push_back({x, y, size, depth});
        return idx;
    }
    QuadtreeNode &getNode(int list_idx, int node_idx) {
        return node_pool[list_idx][node_idx];
    }
    const QuadtreeNode &getNode(int list_idx, int node_idx = 0) const {
        return node_pool[list_idx][node_idx];
    }
    void initQuadtree(int max_depth, int list_idx, int node_idx = 0);
    void subdivideNode(int list_idx, int node_idx = 0);
    void smoothing(int list_idx);
    void recursiveGetLines(
        std::vector<Line> &lines, int list_idx, int node_idx = 0
    ) const;
    void recursiveUpdatePhi(
        float cx, float cy, float radius, int list_idx, int node_idx = 0
    );
    void recursiveBuildTree(float dt, int list_idx, int node_idx = 0);
    void cacheNeighbors(int list_idx);
    void cacheLeaves(int list_idx);
    int getNodeIdxAt(float x, float y, int list_idx, int node_idx = 0) const;
    // void getNodesIn(
    //     QuadtreeNode *node,
    //     float x,
    //     float y,
    //     float radius,
    //     std::vector<QuadtreeNode *> &nodes
    // );
    void getNeighbors(
        std::vector<std::pair<int, int>> &neighbors, int list_idx, int node_idx
    );
    void collectLeafNodes(
        std::vector<int> &leaves_idx, int list_idx, int node_idx = 0
    ) const;

    // util functions
    int IX(int i, int j) const { return i + j * nx; }
    float getVelocity(std::vector<QuadtreeEdge> &field, int id);
    InterpolatedData advect(float x, float y, float dt);
    Particle nearestParticle(float x, float y, float radius);
    void commitQuadtreeS();
    void commitQuadtreePhi();
    float circleSDF(float cx, float cy, float radius, float x, float y) {
        return std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy)) - radius;
    };
    InterpolatedData MLSinterpolate(float x, float y);

  public:
    QTSimulator(int width, int height);

    void update(float dt);                         // update every frame
    void addWater(float x, float y, float radius); // add water and dye
    void delWater(float x, float y, float radius); // add water and dye

    // get functions for renderer
    int getWidth() const override { return nx; }
    int getHeight() const override { return ny; }
    const QuadtreeNode &getNodeAt(float x, float y) const {
        int idx = getNodeIdxAt(x, y, root_list);
        return getNode(root_list, idx);
    };
    const std::vector<Particle> &getParticles() const override {
        return particles;
    }
    std::vector<Line> getLines() const override;
    bool is_water(int x, int y) const override {
        // std::cout << node_pool[root_list].size() << std::endl;
        int idx = getNodeIdxAt(x + 0.5, y + 0.5, root_list);
        return getNode(root_list, idx).phi <= 0;
    }

    // get functions for main loop
    float getGravity() { return G; }
    float getSurfaceTension() { return Sigma; }
    float getMaxVel() { return std::max(max_u, max_v); }

    // set functions
    void setGravity(float gravity) { G = gravity; }
    void setSigma(float sigma) { Sigma = sigma; }
    void reset();
};
