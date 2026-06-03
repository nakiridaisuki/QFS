#pragma once

#include "Base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <cstdint>
#include <limits>
#include <vector>

struct QuadtreeNode {
    // x, y: position
    // size: edge length of node
    float x, y, size;
    int depth;

    bool is_leaf      = true;
    bool known        = false;
    bool has_particle = false;
    int children_idx[4];

    // signed distance to surface, < 0 if inter water
    float phi = std::numeric_limits<float>::infinity(), phi_new;
    // size function
    float S = 0.0f, S_new = 0.0f;
    int ul_id = -1, ur_id = -1, vl_id = -1, vr_id = -1;
    int cached_neighbors_idx[8];
    int cached_neighbors_cnt = 0;
    int neighbor_cnt[4]      = {0};
    int fluid_id             = -1;

    bool operator<(const QuadtreeNode &other) const {
        return depth < other.depth;
    }
};

struct QuadtreeEdge {
    float x, y;
    float length;
    float val            = 0;
    float solid_fraction = 0;
    std::vector<int> adj_cells_idx;
    std::vector<float> grad_coeff;
    float val_old = 0;
};

struct InterpolatedData {
    float S = 0, phi = 0;
    float u = 0, v = 0;
    float u_old = 0, v_old = 0;
};

struct MLSSamplePoint {
    float x, y;         // 採樣點的物理座標（細胞中心或邊中心）
    float h;            // 該點對應的網格尺寸
    float val_1, val_2; // 該點的數值
};

enum InterpOptions : uint32_t {
    OPT_NONE  = 0,
    OPT_S     = 1 << 0,
    OPT_PHI   = 1 << 1,
    OPT_U     = 1 << 2,
    OPT_V     = 1 << 3,
    OPT_U_OLD = 1 << 4,
    OPT_V_OLD = 1 << 5,

    OPT_CELL_ALL = OPT_S | OPT_PHI,
    OPT_U_ALL    = OPT_U_OLD | OPT_U,
    OPT_V_ALL    = OPT_V_OLD | OPT_V,
    OPT_VEL_ALL  = OPT_U_ALL | OPT_V_ALL,
    OPT_ALL      = OPT_CELL_ALL | OPT_VEL_ALL
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
    std::vector<Particle> particles;
    std::vector<std::vector<int>> particle_idx;
    std::vector<int> cached_leaves_idx;
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<float>,
        Eigen::Lower | Eigen::Upper>
        solver;
    std::vector<QuadtreeEdge> QTu, QTv, QTu_new, QTv_new;

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
    void updateParticleIdx();

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
        float cx,
        float cy,
        float radius,
        bool is_delete,
        int list_idx,
        int node_idx = 0
    );
    void recursiveBuildTree(float dt, int list_idx, int node_idx = 0);
    void cacheNeighbors(int list_idx);
    void cacheLeaves(int list_idx);
    int getNodeIdxAt(float x, float y, int list_idx, int node_idx = 0) const;
    void getNodesIdxIn(
        float x, float y, float radius_ratio, std::vector<int> &nodes_idx
    );
    void getNeighbors(
        std::vector<std::pair<int, int>> &neighbors, int list_idx, int node_idx
    );
    void collectLeafNodes(
        std::vector<int> &leaves_idx, int list_idx, int node_idx = 0
    ) const;

    // util functions
    float distance2(float x1, float y1, float x2, float y2);
    float getVelocity(std::vector<QuadtreeEdge> &field, int id);
    InterpolatedData
    advect(float x, float y, float dt, uint32_t opts = OPT_ALL);
    Particle nearestParticle(float x, float y, float radius);
    int IX(int i, int j) const { return i + j * nx; }
    float circleSDF(float cx, float cy, float radius, float x, float y) {
        return std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy)) - radius;
    };
    InterpolatedData MLSinterpolate(float x, float y, uint32_t opts = OPT_ALL);
    void
    MLSMirrorNode(std::vector<MLSSamplePoint> &sample_points, int node_idx);
    void MLSMirrorEdge(
        std::vector<MLSSamplePoint> &sample_points,
        std::vector<int> &visited_faces,
        int face_idx,
        std::vector<QuadtreeEdge> &field,
        bool is_u
    );
    std::pair<float, float>
    solveMLS(float x, float y, const std::vector<MLSSamplePoint> &samples);

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
