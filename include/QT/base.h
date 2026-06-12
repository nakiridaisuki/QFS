#pragma once

#include "Base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

struct QuadtreeNode {
    // x, y: position
    // size: edge length of node
    float x, y, size;
    int depth;

    bool is_leaf = true;
    bool is_new  = false;
    bool known   = false;
    int children_idx[4];

    // signed distance to surface, < 0 if inter water
    float phi     = std::numeric_limits<float>::infinity();
    float phi_new = std::numeric_limits<float>::infinity();
    // size function
    float S = 0.0f, S_new = 0.0f;
    int ul_id = -1, ur_id = -1, vl_id = -1, vr_id = -1;
    int cached_neighbors_idx[8];
    int cached_neighbors_cnt = 0;
    int neighbor_cnt[4]      = {0};
    int fluid_id             = -1;
    int particle_cnt         = 0;
};

struct QuadtreeEdge {
    float x, y;
    float length;
    float val            = 0;
    float solid_fraction = 0;
    std::vector<int> adj_cells_idx;
    std::vector<float> grad_coeff;
    bool is_new   = false;
    float val_old = 0;
};

struct InterpolatedData {
    float S = 0, phi = 0;
    float u = 0, v = 0;
    float u_old = 0, v_old = 0;
};

struct MLSSamplePoint {
    float x, y;
    float h;
    float val_1, val_2;
};

struct NeighborDatas {
    float distance[4];
    float phi[4];
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

class QTSimulatorBase : public BaseSimulator {
  protected:
    // Simulator datas
    int root_list;
    std::vector<QuadtreeNode> node_pool[2];

    // Interact datas
    bool add_water = false;
    bool del_water = false;
    float water_x, water_y, water_radius;

    // Quad Tree datas
    // u for row velocity
    // v for column velocity
    std::vector<int> cached_leaves_idx, leaf_table[2];
    std::vector<Eigen::Triplet<float>> QTtriplets;
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<float>,
        Eigen::Lower | Eigen::Upper>
        solver;
    std::vector<QuadtreeEdge> QTu, QTv, QTu_new, QTv_new;

    // Some local variables
    std::vector<int> leaf_fluid_id_offset;
    std::vector<int> leaf_offset, node_offset;
    std::vector<std::vector<Eigen::Triplet<float>>> thread_triplets;

    // Quad Tree simulation pipeline functions
    void buildNewTree(float dt);
    void applyGravity(float dt);
    void applySurfaceTension(float dt);
    void project();
    void setBoundaries();
    void velExtrapolation();
    void redistancing();
    void FMMSolver(std::vector<std::pair<float, int>> &init_datas);

    // Build New Tree
    void computeSizingFunction(float dt);
    void propagateSizingFunction();
    virtual void recursiveBuildTree(float dt, int node_idx = 0) = 0;
    void smoothing();
    void cacheNeighbors();
    void cacheLeaves();
    void findAllEdges();
    void advectQuadtreeDatas(float dt);

    // Quad Tree manipulate functions
    void initQuadtree(int max_depth, int list_idx, int node_idx = 0);
    int allocate(float x, float y, float size, int depth, int list_idx) {
        int idx = node_pool[list_idx].size();
        node_pool[list_idx].push_back({x, y, size, depth});
        return idx;
    }
    void subdivideNode(int list_idx, int node_idx = 0);
    QuadtreeNode &getNode(int list_idx, int node_idx) {
        return node_pool[list_idx][node_idx];
    }
    const QuadtreeNode &getNode(int list_idx, int node_idx = 0) const {
        return node_pool[list_idx][node_idx];
    }
    int getNodeIdxAt(float x, float y, int list_idx, int node_idx = 0) const;
    void getNodesIdxIn(
        float x, float y, float radius_ratio, std::vector<int> &nodes_idx
    );
    void collectLeafNodes(std::vector<int> &leaves_idx, int list_idx);
    void getNeighbors(
        std::vector<std::pair<int, int>> &neighbors, int list_idx, int node_idx
    );
    NeighborDatas getNeighborDatas(int node_idx);

    // util functions
    inline int IX(int i, int j) const { return i + j * nx; }
    inline float distance2(float x1, float y1, float x2, float y2) {
        return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    }
    inline float getVelocity(std::vector<QuadtreeEdge> &field, int id) {
        return (id == -1 ? 0.0f : field[id].val);
    }
    inline float circleSDF(float x, float y) {
        return std::sqrt(distance2(x, y, water_x, water_y)) - water_radius;
    };
    float computeCurvature(float x, float y, float size);
    float computeWk(QuadtreeEdge &face);

    // MLS interpolate functions
    InterpolatedData
    advect(float x, float y, float dt, uint32_t opts = OPT_ALL);
    InterpolatedData MLSinterpolate(float x, float y, uint32_t opts = OPT_ALL);
    std::pair<float, float>
    solveMLS(float x, float y, const std::vector<MLSSamplePoint> &samples);
    void
    MLSMirrorNode(std::vector<MLSSamplePoint> &sample_points, int node_idx);
    void MLSMirrorEdge(
        std::vector<MLSSamplePoint> &sample_points,
        std::vector<uint64_t> &visited_faces,
        uint64_t visit_epoch,
        int face_idx,
        std::vector<QuadtreeEdge> &field,
        bool is_u
    );

    // Renderer data/functions
    void recursiveGetLines(
        std::vector<Line> &lines, int list_idx, int node_idx = 0
    ) const;

  public:
    QTSimulatorBase(int width, int height);

    void addWater(float x, float y, float radius) override;
    void delWater(float x, float y, float radius) override;

    // get functions for renderer
    const std::vector<QuadtreeEdge> &getUVs(bool is_u = false) const {
        if (is_u)
            return QTu;
        else
            return QTv;
    }
    const QuadtreeNode &getNodeAt(float x, float y) const {
        int idx = getNodeIdxAt(x, y, root_list);
        return getNode(root_list, idx);
    };
    const std::vector<Particle> *getParticles() const override {
        return nullptr;
    }
    std::vector<Line> getLines() const override;
    bool is_water(int x, int y) const override {
        int idx = getNodeIdxAt(x + 0.5, y + 0.5, root_list);
        return getNode(root_list, idx).phi <= 0;
    }

    // set functions
    void reset() override;
};
