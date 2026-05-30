#pragma once

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
    QuadtreeNode *children[4] = {nullptr};

    // signed distance to surface, < 0 if inter water
    float phi = std::numeric_limits<float>::infinity(), phi_new;
    // size function
    float S = 0.0f, S_new = 0.0f;
    float pressure = 0.0f;
    int ul_id = -1, ur_id = -1, vl_id = -1, vr_id = -1;
    std::vector<QuadtreeNode *> cached_neighbors;
    int neighbor_cnt[4] = {0};
};

struct QuadtreeEdge {
    float x, y;
    float length;
    float val = 0;
    float solid_fraction = 0;
    std::vector<QuadtreeNode *> adj_cells;
    std::vector<float> grad_coeff;
    bool advected = false;
};

struct InterpolatedData {
    float S, phi, pressure, u = 0, v = 0;

    // InterpolatedData() = default;
    // InterpolatedData(float s, float p, float pr) : S(s), phi(p), pressure(pr)
    // {}
    //
    // InterpolatedData(const QuadtreeNode *node) {
    //     if (node) {
    //         S = node->S;
    //         phi = node->phi;
    //         pressure = node->pressure;
    //     } else {
    //         S = phi = pressure = 0.0f;
    //     }
    // }
    //
    // InterpolatedData operator+(const InterpolatedData &x) const {
    //     return InterpolatedData{
    //         S + x.S,
    //         phi + x.phi,
    //         pressure + x.pressure,
    //     };
    // }
    //
    // InterpolatedData operator*(const float x) const {
    //     return InterpolatedData{S * x, phi * x, pressure * x};
    // }
};

// struct NeighborData {
//     InterpolatedData datas;
//     float distance;
// };

class QTSimulator : public BaseSimulator {
  private:
    QuadtreeNode *root;
    int nx, ny;
    float G;     // Gravity const
    float Sigma; // surface tension
    float max_u, max_v;

    // QT grid data
    // u for row velocity
    // v for column velocity
    // p for pressure
    std::vector<Particle> particles;
    std::vector<int> cell_type, particles_count, current_count, fluid_map;
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<float>,
        Eigen::Lower | Eigen::Upper>
        solver;

    std::vector<QuadtreeEdge> QTu, QTv, QTu_new, QTv_new;

    void updateParticleIdx();

    // Quad Tree simulation functions
    void QTapplyGravity(float dt);
    void QTproject();

    // Quad Tree functions
    void initQuadtree(QuadtreeNode *node, int max_depth);
    void subdivideNode(QuadtreeNode *node);
    void recursiveGetLines(QuadtreeNode *node, std::vector<Line> &lines) const;
    void
    recursiveUpdatePhi(QuadtreeNode *node, float cx, float cy, float radius);
    void recursiveFree(QuadtreeNode *node);
    void recursiveBuildTree(QuadtreeNode *node, float dt);
    void advectQuadtreeDatas(QuadtreeNode *node, float dt);
    void computeSizingFunction(QuadtreeNode *node, float dt);
    void propagateSizingFunction();
    void redistancing(QuadtreeNode *root);
    void smoothing(QuadtreeNode *new_root);

    void findAllEdges(QuadtreeNode *new_root);

    // util functions
    void cacheNeighbors(QuadtreeNode *new_root);
    float getVelocity(std::vector<QuadtreeEdge> &field, int id);
    InterpolatedData advect(float x, float y, float dt);
    QuadtreeNode *getNodeAt(QuadtreeNode *node, float x, float y) const;
    void getNodesIn(
        QuadtreeNode *node,
        float x,
        float y,
        float radius,
        std::vector<QuadtreeNode *> &nodes
    );
    void getNeighbors(
        QuadtreeNode *root,
        QuadtreeNode *node,
        std::vector<std::pair<int, QuadtreeNode *>> &neighbors
    );
    // NeighborData getNeighborData(QuadtreeNode *node, int direction);
    Particle nearestParticle(float x, float y, float radius);
    void commitQuadtreeS(QuadtreeNode *node);
    void commitQuadtreePhi(QuadtreeNode *node);
    void collectLeafNodes(
        QuadtreeNode *node, std::vector<QuadtreeNode *> &leaves
    ) const;
    float circleSDF(float cx, float cy, float radius, float x, float y) {
        return std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy)) - radius;
    };
    InterpolatedData MLSinterpolate(float x, float y);
    float bilerp(
        const std::vector<float> &field, int w, int h, float x, float y
    ) const;
    void bidistri(
        std::vector<float> &field, int w, int h, float x, float y, float value
    );

  public:
    QTSimulator(int width, int height);

    void update(float dt);                         // update every frame
    void addWater(float x, float y, float radius); // add water and dye
    void delWater(float x, float y, float radius); // add water and dye
    void resampleParticles();

    // get functions for renderer
    int getWidth() const override { return nx; }
    int getHeight() const override { return ny; }
    const std::vector<int> &getCell() const override { return cell_type; }
    const QuadtreeNode *getNodeAt(float x, float y) const {
        return getNodeAt(root, x, y);
    };
    const std::vector<Particle> &getParticles() const override {
        return particles;
    }
    std::vector<Line> getLines() const override;

    // get functions for main loop
    float getGravity() { return G; }
    float getSurfaceTension() { return Sigma; }
    float getMaxVel() { return std::max(max_u, max_v); }

    // set functions
    void setGravity(float gravity) { G = gravity; }
    void setSigma(float sigma) { Sigma = sigma; }
    void reset();
};
