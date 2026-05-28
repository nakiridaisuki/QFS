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

    // signed distance to surface, > 0 if inter water
    float phi = std::numeric_limits<float>::infinity(), phi_new;
    // size function
    float S = 0.0f, S_new = 0.0f;
    float pressure = 0.0f;
    float ul = 0.0f, ur = 0.0f, vl = 0.0f, vr = 0.0f;
};

struct InterpolatedData {
    float S, phi, pressure, ul, ur, vl, vr;

    InterpolatedData() = default;
    InterpolatedData(
        float s, float p, float pr, float u_l, float u_r, float v_l, float v_r
    )
        : S(s), phi(p), pressure(pr), ul(u_l), ur(u_r), vl(v_l), vr(v_r) {}

    InterpolatedData(const QuadtreeNode *node) {
        if (node) {
            S = node->S;
            phi = node->phi;
            pressure = node->pressure;
            ul = node->ul;
            ur = node->ur;
            vl = node->vl;
            vr = node->vr;
        } else {
            S = phi = pressure = ul = ur = vl = vr = 0.0f;
        }
    }

    InterpolatedData operator+(const InterpolatedData &x) const {
        return InterpolatedData{
            S + x.S,
            phi + x.phi,
            pressure + x.pressure,
            ul + x.ul,
            ur + x.ur,
            vl + x.vl,
            vr + x.vr
        };
    }

    InterpolatedData operator*(const float x) const {
        return InterpolatedData{
            S * x, phi * x, pressure * x, ul * x, ur * x, vl * x, vr * x
        };
    }
};

struct NeighborData {
    InterpolatedData datas;
    float distance;
};

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
    std::vector<float> u, v, p, u_old, v_old;
    std::vector<float> weight_u, weight_v;
    std::vector<int> cell_type, particles_count, current_count, fluid_map;
    std::vector<Eigen::Triplet<float>> triplets;
    std::vector<std::vector<int>> particle_idx;
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<float>,
        Eigen::Lower | Eigen::Upper>
        solver;

    // index calculator
    int IX(int i, int j) const { return i + j * nx; }
    int IX_u(int i, int j) const { return i + j * (nx + 1); }
    int IX_v(int i, int j) const { return i + j * nx; }

    // main simulation functions
    void particleToGrid();
    void gridToParticle();
    void velExtrapolation();
    void advectParticles(float dt);
    void buildNewTree(float dt);
    void applyGravity(float dt);
    void applySurfaceTension(float dt);
    void markFluidCells();
    void setBoundaries(std::vector<float> &ufield, std::vector<float> &vfield);
    void project();
    void updateParticleIdx();

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

    // util functions
    QuadtreeNode *getNodeAt(QuadtreeNode *node, float x, float y) const;
    void getNodesIn(
        QuadtreeNode *node,
        float x,
        float y,
        float radius,
        std::vector<QuadtreeNode *> &nodes
    );
    void getNeighbors(
        QuadtreeNode *node,
        std::vector<std::pair<int, QuadtreeNode *>> &neighbors
    );
    NeighborData getNeighborData(QuadtreeNode *node, int direction);
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
