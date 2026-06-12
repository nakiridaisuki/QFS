#pragma once

#include "QT/base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <vector>

class QTEXNBFLIP : public QTSimulatorBase {
  private:
    std::vector<Particle> particles;
    std::vector<int> phash_head, phash_next;

    void recursiveBuildTree(float dt, int node_idx = 0) override;
    void advectParticles(float dt);
    void particleToGrid();
    void gridToParticle();
    void resampleParticles();
    void reconstructSurface();

    // util functions
    float nearestParticleDistance(float x, float y);
    void updateParticleIdx();

  public:
    QTEXNBFLIP(int width, int height);

    void update(float dt) override; // update every frame

    // get functions for renderer
    const std::vector<Particle> *getParticles() const override {
        return &particles;
    }

    // set functions
    void reset() override;
};
