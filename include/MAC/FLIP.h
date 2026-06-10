#pragma once

#include "MAC/base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <vector>

class MACFLIP : public MACSimulatorBase {
  private:
    std::vector<Particle> particles;
    std::vector<int> particles_count;

    // main simulation functions
    void particleToGrid();
    void gridToParticle();
    void advectParticles(float dt);
    void applySurfaceTension(float dt) override;
    void markFluidCells() override;
    void resampleParticles();

  public:
    MACFLIP(int width, int height);

    void update(float dt) override; // update every frame
    void addWater(float x, float y, float radius) override; // add water and dye
    void delWater(float x, float y, float radius) override; // add water and dye

    const std::vector<Particle> &getParticles() const override {
        return particles;
    }

    // set functions
    void reset() override { particles.clear(); }
};
