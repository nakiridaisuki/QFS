#pragma once

#include "MAC/base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <algorithm>
#include <vector>

class MACEulerian : public MACSimulatorBase {
  private:
    std::vector<float> density, density_old;

    // main simulation functions
    void advectDatas(float dt);
    void applySurfaceTension(float dt) override;
    void markFluidCells() override;

  public:
    MACEulerian(int width, int height);

    void update(float dt) override; // update every frame
    void addWater(float x, float y, float radius) override;
    void delWater(float x, float y, float radius) override;

    // set functions
    void reset() override {
        MACSimulatorBase::reset();
        std::fill(density.begin(), density.end(), 0.0);
        std::fill(density_old.begin(), density_old.end(), 0.0);
    }
};
