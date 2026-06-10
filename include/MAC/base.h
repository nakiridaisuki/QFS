#pragma once

#include "Base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <algorithm>
#include <vector>

class MACSimulatorBase : public BaseSimulator {
  protected:
    // MAC grid data
    // u for row velocity
    // v for column velocity
    std::vector<float> u, v;
    std::vector<float> u_old, v_old;
    std::vector<float> weight_u, weight_v;
    std::vector<int> cell_type, current_count, fluid_map;
    std::vector<Eigen::Triplet<float>> triplets;
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<float>,
        Eigen::Lower | Eigen::Upper>
        solver;

    // some local variables
    std::vector<int> valid_u, valid_v;
    std::vector<float> phi[2];
    std::vector<float> nx_n, ny_n;
    std::vector<float> kappa;

    // main simulation functions
    virtual void applyGravity(float dt) final;
    virtual void
    setBoundaries(std::vector<float> &ufield, std::vector<float> &vfield) final;
    virtual void velExtrapolation() final;
    virtual void project() final;
    virtual void applySurfaceTension(float dt) = 0;
    virtual void markFluidCells()              = 0;

    // index calculator
    int IX(int i, int j) const { return i + j * nx; }
    int IX_u(int i, int j) const { return i + j * (nx + 1); }
    int IX_v(int i, int j) const { return i + j * nx; }

    // util functions
    float distance2(float x1, float y1, float x2, float y2) {
        return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    }
    float bilerp(
        const std::vector<float> &field, int w, int h, float x, float y
    ) const;
    void bidistri(
        std::vector<float> &field, int w, int h, float x, float y, float value
    );

  public:
    MACSimulatorBase(int width, int height);

    // get functions for renderer
    std::vector<Line> getLines() const override final;
    bool is_water(int x, int y) const override final {
        return cell_type[IX(x, y)];
    }

    // set functions
    void reset() override {
        std::fill(u.begin(), u.end(), 0.0);
        std::fill(v.begin(), v.end(), 0.0);
        std::fill(u_old.begin(), u_old.end(), 0.0);
        std::fill(v_old.begin(), v_old.end(), 0.0);
        std::fill(cell_type.begin(), cell_type.end(), 0);
    }
};
