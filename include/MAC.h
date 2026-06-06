#pragma once

#include "Base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <algorithm>
#include <vector>

class MACSimulator : public BaseSimulator {
  private:
    int nx, ny;

    // MAC grid data
    // u for row velocity
    // v for column velocity
    std::vector<Particle> particle_place_holder; // just for Renderer
    std::vector<float> u, v, density;
    std::vector<float> u_old, v_old, density_old;
    std::vector<float> weight_u, weight_v;
    std::vector<int> cell_type, current_count, fluid_map;
    std::vector<Eigen::Triplet<float>> triplets;
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<float>,
        Eigen::Lower | Eigen::Upper>
        solver;

    // index calculator
    int IX(int i, int j) const { return i + j * nx; }
    int IX_u(int i, int j) const { return i + j * (nx + 1); }
    int IX_v(int i, int j) const { return i + j * nx; }

    // main simulation functions
    void advectDatas(float dt);
    void applyGravity(float dt);
    void applySurfaceTension(float dt);
    void markFluidCells();
    void setBoundaries(std::vector<float> &ufield, std::vector<float> &vfield);
    void project();
    void velExtrapolation();

    // util functions
    float distance2(float x1, float y1, float x2, float y2);
    float bilerp(
        const std::vector<float> &field, int w, int h, float x, float y
    ) const;
    void bidistri(
        std::vector<float> &field, int w, int h, float x, float y, float value
    );

  public:
    MACSimulator(int width, int height);

    void update(float dt) override; // update every frame
    void addWater(float x, float y, float radius) override;
    void delWater(float x, float y, float radius) override;

    // get functions for renderer
    int getWidth() const override { return nx; }
    int getHeight() const override { return ny; }
    const std::vector<Particle> &getParticles() const override {
        return particle_place_holder;
    }
    std::vector<Line> getLines() const override;
    bool is_water(int x, int y) const override { return cell_type[IX(x, y)]; }

    // set functions
    void reset() override {
        std::fill(u.begin(), u.end(), 0.0);
        std::fill(v.begin(), v.end(), 0.0);
        std::fill(u_old.begin(), u_old.end(), 0.0);
        std::fill(v_old.begin(), v_old.end(), 0.0);
        std::fill(density.begin(), density.end(), 0.0);
        std::fill(density_old.begin(), density_old.end(), 0.0);
        std::fill(cell_type.begin(), cell_type.end(), 0);
    }
};
