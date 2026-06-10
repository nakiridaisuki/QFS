#pragma once

#include "Base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <vector>

class MACSimulator : public BaseSimulator {
  private:
    int nx, ny;
    float G;     // Gravity const
    float Sigma; // surface tension
    float max_u, max_v;

    // MAC grid data
    // u for row velocity
    // v for column velocity
    // p for pressure
    std::vector<Particle> particles;
    std::vector<float> u, v, p, u_old, v_old;
    std::vector<float> weight_u, weight_v;
    std::vector<int> cell_type, particles_count, current_count, fluid_map;
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
    void particleToGrid();
    void gridToParticle();
    void velExtrapolation();
    void advectParticles(float dt);
    void applyGravity(float dt);
    void applySurfaceTension(float dt);
    void markFluidCells();
    void setBoundaries(std::vector<float> &ufield, std::vector<float> &vfield);
    void project();
    void resampleParticles();

    // util functions
    float bilerp(
        const std::vector<float> &field, int w, int h, float x, float y
    ) const;
    void bidistri(
        std::vector<float> &field, int w, int h, float x, float y, float value
    );

  public:
    MACSimulator(int width, int height);

    void update(float dt);                         // update every frame
    void addWater(float x, float y, float radius); // add water and dye
    void delWater(float x, float y, float radius); // add water and dye

    // get functions for renderer
    int getWidth() const override { return nx; }
    int getHeight() const override { return ny; }
    const std::vector<Particle> &getParticles() const override {
        return particles;
    }
    std::vector<Line> getLines() const override;
    bool is_water(int x, int y) const override { return cell_type[IX(x, y)]; }

    // get functions for main loop
    float getGravity() { return G; }
    float getSurfaceTension() { return Sigma; }
    float getMaxVel() { return std::max(max_u, max_v); }

    // set functions
    void setGravity(float gravity) { G = gravity; }
    void setSigma(float sigma) { Sigma = sigma; }
    void reset() {
        particles.clear();
        std::fill(u.begin(), u.end(), 0.0);
        std::fill(v.begin(), v.end(), 0.0);
        std::fill(u_old.begin(), u_old.end(), 0.0);
        std::fill(v_old.begin(), v_old.end(), 0.0);
        std::fill(p.begin(), p.end(), 0.0);
        std::fill(cell_type.begin(), cell_type.end(), 0);
    }
};
