#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <vector>

struct Particle {
    float x, y;
    float u, v;
};

class MACSimulator {
  private:
    int nx, ny;
    int fps;
    float G; // Gravity const

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
    void advectParticles(float dt);
    void applyGravity(float dt);
    void markFluidCells();
    void setBoundaries(std::vector<float> &ufield, std::vector<float> &vfield);
    void project();

    // util functions
    float bilerp(
        const std::vector<float> &field, int w, int h, float x, float y
    ) const;
    void bidistri(
        std::vector<float> &field, int w, int h, float x, float y, float value
    );

  public:
    MACSimulator(int width, int height, int fps);

    void update(float dt);                         // update every frame
    void addWater(float x, float y, float radius); // add water and dye
    void delWater(float x, float y, float radius); // add water and dye
    void resampleParticles();

    // get functions for renderer
    int getWidth() const { return nx; }
    int getHeight() const { return ny; }
    float getGravity() { return G; }
    const std::vector<int> &getCell() const { return cell_type; }

    // set functions
    void setGravity(float gravity) { G = gravity; }
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
