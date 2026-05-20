#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <vector>

class MACSimulator {
  private:
    int nx, ny;

    // MAC grid data
    // u for row velocity
    // v for column velocity
    // p for pressure
    // dye for renderer
    std::vector<double> u, v, p, dye;

    // // Eigen solvers
    // Eigen::SparseMatrix<double> A;
    // Eigen::ConjugateGradient<
    //     Eigen::SparseMatrix<double>,
    //     Eigen::Lower | Eigen::Upper>
    //     solver;

    // index calculator
    int IX(int i, int j) const { return i + j * nx; }
    int IX_u(int i, int j) const { return i + j * (nx + 1); }
    int IX_v(int i, int j) const { return i + j * nx; }

    // main simulation functions
    void setBoundaries();
    void advect(float dt);
    void project();
    void applyGravity(float dt);

    // util functions
    double bilerp(
        const std::vector<double> &field, int w, int h, double x, double y
    ) const;

  public:
    MACSimulator(int width, int height);

    void update(float dt); // update every frame
    void addForce(
        float x, float y, float dx, float dy, float radius
    ); // add force and dye
    void initCircle(float cx, float cy, float radis);

    // get functions for renderer
    int getWidth() const { return nx; }
    int getHeight() const { return ny; }
    const std::vector<double> &getDye() const { return dye; }
};
