#pragma once

#include "QT/base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

class QTEulerian : public QTSimulatorBase {
  private:
    void recursiveBuildTree(float dt, int node_idx = 0) override;

  public:
    QTEulerian(int width, int height);

    void update(float dt) override; // update every frame
    void reset() override;
};
