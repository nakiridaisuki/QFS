#include "QT/Eulerian.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <omp.h>
#include <vector>

// PUBLIC
QTEulerian::QTEulerian(int width, int height) : QTSimulatorBase(width, height) {
    simulator_type = SimType::QT_Eulerian;
    root_list      = 0;
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, root_list);
    initQuadtree(3, root_list);
}
void QTEulerian::reset() {
    QTSimulatorBase::reset();

    root_list = 0;
    allocate((float)nx / 2, (float)ny / 2, (float)nx, 0, root_list);
    initQuadtree(3, root_list);
}

void QTEulerian::update(float dt) {

    buildNewTree(dt);

    // Start solve new frame data
    // 1. apply gravity on new tree
    applyGravity(dt);
    applySurfaceTension(dt);
    // 2. solve the pressure
    setBoundaries();
    project();
    setBoundaries();
    // 3. extrapolate velocity
    velExtrapolation();
    // 4. rebuild level set(phi)
    redistancing();

    // Reset interact state
    add_water = false;
    del_water = false;
}

void QTEulerian::recursiveBuildTree(float dt, int node_idx) {
    /*
     * Build new quadtree -> new tree list id = 1 - root_list
     */

    int new_list = 1 - root_list;
    auto &node   = getNode(new_list, node_idx);

    if (node.size < 1.5f)
        return;

    auto advected_data = advect(node.x, node.y, dt, OPT_PHI);
    auto data          = MLSinterpolate(node.x, node.y, OPT_S);

    float exp_phi = advected_data.phi;
    float exp_S   = data.S;

    if (add_water || del_water) {
        float water_phi = circleSDF(node.x, node.y);

        if (add_water && water_phi < exp_phi)
            node.is_new = true;
        if (del_water && -water_phi > exp_phi)
            node.is_new = true;

        if (add_water)
            exp_phi = std::min(exp_phi, water_phi);
        if (del_water)
            exp_phi = std::max(exp_phi, -water_phi);
        exp_S = std::max(exp_S, 1.f / node.size + 1);
    }

    if (std::abs(exp_phi) < node.size && exp_S > (1.f / node.size)) {
        subdivideNode(new_list, node_idx);

        for (int i = 0; i < 4; i++) {
            recursiveBuildTree(dt, getNode(new_list, node_idx).children_idx[i]);
        }
    }
}
