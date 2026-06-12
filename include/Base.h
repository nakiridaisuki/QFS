#pragma once

#include <vector>

namespace SimType {
enum SimulatorType {
    MAC_Eulerian = 0,
    MAC_FLIP     = 1,
    QT_Eulerian  = 2,
    QT_EXNBFLIP  = 3
};
}

struct Particle {
    float x, y;
    float u, v;
};

struct Line {
    float x1, y1, x2, y2;
};

class BaseSimulator {
  protected:
    int nx, ny;
    int PCGItertimes = 0;
    float G;     // Gravity const
    float Sigma; // surface tension
    int simulator_type;

  public:
    BaseSimulator(int width, int height) : nx(width), ny(height) {}
    virtual ~BaseSimulator() = default;
    virtual const std::vector<Particle> *getParticles() const {
        return nullptr;
    }
    virtual std::vector<Line> getLines() const = 0;
    virtual bool is_water(int x, int y) const  = 0;
    virtual void update(float dt)              = 0;

    // get functions for main loop
    float getGravity() const { return G; }
    float getSurfaceTension() const { return Sigma; }
    int getPCGIter() const { return PCGItertimes; }
    int getWidth() const { return nx; }
    int getHeight() const { return ny; }
    int getSimType() const { return simulator_type; }

    void setGravity(float gravity) { G = gravity; }
    void setSigma(float sigma) { Sigma = sigma; }

    virtual void
    addWater(float x, float y, float radius) = 0; // add water and dye
    virtual void
    delWater(float x, float y, float radius) = 0; // delete water and dye

    virtual void reset() = 0;
};
