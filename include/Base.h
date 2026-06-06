#pragma once

#include <vector>

struct Particle {
    float x, y;
    float u, v;
};

struct Line {
    float x1, y1, x2, y2;
};

class BaseSimulator {
  protected:
    float G;     // Gravity const
    float Sigma; // surface tension

  public:
    virtual ~BaseSimulator()                                  = default;
    virtual int getWidth() const                              = 0;
    virtual int getHeight() const                             = 0;
    virtual const std::vector<Particle> &getParticles() const = 0;
    virtual std::vector<Line> getLines() const                = 0;
    virtual bool is_water(int x, int y) const                 = 0;
    virtual void update(float dt)                             = 0;

    // get functions for main loop
    float getGravity() { return G; }
    float getSurfaceTension() { return Sigma; }

    void setGravity(float gravity) { G = gravity; }
    void setSigma(float sigma) { Sigma = sigma; }

    virtual void
    addWater(float x, float y, float radius) = 0; // add water and dye
    virtual void
    delWater(float x, float y, float radius) = 0; // delete water and dye

    virtual void reset() = 0;
};
