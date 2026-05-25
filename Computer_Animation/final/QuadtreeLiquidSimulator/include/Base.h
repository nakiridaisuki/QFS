#pragma once

#include <vector>

struct Particle {
    float x, y;
    float u, v;
};

struct Line {
    int x1, y1, x2, y2;
};

class BaseSimulator {
  public:
    virtual ~BaseSimulator() = default;
    virtual int getWidth() const = 0;
    virtual int getHeight() const = 0;
    virtual const std::vector<int> &getCell() const = 0;
    virtual const std::vector<Particle> &getParticles() const = 0;
    virtual std::vector<Line> getLines() const = 0;
};
