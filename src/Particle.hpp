#ifndef PARTICLE_H
#define PARTICLE_H

#include <glm/glm.hpp>
#include <vector>

class Particle
{
public:
    glm::dvec2 position;
    glm::dvec2 velocity;
    glm::dvec2 Fpress;
    glm::dvec2 Fvisc;
    glm::dvec2 Fother;
    glm::dvec2 Ftot;
    glm::ivec2 gridPos;
    double mass;
    double density;
    double ddensity;
    double densityErr;
    double pressure;
    bool   active;
    std::vector<int> nei;

public:
    Particle():
    position(-100.0),
    velocity(0.0),
    Fpress(0.0),
    Fvisc(0.0),
    Fother(0.0),
    Ftot(0.0),
    gridPos(0),
    mass(0.0),
    density(0.0),
    ddensity(0.0),
    densityErr(0.0),
    pressure(0.0),
    active(true),
    nei()
    {}
};

class ScreenPoint
{
public:
    glm::vec2 position;
    glm::vec3 color;
    glm::vec2 force;

public:
    ScreenPoint():
    position(0.0f),
    color(0.0f),
    force(0.0f)
    {}
};

#endif
