
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Particle 
 
Description
    All parameters required to set-up simulation

SourceFiles
    -

\************************************************************************/

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
    glm::ivec2 gridPos; //in background grid
    double mass;
    double density;
    double ddensity;    //1./density
    double densityErr;
    double pressure;
    std::vector<int> nei; //list of neighbours

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
    nei()
    {}
};

#endif
