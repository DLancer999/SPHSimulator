
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

class Neigbhor
{
public:
    glm::dvec2 dir;
    double     dist;
    unsigned   ID;

public:
    Neigbhor():
    dir(0.0),
    dist(0.0),
    ID(0)
    {}

    Neigbhor(const unsigned neiID):
    dir(0.0),
    dist(0.0),
    ID(neiID)
    {}
};

class Particle
{
public:
    glm::dvec2 position = glm::dvec2(-100.0);
    glm::dvec2 velocity = glm::dvec2(0.0);
    glm::dvec2 normal   = glm::dvec2(0.0);
    glm::dvec2 Fpress   = glm::dvec2(0.0);
    glm::dvec2 Fvisc    = glm::dvec2(0.0);
    glm::dvec2 Fsurf    = glm::dvec2(0.0);
    glm::dvec2 Fother   = glm::dvec2(0.0);
    glm::dvec2 Ftot     = glm::dvec2(0.0);
    double mass       = 0.0;
    double density    = 0.0;
    double ddensity   = 0.0;
    double densityErr = 0.0;
    double pressure   = 0.0;
    std::vector<Neigbhor> nei = std::vector<Neigbhor>();
};

#endif
