
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    HashTable
 
Description
    Class for particle grouping
    -- Not actual hash table yet... just AABB background grid
    -- AABB is memory hungry... brute force approach

SourceFiles
    -

\************************************************************************/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <glm/glm.hpp>
#include <string>

#include "Particle.hpp"
#include "ParticleCloud.hpp"
#include "Util/Matrix.hpp"

class HashTable
{

private:
    Util::Matrix<std::vector<unsigned>> particlesIn_;
    std::vector<glm::uvec2> particleGridPos_;
    glm::dvec2 minPos_;
    glm::uvec2 gridSize_;

public:
    HashTable():
    particlesIn_(),
    particleGridPos_(),
    minPos_(0.0),
    gridSize_(0,0)
    {
    }

    void setHashTable(const glm::dvec2& minPos, const glm::dvec2& dx);

    void writeGridRAW(std::string fileName) const;

    void write() const;

    void clear();

    glm::uvec2 findGridPos(glm::dvec2 pos) const;

    void findNei(ParticleCloud& cloud);

    //return particleList of particles around grid cell
    std::vector<unsigned> neiParticlesFor(glm::uvec2 gridPos) const;
};

#endif
