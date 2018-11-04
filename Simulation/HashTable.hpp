
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
    -- both AABB and ZSorting implementations are memory hungry... brute force approach

SourceFiles
    -

\************************************************************************/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <glm/glm.hpp>
#include <string>

#include "Particle.hpp"
#include "Util/Matrix.hpp"

typedef std::vector< glm::ivec2 > ZMap;

class HashTable
{
public:
    static int JoinBits(int a, int b);

private:
    Util::Matrix<std::vector<unsigned>> particlesIn_;
    Util::Matrix<int> gridZindex_;
    ZMap   gridZMap_;
    glm::dvec2 minPos_;
    glm::uvec2 gridSize_;

public:
    HashTable():
    particlesIn_(),
    gridZindex_(),
    gridZMap_(),
    minPos_(0.0),
    gridSize_(0,0)
    {
    }

    void setHashTable(const glm::dvec2& minPos, const glm::dvec2& dx);

    void writeGridRAW(std::string fileName) const;

    void write() const;

    void clear();

    glm::ivec2 findGridPos(glm::dvec2 pos) const;

    void findNei(std::vector<Particle>& cloud, const unsigned NParticles);

    void reorderCloud(std::vector<Particle>& cloud, const unsigned NParticles);

    //return particleList of given grid cell
    std::vector<unsigned>& neiParticlesFor(glm::ivec2 gridPos);
    const std::vector<unsigned>& neiParticlesFor(glm::ivec2 gridPos) const;
};

#endif
