
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

typedef std::vector< glm::ivec2 > ZMap;

template <typename T>
class Matrix
{
public:
  using iterator       = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  iterator       begin()        { return _data.begin(); }
  const_iterator begin()  const { return _data.begin(); }
  const_iterator cbegin() const { return _data.cbegin(); }
  iterator       end()        { return _data.end(); }
  const_iterator end()  const { return _data.end(); }
  const_iterator cend() const { return _data.cend(); }

  Matrix() :_data() ,_sizeX(0) ,_sizeY(0) {};
  Matrix(const Matrix&) = default;
  Matrix(Matrix&&)      = default;
  ~Matrix() = default;
  Matrix& operator=(const Matrix&) = default;
  Matrix& operator=(Matrix&&)      = default;

  void resize(size_t sizeX, size_t sizeY)
  {
    _sizeX = sizeX;
    _sizeY = sizeY;
    _data.resize(sizeX*sizeY, T{});
  }
  void clear()
  {
    _data.clear();
  }

  T& operator()(size_t i, size_t j) {
    return _data[i*_sizeY + j];
  }

  const T& operator()(size_t i, size_t j) const {
    return _data[i*_sizeY + j];
  }

private:
  std::vector<T> _data;
  size_t _sizeX;
  size_t _sizeY;
};

class HashTable
{
public:
    static int JoinBits(int a, int b);

private:
    Matrix<std::vector<int>> particlesIn_;
    Matrix<int> gridZindex_;
    ZMap   gridZMap_;
    glm::dvec2 minPos_;
    glm::ivec2 gridSize_;

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

    void findNei(std::vector<Particle>& cloud, const int NParticles);

    void reorderCloud(std::vector<Particle>& cloud, const int NParticles);

    //return particleList of given grid cell
    std::vector<int>& neiParticlesFor(glm::ivec2 gridPos);
    const std::vector<int>& neiParticlesFor(glm::ivec2 gridPos) const;
};

#endif
