
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
#include <glm/gtx/norm.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "Kernels.hpp"
#include "Particle.hpp"

#include "Statistics/Statistics.hpp"

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
    static int JoinBits(int a, int b) 
    {
        long int result = 0;
        for(int ii = 8; ii >= 0; ii--)
        {
            result |= (a >> ii) & 1;
            result <<= 1;
            result |= (b >> ii) & 1;
            if(ii != 0){
                result <<= 1;
            }
        }
        return (int)result;
    }
    

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

    void setHashTable(const glm::dvec2& minPos, const glm::dvec2& dx)
    {
        minPos_ = minPos-glm::dvec2(Kernel::SmoothingLength::h);
        gridSize_.x = int(floor(dx.x*Kernel::SmoothingLength::dh))+2;
        gridSize_.y = int(floor(dx.y*Kernel::SmoothingLength::dh))+2;

        particlesIn_.resize(gridSize_.x, gridSize_.y);

        int ZMapSize = std::max(gridSize_.x,gridSize_.y);
        int pow2 = 0;
        while (ZMapSize>0)
        {
            ZMapSize>>=1;
            pow2++;
        }
        ZMapSize=1;
        for (int i=0;i<pow2;i++) ZMapSize<<=1;
        ZMapSize*=ZMapSize;

        gridZMap_.resize(ZMapSize);
        for (int i=0;i<ZMapSize;i++) gridZMap_[i] = glm::ivec2(-1,-1);

        gridZindex_.resize(gridSize_.x, gridSize_.y);
        for (int i=0;i<gridSize_.x;i++)
        {
            for (int j=0;j<gridSize_.y;j++)
            {
                int Zindex = JoinBits(i,j);
                gridZindex_(i,j)=Zindex;
                gridZMap_[Zindex] = glm::ivec2(i,j);
            }
        }
    }

    void writeGridRAW(std::string fileName) const
    {
        std::string fileNameGNU = fileName+".raw";
        std::cout<<"#writing file "<<fileNameGNU<<std::endl;
        std::ofstream outfile (fileNameGNU.c_str(),std::ofstream::binary);
        for (int i=0;i<gridSize_.x;i++)
        {
            for (int j=0;j<gridSize_.y;j++)
            {
                glm::dvec2 point1 = minPos_+glm::dvec2((i  )*Kernel::SmoothingLength::h,(j  )*Kernel::SmoothingLength::h);
                glm::dvec2 point2 = minPos_+glm::dvec2((i+1)*Kernel::SmoothingLength::h,(j  )*Kernel::SmoothingLength::h);
                glm::dvec2 point3 = minPos_+glm::dvec2((i+1)*Kernel::SmoothingLength::h,(j+1)*Kernel::SmoothingLength::h);
                glm::dvec2 point4 = minPos_+glm::dvec2((i  )*Kernel::SmoothingLength::h,(j+1)*Kernel::SmoothingLength::h);
                glm::dvec2 point5 = minPos_+glm::dvec2((i  )*Kernel::SmoothingLength::h,(j  )*Kernel::SmoothingLength::h);
                outfile<<point1.x<<"\t"<<point1.y<<"\n" 
                       <<point2.x<<"\t"<<point2.y<<"\n" 
                       <<point3.x<<"\t"<<point3.y<<"\n" 
                       <<point4.x<<"\t"<<point4.y<<"\n" 
                       <<point5.x<<"\t"<<point5.y<<"\n\n";
            }
        }
        outfile.close();
    }

    void write() const
    {
        for (int j=gridSize_.y-1;j>=0;j--)
        {
            for (int i=0;i<gridSize_.x;i++)
            {
                std::cout<<" "<<particlesIn_(i,j).size();
            }
            std::cout<<std::endl;
        }
    }

    void clear()
    {
        #pragma omp parallel for
        for (int i=0;i<gridSize_.x;i++)
        {
            #pragma omp parallel for
            for (int j=0;j<gridSize_.y;j++)
            {
                if(particlesIn_(i,j).size())
                    particlesIn_(i,j).clear();
            }
        }
    }

    glm::ivec2 findGridPos(glm::dvec2 pos) const
    {
        static auto findGridPosTimerID= Statistics::createTimer("HashTable::findGridPos");
        Statistics::TimerGuard findGridPosTimerGuard(findGridPosTimerID);

        glm::dvec2 dgrdPos = (pos-minPos_)*Kernel::SmoothingLength::dh;
        glm::ivec2 gridPos = glm::ivec2(int(floor(dgrdPos.x)),int(floor(dgrdPos.y)));
        while (gridPos.x<      0     ){ gridPos.x+=gridSize_.x; }
        while (gridPos.x>=gridSize_.x){ gridPos.x-=gridSize_.x; }
        while (gridPos.y<      0     ){ gridPos.y+=gridSize_.y; }
        while (gridPos.y>=gridSize_.y){ gridPos.y-=gridSize_.y; }

        return gridPos;
    }

    void findNei(std::vector<Particle>& cloud, const int NParticles)
    {
        static auto findNeiTimerID = Statistics::createTimer("HashTable::findNei");
        Statistics::TimerGuard findNeiTimerGuard(findNeiTimerID);

        //find grid position of each particle
        for (int i=0;i<NParticles;i++)
        {
            glm::dvec2 pos = cloud[i].position;
            glm::dvec2 dgrdPos = (pos-minPos_)*Kernel::SmoothingLength::dh;
            glm::ivec2 gridPos = glm::ivec2(int(floor(dgrdPos.x)),int(floor(dgrdPos.y)));
            while (gridPos.x<      0     ){ gridPos.x+=gridSize_.x; }
            while (gridPos.x>=gridSize_.x){ gridPos.x-=gridSize_.x; }
            while (gridPos.y<      0     ){ gridPos.y+=gridSize_.y; }
            while (gridPos.y>=gridSize_.y){ gridPos.y-=gridSize_.y; }

            particlesIn_(gridPos.x,gridPos.y).push_back(i);
            cloud[i].gridPos = gridPos;
        }

        //find nei of each particle
        #pragma omp parallel for
        for (int i=0;i<NParticles;i++)
        {
            auto& iParticle = cloud[i];
            for (int iGrid=-1;iGrid<2;iGrid++)
            {
                for (int jGrid=-1;jGrid<2;jGrid++)
                {
                    glm::ivec2 gridPos = iParticle.gridPos + glm::ivec2(iGrid,jGrid);

                    if      (gridPos.x<      0     ){ gridPos.x+=gridSize_.x; }
                    else if (gridPos.x>=gridSize_.x){ gridPos.x-=gridSize_.x; }
                    if      (gridPos.y<      0     ){ gridPos.y+=gridSize_.y; }
                    else if (gridPos.y>=gridSize_.y){ gridPos.y-=gridSize_.y; }

                    const size_t nNei = particlesIn_(gridPos.x,gridPos.y).size();
                    for (size_t iNei=0;iNei<nNei;iNei++)
                    {
                        const int neiPos = particlesIn_(gridPos.x,gridPos.y)[iNei];
                        const double dist2 = glm::length2(iParticle.position - cloud[neiPos].position);
                        if (dist2 < Kernel::SmoothingLength::h2)
                        {
                            iParticle.nei.push_back(Neigbhor(neiPos));
                        }
                    }
                }
            }
        }
    }

    void reorderCloud(std::vector<Particle>& cloud, const int NParticles)
    {
        static auto reorTimerID = Statistics::createTimer("HashTable::reorderTimer");
        Statistics::TimerGuard reorderTimerGuard(reorTimerID);

        std::vector<Particle> newCloud(NParticles);
        std::vector<int>      oldToNewMap;

        const int NZcells = (int)gridZMap_.size();
        for (int iZ=0;iZ<NZcells;iZ++)
        {
            glm::ivec2 gridPos = gridZMap_[iZ];
            if (gridPos.x<0) continue;
            const int NPartsInCell = (int)particlesIn_(gridPos.x,gridPos.y).size();
            for (int iPart=0;iPart<NPartsInCell;iPart++)
            {
                oldToNewMap.push_back(particlesIn_(gridPos.x,gridPos.y)[iPart]);
            }
        }

        for (int iPart=0;iPart<NParticles;iPart++)
        {
            newCloud[iPart] = cloud[oldToNewMap[iPart]];
        }

        for (int iPart=0;iPart<NParticles;iPart++)
        {
            cloud[iPart] = newCloud[iPart];
        }
    }

    //return particleList of given grid cell
    std::vector<int>& neiParticlesFor(glm::ivec2 gridPos)
    {
        if      (gridPos.x<      0     ){ gridPos.x+=gridSize_.x; }
        else if (gridPos.x>=gridSize_.x){ gridPos.x-=gridSize_.x; }
        if      (gridPos.y<      0     ){ gridPos.y+=gridSize_.y; }
        else if (gridPos.y>=gridSize_.y){ gridPos.y-=gridSize_.y; }
        return particlesIn_(gridPos.x,gridPos.y);
    }

    //return particleList of given grid cell
    const std::vector<int>& neiParticlesFor(glm::ivec2 gridPos) const
    {
        if      (gridPos.x<      0     ){ gridPos.x+=gridSize_.x; }
        else if (gridPos.x>=gridSize_.x){ gridPos.x-=gridSize_.x; }
        if      (gridPos.y<      0     ){ gridPos.y+=gridSize_.y; }
        else if (gridPos.y>=gridSize_.y){ gridPos.y-=gridSize_.y; }
        return particlesIn_(gridPos.x,gridPos.y);
    }
};

#endif
