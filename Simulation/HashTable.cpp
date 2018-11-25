
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include <iostream>
#include <fstream>
#include <algorithm>

#include <glm/gtx/norm.hpp>

#include "HashTable.hpp"

#include "Statistics/Statistics.hpp"
#include "Util/ZValue.hpp"
#include "Kernels.hpp"

//********************************************************************************
void HashTable::setHashTable(const glm::dvec2& minPos, const glm::dvec2& dx)
//********************************************************************************
{
    minPos_ = minPos-glm::dvec2(Kernel::SmoothingLength::h);
    gridSize_.x = unsigned(floor(dx.x*Kernel::SmoothingLength::dh))+2;
    gridSize_.y = unsigned(floor(dx.y*Kernel::SmoothingLength::dh))+2;

    particlesIn_.resize(gridSize_.x, gridSize_.y);
}

//********************************************************************************
void HashTable::writeGridRAW(std::string fileName) const
//********************************************************************************
{
    std::string fileNameGNU = fileName+".raw";
    std::cout<<"#writing file "<<fileNameGNU<<std::endl;
    std::ofstream outfile (fileNameGNU.c_str(),std::ofstream::binary);
    for (unsigned i=0;i<gridSize_.x;i++)
    {
        for (unsigned j=0;j<gridSize_.y;j++)
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

//********************************************************************************
void HashTable::write() const
//********************************************************************************
{
    for (unsigned j=0;j<gridSize_.y;j++)
    {
        for (unsigned i=0;i<gridSize_.x;i++)
        {
            std::cout<<" "<<particlesIn_(i,j).size();
        }
        std::cout<<std::endl;
    }
}

//********************************************************************************
void HashTable::clear()
//********************************************************************************
{
    #pragma omp parallel for
    for (unsigned i=0;i<gridSize_.x;i++)
    {
        #pragma omp parallel for
        for (unsigned j=0;j<gridSize_.y;j++)
        {
            if(particlesIn_(i,j).size())
                particlesIn_(i,j).clear();
        }
    }
}

//********************************************************************************
glm::ivec2 HashTable::findGridPos(glm::dvec2 pos) const
//********************************************************************************
{
    static auto findGridPosTimerID= Statistics::createTimer("HashTable::findGridPos");
    Statistics::TimerGuard findGridPosTimerGuard(findGridPosTimerID);

    glm::dvec2 dgrdPos = (pos-minPos_)*Kernel::SmoothingLength::dh;
    glm::ivec2 gridPos = glm::ivec2(int(floor(dgrdPos.x)),int(floor(dgrdPos.y)));
    while (gridPos.x<      0          ){ gridPos.x+=gridSize_.x; }
    while (gridPos.x>=int(gridSize_.x)){ gridPos.x-=gridSize_.x; }
    while (gridPos.y<      0          ){ gridPos.y+=gridSize_.y; }
    while (gridPos.y>=int(gridSize_.y)){ gridPos.y-=gridSize_.y; }

    return gridPos;
}

//********************************************************************************
void HashTable::findNei(std::vector<Particle>& cloud, const unsigned NParticles)
//********************************************************************************
{
    static auto findNeiTimerID = Statistics::createTimer("HashTable::findNei");
    Statistics::TimerGuard findNeiTimerGuard(findNeiTimerID);

    //find grid position of each particle
    for (unsigned i=0;i<NParticles;i++)
    {
        glm::dvec2 pos = cloud[i].position;
        glm::dvec2 dgrdPos = (pos-minPos_)*Kernel::SmoothingLength::dh;
        glm::ivec2 gridPos = glm::ivec2(int(floor(dgrdPos.x)),int(floor(dgrdPos.y)));
        while (gridPos.x<      0          ){ gridPos.x+=gridSize_.x; }
        while (gridPos.x>=int(gridSize_.x)){ gridPos.x-=gridSize_.x; }
        while (gridPos.y<      0          ){ gridPos.y+=gridSize_.y; }
        while (gridPos.y>=int(gridSize_.y)){ gridPos.y-=gridSize_.y; }

        particlesIn_(gridPos.x,gridPos.y).push_back(i);
        cloud[i].gridPos = gridPos;
    }

    //find nei of each particle
    #pragma omp parallel for
    for (unsigned i=0;i<NParticles;i++)
    {
        auto& iParticle = cloud[i];
        for (int iGrid=-1;iGrid<2;iGrid++)
        {
            for (int jGrid=-1;jGrid<2;jGrid++)
            {
                glm::ivec2 gridPos = iParticle.gridPos + glm::ivec2(iGrid,jGrid);

                if      (gridPos.x<      0          ){ gridPos.x+=gridSize_.x; }
                else if (gridPos.x>=int(gridSize_.x)){ gridPos.x-=gridSize_.x; }
                if      (gridPos.y<      0          ){ gridPos.y+=gridSize_.y; }
                else if (gridPos.y>=int(gridSize_.y)){ gridPos.y-=gridSize_.y; }

                const size_t nNei = particlesIn_(gridPos.x,gridPos.y).size();
                for (size_t iNei=0;iNei<nNei;iNei++)
                {
                    const unsigned neiPos = particlesIn_(gridPos.x,gridPos.y)[iNei];
                    const double dist2 = glm::length2(iParticle.position - cloud[neiPos].position);
                    if (dist2 < Kernel::SmoothingLength::h2)
                    {
                        iParticle.nei.push_back(Neigbhor(neiPos));
                    }
                }
            }
        }
    }

    auto comp = [&cloud](const Neigbhor& i, const Neigbhor& j){
      return cloud[i.ID].position.x < cloud[j.ID].position.x;
    };
    #pragma omp parallel for
    for (unsigned i=0;i<NParticles;i++)
    {
      auto& nei = cloud[i].nei;
      std::sort(nei.begin(), nei.end(), comp);
    }
}

//********************************************************************************
std::vector<unsigned>& HashTable::neiParticlesFor(glm::ivec2 gridPos)
//********************************************************************************
{
    if      (gridPos.x<      0          ){ gridPos.x+=gridSize_.x; }
    else if (gridPos.x>=int(gridSize_.x)){ gridPos.x-=gridSize_.x; }
    if      (gridPos.y<      0          ){ gridPos.y+=gridSize_.y; }
    else if (gridPos.y>=int(gridSize_.y)){ gridPos.y-=gridSize_.y; }
    return particlesIn_(gridPos.x,gridPos.y);
}

//********************************************************************************
const std::vector<unsigned>& HashTable::neiParticlesFor(glm::ivec2 gridPos) const
//********************************************************************************
{
    if      (gridPos.x<      0          ){ gridPos.x+=gridSize_.x; }
    else if (gridPos.x>=int(gridSize_.x)){ gridPos.x-=gridSize_.x; }
    if      (gridPos.y<      0          ){ gridPos.y+=gridSize_.y; }
    else if (gridPos.y>=int(gridSize_.y)){ gridPos.y-=gridSize_.y; }
    return particlesIn_(gridPos.x,gridPos.y);
}
