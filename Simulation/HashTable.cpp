
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

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm_ext.hpp>

using boost::adaptors::indexed;

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
    for (unsigned i=0;i<gridSize_.x;i++)
    {
        for (unsigned j=0;j<gridSize_.y;j++)
        {
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
    if ( dgrdPos.x < 0. ) dgrdPos.x = 0.;
    if ( dgrdPos.y < 0. ) dgrdPos.y = 0.;
    glm::ivec2 gridPos = glm::ivec2(int(floor(dgrdPos.x)),int(floor(dgrdPos.y)));
    if ( gridPos.x >= int(gridSize_.x) ) gridPos.x = int(gridSize_.x);
    if ( gridPos.y >= int(gridSize_.y) ) gridPos.y = int(gridSize_.y);

    return gridPos;
}

//********************************************************************************
void HashTable::findNei(ParticleCloud& cloud)
//********************************************************************************
{
    static auto findNeiTimerID = Statistics::createTimer("HashTable::findNei");
    Statistics::TimerGuard findNeiTimerGuard(findNeiTimerID);
    //static auto findNeiTimerID_1 = Statistics::createTimer("HashTable::findNei_P1");
    //static auto findNeiTimerID_2 = Statistics::createTimer("HashTable::findNei_P2");
    //static auto findNeiTimerID_3 = Statistics::createTimer("HashTable::findNei_P3");

    const auto& particlePos = cloud.get<Attr::ePosition>();

    auto& particleNei = cloud.get<Attr::eNei>();

    std::vector<glm::ivec2> particleGridPos(particlePos.size(), glm::ivec2(0));

    const size_t nPart = cloud.size();

    {
        //Statistics::TimerGuard findNeiTimerGuard_1(findNeiTimerID_1);
        const double dh = Kernel::SmoothingLength::dh;
        //find grid position of each particle
        for (size_t iPart=0; iPart<nPart; ++iPart)
        {
            const glm::dvec2 pos = particlePos[iPart];
            const glm::dvec2 dgrdPos = (pos-minPos_)*dh;
            glm::ivec2 gridPos = glm::ivec2(int(floor(dgrdPos.x)),int(floor(dgrdPos.y)));

            particlesIn_(gridPos.x,gridPos.y).push_back(unsigned(iPart));
            particleGridPos[iPart] = gridPos;
        }
    }

    {
        //Statistics::TimerGuard findNeiTimerGuard_2(findNeiTimerID_2);
        const double h2 = Kernel::SmoothingLength::h2;
        //find nei of each particle
        #pragma omp parallel for
        for (size_t iPart=0; iPart<nPart; ++iPart)
        {
            const glm::dvec2 iPos = particlePos[iPart];
            const glm::ivec2 iGridPos = particleGridPos[iPart];
            auto& iNei = particleNei[iPart];

            auto iRange = boost::irange(
                (iGridPos.x > 0) ? -1 : 0,
                (iGridPos.x == int(gridSize_.x-1)) ? 1 : 2
            );

            auto jRange = boost::irange(
                (iGridPos.y > 0) ? -1 : 0,
                (iGridPos.y == int(gridSize_.y-1)) ? 1 : 2
            );

            for (int iGrid : iRange)
            {
                for (int jGrid : jRange)
                {
                    glm::ivec2 gridPos = iGridPos + glm::ivec2(iGrid,jGrid);

                    for ( unsigned neiPos : particlesIn_(gridPos.x,gridPos.y) ) {
                        const double dist2 = glm::length2(iPos - particlePos[neiPos]);
                        if (dist2 < h2)
                        {
                            iNei.push_back(Neigbhor(neiPos));
                        }
                    }
                }
            }
        }
    }

    {
        //Statistics::TimerGuard findNeiTimerGuard_3(findNeiTimerID_3);
        auto comp = [&particlePos](const Neigbhor& i, const Neigbhor& j){
          return particlePos[i.ID].x < particlePos[j.ID].x;
        };

        #pragma omp parallel for
        for (size_t iPart=0; iPart<nPart; ++iPart)
        {
          auto& nei = particleNei[iPart];
          std::sort(nei.begin(), nei.end(), comp);
        }
    }
}

//********************************************************************************
std::vector<unsigned> HashTable::neiParticlesFor(glm::ivec2 gridPos) const
//********************************************************************************
{
    std::vector<unsigned> ret;

    auto iRange = boost::irange(
        (gridPos.x > 0) ? -1 : 0,
        (gridPos.x == int(gridSize_.x-1)) ? 1 : 2
    );
    auto jRange = boost::irange(
        (gridPos.y > 0) ? -1 : 0,
        (gridPos.y == int(gridSize_.y-1)) ? 1 : 2
    );

    for (int iGrid : iRange)
    {
        for (int jGrid : jRange)
        {
            glm::ivec2 sideGridPos = gridPos + glm::ivec2(iGrid,jGrid);
            boost::push_back(ret, particlesIn_(sideGridPos.x,sideGridPos.y));
        }
    }
    return ret;
}
