
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

SourceFiles
    -

\************************************************************************/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <omp.h>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <fstream>

#include "Kernels.hpp"

typedef std::vector< std::vector< std::vector<int> > > grid;

class HashTable
{
private:
    grid  particlesIn_;
    glm::dvec2 minPos_;
    glm::ivec2 gridSize_;

public:
    HashTable():
    particlesIn_(),
    minPos_(0.0),
    gridSize_(0,0)
    {
    }

    void setHashTable(const glm::dvec2& minPos, const glm::dvec2& dx)
    {
        minPos_ = minPos-glm::dvec2(Kernel::SmoothingLength::h);
        gridSize_.x = int(floor(dx.x*Kernel::SmoothingLength::dh))+2;
        gridSize_.y = int(floor(dx.y*Kernel::SmoothingLength::dh))+2;

        particlesIn_.resize(gridSize_.x);
        for (int i=0;i<gridSize_.x;i++)
        {
            particlesIn_[i].resize(gridSize_.y);
        }
    }

    void writeGridRAW(std::string fileName)
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

    void write()
    {
        for (int j=gridSize_.y-1;j>=0;j--)
        {
            for (int i=0;i<gridSize_.x;i++)
            {
                std::cout<<" "<<particlesIn_[i][j].size();
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
                if(particlesIn_[i][j].size())
                    particlesIn_[i][j].clear();
            }
        }
    }

    glm::ivec2 findGridPos(glm::dvec2 pos)
    {
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
        //clear lists
        clear();
        #pragma omp parallel for
        for (int i=0;i<NParticles;i++)
        {
            if(cloud[i].nei.size())
                cloud[i].nei.clear();
        }

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

            particlesIn_[gridPos.x][gridPos.y].push_back(i);
            cloud[i].gridPos = gridPos;
        }

        //find nei of each particle
        #pragma omp parallel for
        for (int i=0;i<NParticles;i++)
        {
            for (int iGrid=-1;iGrid<2;iGrid++)
            {
                for (int jGrid=-1;jGrid<2;jGrid++)
                {
                    glm::ivec2 gridPos = cloud[i].gridPos + glm::ivec2(iGrid,jGrid);

                    if      (gridPos.x<      0     ){ gridPos.x+=gridSize_.x; }
                    else if (gridPos.x>=gridSize_.x){ gridPos.x-=gridSize_.x; }
                    if      (gridPos.y<      0     ){ gridPos.y+=gridSize_.y; }
                    else if (gridPos.y>=gridSize_.y){ gridPos.y-=gridSize_.y; }

                    int nNei = (int)particlesIn_[gridPos.x][gridPos.y].size();
                    for (int iNei=0;iNei<nNei;iNei++)
                    {
                        int neiPos = particlesIn_[gridPos.x][gridPos.y][iNei];
                        double dist2 = glm::length2(cloud[i].position - cloud[neiPos].position);
                        if (dist2 < Kernel::SmoothingLength::h2)
                        {
                            cloud[i].nei.push_back(neiPos);
                        }
                    }
                }
            }
        }
    }

    //return particleList of given grid cell
    std::vector<int>& neiParticlesFor(glm::ivec2 gridPos)
    {
        if      (gridPos.x<      0     ){ gridPos.x+=gridSize_.x; }
        else if (gridPos.x>=gridSize_.x){ gridPos.x-=gridSize_.x; }
        if      (gridPos.y<      0     ){ gridPos.y+=gridSize_.y; }
        else if (gridPos.y>=gridSize_.y){ gridPos.y-=gridSize_.y; }
        return particlesIn_[gridPos.x][gridPos.y];
    }
};

#endif
