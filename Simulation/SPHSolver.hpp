
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    SPHSolver
 
Description
    Class to update position/velocity of each particle

SourceFiles
    SPHSolver.cpp

\************************************************************************/

#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include <vector>

#include "Particle.hpp"
#include "HashTable.hpp"
#include "ParticleCloud.hpp"

class SPHSolver
{
protected:
    ParticleCloud cloud_;  //list of particles
    HashTable     neibhs_; 

    double delta_;           //used for PCISPH

    void generateParticles(); //activate and initialize new particles when needed
    void updateNei();         //update neibghor values
    void calcDensity();       //calculate density
    void calcDensityErr();    //calculate density error
    void initPressure();      //set pressure to 0
    void calcPressure();      //calc pressure - EOS
    void calcNormal();        //calc normal
    void calcPressForces();   //calc pressure forces
    void calcViscForces();    //calc viscous forces
    void calcSurfForces();    //calc viscous forces
    void calcOtherForces();   //calc gravity/boundary conditions (repulsive forces)
    void updatePressure();    //update pressure in PCISPH

    void WCSPHStep();  //WCSPH  solver
    void PCISPHStep(); //PCISPH solver
public:
    SPHSolver():
    cloud_(),
    neibhs_(),
    delta_(0.0)
    {}

    void init();
    bool step();

    double calcCFL() const;
    double calcKineticEnergy() const;

          ParticleCloud& cloud()       { return cloud_; }
    const ParticleCloud& cloud() const { return cloud_; }
          HashTable& neibhs()       {return neibhs_;}
    const HashTable& neibhs() const {return neibhs_;}
};
#endif
