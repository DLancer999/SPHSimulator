
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Various Settings classes
 
Description
    All parameters required to set-up simulation

SourceFiles
    Settings.cpp

\************************************************************************/

#ifndef SPHSETTINGS_H
#define SPHSETTINGS_H

#include "BoundingBox.hpp"

class SPHSettings
{
public:
    enum Solver
    {
        WCSPH,
        PCISPH
    };
    
    static Solver SPHstep;
    
    //field Particles
    static int LParticles; //used to control particle generation
    static int NParticles; //total particles
          
    static double initDx; //initial particle distance

    static double particleDensity;
    static double dParticleDensity;
    static double particleMass;
    static double stiffness;
    static double densityErr;
    static double viscosity;
    static double surfTension;
    static glm::dvec2 grav;
};

class SimulationSettings
{
public:
    static double dt;                //simulation timestep
    static double simTimeEnd;        //time on which simulation ends
    static double simTime;           //running time
    static int    showProgressEvery; //print progress on screen

    static void updateSimTime()
    {
        simTime+=dt;
    }

    static bool breakLoop()
    {
        if (simTime>simTimeEnd) return true;
        else                    return false;
    }
};

class InitialConditions
{
public:
    enum ParticleGeneration
    {
        ALLIN,    //column of fluid
        FAUCET,   //tap like flow
        DRIPPING  //dropplets generation
    };
    static ParticleGeneration particleGeneration; //type of particle generation
    static glm::dvec2 particleInitPos;            //initial particle position
    static glm::dvec2 particleInitVel;            //initial particle velocity
    static double     particleGenTime;            //particles generated every ... seconds, used only for DRIPPING scheme
};

class BoundaryConditions
{
public:
    //background grid
    static BoundingBox<glm::dvec2> bndBox; //boundary wall box
    static double bndCoeff;                //multiplier to boundary force
    static double bndConditionRange;       //range of boundary forces (times h)
};

//rendering settings
class RenderSettings
{
public:
    enum FileRenderType
    {
        RAWDATA,  //write particle attributes to file
        DISCRETE, //render to png as particles
        METABALL  //render to png as metaball
    };
    enum DisplayRenderType
    {
        SIMPLE,      //render just the particles
        PRESSFORCES, //render the particles + pressureForce
        VISCFORCES,  //render the particles + viscForce
        SURFFORCES,  //render the particles + surfForce
        OTHERFORCES, //render the particles + otherForces
        ALLFORCES    //render the particles + totalForce
    };
    static BoundingBox<glm::vec2> displayBox; //rendering window
    static int width;                         //width resolution
    static int height;                        //height resolution
    static int printEvr;                      //write file every ... timesteps
    static FileRenderType fileRender;
    static DisplayRenderType displayRender;
};

#endif
