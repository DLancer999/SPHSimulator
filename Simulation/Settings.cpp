
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include "Settings.hpp"

//----------------SPHSettings-----------------//
//SPHSettings::Solver SPHSettings::SPHstep = SPHSettings::Solver::WCSPH;
SPHSettings::Solver SPHSettings::SPHstep = SPHSettings::Solver::PCISPH;

//field Particles
int SPHSettings::LParticles = 15;
int SPHSettings::NParticles = 2000;
      
double SPHSettings::initDx = 0.03;

double     SPHSettings::particleDensity = 998.29;
double     SPHSettings::dParticleDensity = 1.0/particleDensity;
double     SPHSettings::particleMass = particleDensity*initDx*initDx*0.98;
double     SPHSettings::stiffness  = 1000.;// for WCSPH
double     SPHSettings::densityErr = 0.01; // for PCISPH - %
double     SPHSettings::viscosity = 10.;   
double     SPHSettings::surfTension = 0.3;
glm::dvec2 SPHSettings::grav = glm::dvec2(0.0,-9.81e0);

//----------------SimulationSettings-----------------//
double SimulationSettings::dt         = 1e-3;
double SimulationSettings::simTimeEnd = 10.;
double SimulationSettings::simTime    = 0.0;
int    SimulationSettings::showProgressEvery = 20000;

//----------------InitialConditions-----------------//
InitialConditions::ParticleGeneration 
  //InitialConditions::particleGeneration = InitialConditions::ALLIN;
    InitialConditions::particleGeneration = InitialConditions::FAUCET;
  //InitialConditions::particleGeneration = InitialConditions::DRIPPING;
double     InitialConditions::particleGenTime = 0.2; //for DRIPPING only
glm::dvec2 InitialConditions::particleInitVel = glm::dvec2(1.00,-2.00);
glm::dvec2 InitialConditions::particleInitPos = glm::dvec2(0.15,1.65);

//----------------BoundaryConditions-----------------//
BoundingBox<glm::dvec2> BoundaryConditions::bndBox(glm::dvec2(0.0,0.0),glm::dvec2(2.0,5.0));
double                  BoundaryConditions::bndCoeff = 2.;
double                  BoundaryConditions::bndConditionRange = 4.0;


//----------------RenderSettings-----------------//
int RenderSettings::width  = 640;
int RenderSettings::height = 640;
int RenderSettings::printEvr= 1500; //print to file
RenderSettings::FileRenderType    RenderSettings::fileRender= RenderSettings::METABALL;
RenderSettings::DisplayRenderType RenderSettings::displayRender= RenderSettings::SIMPLE;

BoundingBox<glm::vec2> RenderSettings::displayBox
(
    glm::vec2
    ( 
        float(BoundaryConditions::bndBox.minX()-0.1*BoundaryConditions::bndBox.dx()), 
        float(BoundaryConditions::bndBox.minY()-0.1*BoundaryConditions::bndBox.dx())
    ),
    glm::vec2
    (
        float(BoundaryConditions::bndBox.maxX()+0.1*BoundaryConditions::bndBox.dx()),  
        float(double(height)/double(width)*BoundaryConditions::bndBox.maxX()+0.1*BoundaryConditions::bndBox.dx())
    )
);

