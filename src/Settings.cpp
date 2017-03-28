
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
int SPHSettings::LParticles = 20;
int SPHSettings::NParticles = 3000;
      
double SPHSettings::initDx = 0.02;

double     SPHSettings::particleDensity = 998.29;
double     SPHSettings::dParticleDensity = 1.0/particleDensity;
double     SPHSettings::particleMass = particleDensity*initDx*initDx*0.98;
double     SPHSettings::stiffness  = 1000.;// for WCSPH
double     SPHSettings::densityErr = 0.01; // for PCISPH - %
double     SPHSettings::viscosity = 10.;   
double     SPHSettings::surfTension = 0.2;
glm::dvec2 SPHSettings::grav = glm::dvec2(0.0,-1.e1);
//glm::dvec2 SPHSettings::grav = glm::dvec2(0.0,0.e1);

//----------------SimulationSettings-----------------//
double SimulationSettings::dt         = 1e-3;
double SimulationSettings::simTimeEnd = 15.;
double SimulationSettings::simTime    = 0.0;
int    SimulationSettings::showProgressEvery = 20000;

//----------------InitialConditions-----------------//
InitialConditions::ParticleGeneration 
  //InitialConditions::particleGeneration = InitialConditions::ALLIN;
  //InitialConditions::particleGeneration = InitialConditions::FAUCET;
    InitialConditions::particleGeneration = InitialConditions::DRIPPING;
double     InitialConditions::particleGenTime = 0.8; //for DRIPPING only
glm::dvec2 InitialConditions::particleInitVel = glm::dvec2(0.00,0.00);
glm::dvec2 InitialConditions::particleInitPos = glm::dvec2(0.30,1.55);

//----------------BoundaryConditions-----------------//
BoundingBox<glm::dvec2> BoundaryConditions::bndBox(glm::dvec2(0.0,0.0),glm::dvec2(1.0,5.0));
double                  BoundaryConditions::bndCoeff = 0.5;
double                  BoundaryConditions::bndConditionRange = 2.0;


//----------------RenderSettings-----------------//
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
        float(2.*BoundaryConditions::bndBox.maxX()+0.1*BoundaryConditions::bndBox.dx())
    )
);

int RenderSettings::width  = 512;
int RenderSettings::height = 1024;
int RenderSettings::printEvr= 16000; //print to file
RenderSettings::FileRenderType    RenderSettings::fileRender= RenderSettings::DISCRETE;
//RenderSettings::DisplayRenderType RenderSettings::displayRender= RenderSettings::INDEX;
RenderSettings::DisplayRenderType RenderSettings::displayRender= RenderSettings::SIMPLE;
