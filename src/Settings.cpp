#include "Settings.hpp"

//----------------SPHSettings-----------------//
SPHSettings::Solver SPHSettings::SPHstep = SPHSettings::Solver::PCISPH;

//field Particles
int SPHSettings::LParticles = 30;
int SPHSettings::NParticles = 1500;
      
double SPHSettings::initDx = 0.02;

double     SPHSettings::particleDensity = 998.29;
double     SPHSettings::dParticleDensity = 1.0/particleDensity;
double     SPHSettings::particleMass = particleDensity*initDx*initDx*0.98;
double     SPHSettings::stiffness  = 2500.;// for WCSPH
double     SPHSettings::densityErr = 0.01; // for PCISPH - %
double     SPHSettings::viscosity = 5;   //molecular
glm::dvec2 SPHSettings::grav = glm::dvec2(0.0,-1.e1);

//----------------SimulationSettings-----------------//
double SimulationSettings::dt = 1e-3;
double SimulationSettings::simTimeEnd = 1;
double SimulationSettings::simTime = 0.0;

//----------------InitialConditions-----------------//
InitialConditions::ParticleGeneration 
    InitialConditions::particleGeneration = InitialConditions::ALLIN;
double     InitialConditions::particleGenTime = 0.4; //for DRIPPING only
glm::dvec2 InitialConditions::particleInitVel = glm::dvec2(0.00,-0.00);
glm::dvec2 InitialConditions::particleInitPos = glm::dvec2(0.02,0.02);

//----------------BoundaryConditions-----------------//
BoundingBox<glm::dvec2> BoundaryConditions::bndBox(glm::dvec2(0.0,0.0),glm::dvec2(1.0,5.0));
double                  BoundaryConditions::bndCoeff = 2.;
double                  BoundaryConditions::bndConditionRange = 3.0;


//----------------RenderSettings-----------------//
BoundingBox<glm::vec2> RenderSettings::displayBox(
                        glm::vec2( 
                            float(BoundaryConditions::bndBox.minX()-0.1*BoundaryConditions::bndBox.dx()), 
                            float(BoundaryConditions::bndBox.minY()-0.1*BoundaryConditions::bndBox.dx())
                            ),
                        glm::vec2(
                            float(BoundaryConditions::bndBox.maxX()+0.1*BoundaryConditions::bndBox.dx()),  
                            float(2.*BoundaryConditions::bndBox.maxX()+0.1*BoundaryConditions::bndBox.dx()))
                       );

int RenderSettings::width  = 512;
int RenderSettings::height = 1024;
int RenderSettings::printEvr= 100;
RenderSettings::FileRenderType RenderSettings::fileRender= RenderSettings::DISCRETE;
RenderSettings::DisplayRenderType RenderSettings::displayRender= RenderSettings::SIMPLE;
