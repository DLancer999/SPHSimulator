
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include "Settings.hpp"

#include <fstream>
#include <boost/program_options.hpp>

//----------------SPHSettings-----------------//
SPHSettings::Solver SPHSettings::SPHstep = SPHSettings::Solver::WCSPH;

//field Particles
unsigned SPHSettings::LParticles = 0;
unsigned SPHSettings::NParticles = 0;
      
double SPHSettings::initDx = 0.;

double     SPHSettings::particleDensity = 0.;
double     SPHSettings::dParticleDensity = 0.;
double     SPHSettings::particleMass = 0.;
double     SPHSettings::stiffness  = 0.;
double     SPHSettings::densityErr = 0.;
double     SPHSettings::viscosity = 0.;
double     SPHSettings::surfTension = 0.;
glm::dvec2 SPHSettings::grav{0.,0.};

//----------------SimulationSettings-----------------//
double SimulationSettings::dt         = 0.;
double SimulationSettings::simTimeEnd = 0.;
double SimulationSettings::simTime    = 0.;
int    SimulationSettings::showProgressEvery = 0;

//----------------InitialConditions-----------------//
InitialConditions::ParticleGeneration 
    InitialConditions::particleGeneration = InitialConditions::ALLIN;
double     InitialConditions::particleGenTime = 0.; //for DRIPPING only
glm::dvec2 InitialConditions::particleInitVel{0.,0.0};
glm::dvec2 InitialConditions::particleInitPos{0.,0.0};

//----------------BoundaryConditions-----------------//
BoundingBox<glm::dvec2> BoundaryConditions::bndBox(glm::dvec2(0.0,0.0),glm::dvec2(0.0,0.0));
double                  BoundaryConditions::bndCoeff = 0.;
double                  BoundaryConditions::bndConditionRange = 0.;


//----------------RenderSettings-----------------//
unsigned RenderSettings::width  = 0;
unsigned RenderSettings::height = 0;
int RenderSettings::printEvr= 0; //print to file
RenderSettings::FileRenderType    RenderSettings::fileRender= RenderSettings::RAWDATA;
RenderSettings::DisplayRenderType RenderSettings::displayRender= RenderSettings::SIMPLE;

BoundingBox<glm::vec2> RenderSettings::displayBox;

using namespace boost::program_options;

void SPHSettings::readSettings()
{
    std::string SPHstepMethod;
    double gravityForce = 0.;

    std::ifstream Config_File("SPHConfig.ini");
    options_description SimSet("Settings");
    SimSet.add_options()
        ("SPHSettings.StepSolver",      value<std::string>(&SPHstepMethod))
        ("SPHSettings.LParticles",      value<unsigned>(&LParticles))
        ("SPHSettings.NParticles",      value<unsigned>(&NParticles))
        ("SPHSettings.initDx",          value<double>(&initDx))
        ("SPHSettings.particleDensity", value<double>(&particleDensity))
        ("SPHSettings.stiffness",       value<double>(&stiffness))
        ("SPHSettings.densityErr",      value<double>(&densityErr))
        ("SPHSettings.viscosity",       value<double>(&viscosity))
        ("SPHSettings.surfTension",     value<double>(&surfTension))
        ("SPHSettings.gravity",         value<double>(&gravityForce));
    
    variables_map vm;
    store(parse_config_file(Config_File, SimSet), vm);
    notify(vm);

    //process...
    if (SPHstepMethod == "WCSPH") {
      SPHstep = Solver::WCSPH;
    }
    else if (SPHstepMethod == "PCISPH") {
      SPHstep = Solver::PCISPH;
    }
    else {
      std::string errorMsg = "Unrecognized StepSolver method - AvailableMethods::\n";
      errorMsg += "  WCSPH\n";
      errorMsg += "  PCISPH\n";
      throw std::runtime_error(errorMsg);
    }

    dParticleDensity = 1.0/particleDensity;
    particleMass = particleDensity*initDx*initDx*0.98;
    grav = glm::dvec2(0.0,gravityForce);
}

void SimulationSettings::readSettings()
{
    std::ifstream Config_File("SimConfig.ini");
    options_description SimSet("Settings");
    SimSet.add_options()
        ("SimulationSettings.stepDt",      value<double>(&dt))
        ("SimulationSettings.StartTime",   value<double>(&simTime))
        ("SimulationSettings.EndTime",     value<double>(&simTimeEnd))
        ("SimulationSettings.ReportEvery", value<int>(&showProgressEvery));

    variables_map vm;
    store(parse_config_file(Config_File, SimSet), vm);
    notify(vm);
}

void InitialConditions::readSettings()
{
    std::string particleGenerator;
    double initPosX;
    double initPosY;
    double initVelX;
    double initVelY;

    std::ifstream Config_File("InitConfig.ini");
    options_description ICSet("Settings");
    ICSet.add_options()
        ("InitialConditions.particleGenerator", value<std::string>(&particleGenerator))
        ("InitialConditions.initPosX",    value<double>(&initPosX))
        ("InitialConditions.initPosY",    value<double>(&initPosY))
        ("InitialConditions.initVelX",    value<double>(&initVelX))
        ("InitialConditions.initVelY",    value<double>(&initVelY))
        ("InitialConditions.genInterval", value<double>(&particleGenTime));

    variables_map vm;
    store(parse_config_file(Config_File, ICSet), vm);
    notify(vm);

    //process...
    if (particleGenerator == "ALLIN") {
      particleGeneration = ParticleGeneration::ALLIN;
    }
    else if (particleGenerator == "FAUCET") {
      particleGeneration = ParticleGeneration::FAUCET;
    }
    else if (particleGenerator == "DRIPPING") {
      particleGeneration = ParticleGeneration::DRIPPING;
    }
    else {
      std::string errorMsg = "Unrecognized ParticleGenerator - Available Generators::\n";
      errorMsg += "  ALLIN\n";
      errorMsg += "  FAUCET\n";
      errorMsg += "  DRIPPING\n";
      throw std::runtime_error(errorMsg);
    }
    
    particleInitPos = glm::dvec2{initPosX, initPosY};
    particleInitVel = glm::dvec2{initVelX, initVelY};
}

void BoundaryConditions::readSettings()
{
    double bboxMinX;
    double bboxMinY;
    double bboxMaxX;
    double bboxMaxY;

    std::ifstream Config_File("BCConfig.ini");
    options_description BCSet("Settings");
    BCSet.add_options()
        ("BoundaryConditions.bboxMinX",    value<double>(&bboxMinX))
        ("BoundaryConditions.bboxMinY",    value<double>(&bboxMinY))
        ("BoundaryConditions.bboxMaxX",    value<double>(&bboxMaxX))
        ("BoundaryConditions.bboxMaxY",    value<double>(&bboxMaxY))
        ("BoundaryConditions.boundaryCoeff", value<double>(&bndCoeff))
        ("BoundaryConditions.boundaryRange", value<double>(&bndConditionRange));

    variables_map vm;
    store(parse_config_file(Config_File, BCSet), vm);
    notify(vm);

    //process...
    bndBox = BoundingBox<glm::dvec2>(glm::dvec2(bboxMinX, bboxMinY),
                                     glm::dvec2(bboxMaxX, bboxMaxY));
}

void RenderSettings::readSettings()
{
    std::string fileRenderMode;
    std::string displayRenderMode;

    std::ifstream Config_File("RenderConfig.ini");
    options_description RenderSet("Settings");
    RenderSet.add_options()
        ("Rendering.width",      value<unsigned>(&width))
        ("Rendering.height",     value<unsigned>(&height))
        ("Rendering.printEvery", value<int>(&printEvr))
        ("Rendering.renderToFileAs",     value<std::string>(&fileRenderMode))
        ("Rendering.initRenderToScreen", value<std::string>(&displayRenderMode));

    variables_map vm;
    store(parse_config_file(Config_File, RenderSet), vm);
    notify(vm);

    //process...
    if (fileRenderMode == "RAWDATA") {
      fileRender= RAWDATA;
    }
    else if (fileRenderMode == "DISCRETE") {
      fileRender = DISCRETE;
    }
    else if (fileRenderMode == "METABALL") {
      fileRender = METABALL;
    }
    else {
      std::string errorMsg = "Unrecognized fileRender mode - Available modes::\n";
      errorMsg += "  RAWDATA\n";
      errorMsg += "  DISCRETE\n";
      errorMsg += "  METABALL\n";
      throw std::runtime_error(errorMsg);
    }

    if (displayRenderMode == "INDEX") {
      displayRender= INDEX;
    }
    else if (displayRenderMode == "SIMPLE") {
      displayRender= SIMPLE;
    }
    else if (displayRenderMode == "PRESSFORCES") {
      displayRender= PRESSFORCES;
    }
    else if (displayRenderMode == "VISCFORCES") {
      displayRender= VISCFORCES;
    }
    else if (displayRenderMode == "SURFFORCES") {
      displayRender= SURFFORCES;
    }
    else if (displayRenderMode == "OTHERFORCES") {
      displayRender= OTHERFORCES;
    }
    else if (displayRenderMode == "ALLFORCES") {
      displayRender= ALLFORCES;
    }
    else {
      std::string errorMsg = "Unrecognized initial display Render mode - Available modes::\n";
      errorMsg += "  INDEX\n";
      errorMsg += "  SIMPLE\n";
      errorMsg += "  PRESSFORCES\n";
      errorMsg += "  VISCFORCES\n";
      errorMsg += "  SURFFORCES\n";
      errorMsg += "  OTHERFORCES\n";
      errorMsg += "  ALLFORCES\n";
      throw std::runtime_error(errorMsg);
    }

    displayBox = BoundingBox<glm::vec2> (
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
}
