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
    static int LParticles;
    static int NParticles;
          
    static double initDx;

    static double particleDensity;
    static double dParticleDensity;
    static double particleMass;
    static double stiffness;
    static double densityErr;
    static double viscosity;
    static glm::dvec2 grav;
};

class SimulationSettings
{
public:
    static double dt;
    static double simTimeEnd;
    static double simTime;

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
        ALLIN,
        FAUCET,
        DRIPPING
    };
    static ParticleGeneration particleGeneration;
    static glm::dvec2 particleInitPos;
    static glm::dvec2 particleInitVel;
    static double     particleGenTime;
};

class BoundaryConditions
{
public:
    //background grid
    static BoundingBox<glm::dvec2> bndBox;
    static double bndCoeff;
    static double bndConditionRange;
};

//rendering settings
class RenderSettings
{
public:
    enum FileRenderType
    {
        GNUPLOT,
        DISCRETE,
        METABALL
    };
    enum DisplayRenderType
    {
        SIMPLE,
        PRESSFORCES,
        VISCFORCES,
        OTHERFORCES,
        ALLFORCES
    };
    static BoundingBox<glm::vec2> displayBox;
    //static glm::vec3 displayMinPos(float(bndBox.minX()-0.1*bndBox.dx()),  float(   bndBox.minY()-0.1*bndBox.dx()),  0.f);
    //static glm::vec3 displayMaxPos(float(bndBox.maxX()+0.1*bndBox.dx()),  float(2.*bndBox.maxX()+0.1*bndBox.dx()), 10.f);
    static int width;
    static int height;
    static int printEvr;
    static FileRenderType fileRender;
    static DisplayRenderType displayRender;
};

#endif
