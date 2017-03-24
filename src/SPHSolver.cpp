
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "SPHSolver.hpp"
#include "WriteFunctions.hpp"
#include "Settings.hpp"
#include "Statistics.hpp"

//********************************************************************************
void SPHSolver::init()
//********************************************************************************
{
    neibhs_.setHashTable(BoundaryConditions::bndBox.minPos(),BoundaryConditions::bndBox.delta());
    //neibhs_.writeGridRAW("neiGrid");

    cloud_.resize(SPHSettings::NParticles);

    if (SPHSettings::SPHstep==SPHSettings::Solver::PCISPH)
    {
        //calculate delta
        glm::dvec2 cntParticle(0.0);
        glm::dvec2 offParticle(0.0);
        glm::dvec2 gradWij(0.0);
        double     dotGradWij(0.0);
        for (int i=-2;i<=2;i++)
        {
            for (int j=-2;j<=2;j++)
            {
                if (i==0&&j==0) continue;
                offParticle = glm::dvec2(i*SPHSettings::initDx,j*SPHSettings::initDx);
                glm::dvec2 tmpGradWij = Kernel::poly6::gradW(cntParticle,offParticle);
                gradWij += tmpGradWij;
                dotGradWij += glm::dot(tmpGradWij,tmpGradWij);
            }
        }
        double beta = SPHSettings::particleMass * SimulationSettings::dt * SPHSettings::dParticleDensity;
        beta *= beta;
        beta *= 2.0;

        delta_ =  - 1./(beta*(-glm::dot(gradWij,gradWij)-dotGradWij));
    }

    if (InitialConditions::particleGeneration==InitialConditions::ALLIN)
    {
        //initialize field Particles
        int iPart = 0;
        int i = 0;
        //for (int i=0;i<YParticles;i++)
        while (iPart<SPHSettings::NParticles)
        {
            for (int j=0;j<SPHSettings::LParticles;j++)
            {
                //int iPart = i*XParticles+j;
                cloud_[iPart].position = glm::dvec2
                (
                    InitialConditions::particleInitPos.x+(SPHSettings::initDx*(double(j)+0.5)),
                    InitialConditions::particleInitPos.y+SPHSettings::initDx*(double(i)+0.5)
                ); 
                cloud_[iPart].velocity = InitialConditions::particleInitVel; 
                cloud_[iPart].density  = SPHSettings::particleDensity;
                cloud_[iPart].ddensity = SPHSettings::dParticleDensity;
                cloud_[iPart].mass     = SPHSettings::particleMass;
                iPart++;
                activeParticles_++;
            }
            i++;
        }
    }
    else
    {
        //particles will be generated in generateParticles()
    }

    neibhs_.findNei(cloud_, activeParticles_);
    calcDensity();

    if (RenderSettings::fileRender==RenderSettings::RAWDATA)
    {
        writeRAWfile("step", cloud_);
    }
    else
    {
        renderImage("render", cloud_, neibhs_);
    }
}

//********************************************************************************
void SPHSolver::WCSPHStep() 
//********************************************************************************
{
    static int densCalcTimerID   = Statistics::createTimer("SPHSolver::Step::WCSPH::densCalcTimer  ");
    static int pressCalcTimerID  = Statistics::createTimer("SPHSolver::Step::WCSPH::pressCalcTimer ");
    static int pressForceTimerID = Statistics::createTimer("SPHSolver::Step::WCSPH::pressForceTimer");
    static int viscForceTimerID  = Statistics::createTimer("SPHSolver::Step::WCSPH::viscForceTimer ");
    static int otherForceTimerID = Statistics::createTimer("SPHSolver::Step::WCSPH::otherForceTimer");
    static int updatePosTimerID  = Statistics::createTimer("SPHSolver::Step::WCSPH::updatePosTimer ");

    COUNT_TIME(calcDensity(), densCalcTimerID);
    COUNT_TIME(calcPressure(),pressCalcTimerID);

    //calc forces
    COUNT_TIME(calcPressForces(), pressForceTimerID);
    COUNT_TIME(calcViscForces(),  viscForceTimerID);
    COUNT_TIME(calcOtherForces(), otherForceTimerID);

    //combine forces and update values
    Statistics::timers[updatePosTimerID].start();
    #pragma omp parallel for
    for (int iPart = 0;iPart<activeParticles_;iPart++) 
    {
        cloud_[iPart].Ftot = cloud_[iPart].Fpress+cloud_[iPart].Fvisc+cloud_[iPart].Fother;

        cloud_[iPart].velocity += SimulationSettings::dt*cloud_[iPart].ddensity*cloud_[iPart].Ftot;
        cloud_[iPart].position += SimulationSettings::dt*cloud_[iPart].velocity;
    }
    Statistics::timers[updatePosTimerID].end();
    Statistics::timers[updatePosTimerID].addTime();
}

//********************************************************************************
void SPHSolver::PCISPHStep() 
//********************************************************************************
{
    static int densCalcTimerID   = Statistics::createTimer("SPHSolver::Step::PCISPH::densCalcTimer  ");
    static int viscForceTimerID  = Statistics::createTimer("SPHSolver::Step::PCISPH::viscForceTimer ");
    static int otherForceTimerID = Statistics::createTimer("SPHSolver::Step::PCISPH::otherForceTimer");
    static int pressForceTimerID = Statistics::createTimer("SPHSolver::Step::PCISPH::pressForceTimer");
    static int updatePosTimerID  = Statistics::createTimer("SPHSolver::Step::PCISPH::updatePosTimer ");

    COUNT_TIME(calcDensity(), densCalcTimerID);
    
    //calc forces
    COUNT_TIME(calcViscForces(),  viscForceTimerID);
    COUNT_TIME(calcOtherForces(), otherForceTimerID);

    Statistics::timers[updatePosTimerID].start();
    #pragma omp parallel for
    for (int iPart = 0;iPart<activeParticles_;iPart++) 
    {
        cloud_[iPart].Ftot = cloud_[iPart].Fvisc+cloud_[iPart].Fother;

        cloud_[iPart].velocity += SimulationSettings::dt*cloud_[iPart].ddensity*cloud_[iPart].Ftot;
        cloud_[iPart].position += SimulationSettings::dt*cloud_[iPart].velocity;
    }
    Statistics::timers[updatePosTimerID].end();
    Statistics::timers[updatePosTimerID].addTime();

    double densErr=SPHSettings::particleDensity;
    int iter=0;

    //calculate pressure and pressure force
    Statistics::timers[pressForceTimerID].start();
    initPressure();
    const double targetDensErr = SPHSettings::densityErr*SPHSettings::particleDensity;
    while(densErr>targetDensErr && iter<500)
    {
        calcDensity();
        calcDensityErr();
        updatePressure();

        densErr = 0.;
        for (int iPart = 0;iPart<activeParticles_;iPart++) 
        {
            if (densErr<cloud_[iPart].densityErr)
                densErr = cloud_[iPart].densityErr;
        }

        if (iter==0)
        iter++;

        calcPressForces();
        //#pragma omp parallel for
        for (int iPart = 0;iPart<activeParticles_;iPart++) 
        {
            Particle& iParticle = cloud_[iPart];
            iParticle.Ftot = iParticle.Fpress;

            glm::dvec2 update = SimulationSettings::dt*iParticle.ddensity*iParticle.Ftot;
            iParticle.velocity += update;
            iParticle.position += SimulationSettings::dt*update;
        }
    }
    Statistics::timers[pressForceTimerID].end();
    Statistics::timers[pressForceTimerID].addTime();
}

//********************************************************************************
bool SPHSolver::step() 
//********************************************************************************
{
    static int neibTimerID = Statistics::createTimer("SPHSolver::Step::neibTime");
    static int stepTimerID = Statistics::createTimer("SPHSolver::Step::SPHsolver");

    generateParticles();

    COUNT_TIME(neibhs_.findNei(cloud_, activeParticles_), neibTimerID);

    Statistics::timers[stepTimerID].start();
    switch (SPHSettings::SPHstep)
    {
        case (SPHSettings::Solver::WCSPH):
        {
            WCSPHStep();
            break;
        }
        case (SPHSettings::Solver::PCISPH):
        {
            PCISPHStep();
            break;
        }
        default:
        {
            std::cerr<<"Invalid SPHSettings::Solver"<<std::endl;
            exit(1);
        }
    }
    SimulationSettings::updateSimTime();

    Statistics::timers[stepTimerID].end();
    Statistics::timers[stepTimerID].addTime();

    return (!SimulationSettings::breakLoop());
}

//********************************************************************************
void SPHSolver::generateParticles()
//********************************************************************************
{
    static double genTime = 0;
    static double dtGenFaucet = 0.5*Kernel::SmoothingLength::h
                     /(glm::length(InitialConditions::particleInitVel)+1.e-8);

    if (activeParticles_>=SPHSettings::NParticles) return;
     
    switch(InitialConditions::particleGeneration)
    {
        case (InitialConditions::FAUCET):
        {
            if (SimulationSettings::simTime-genTime>=dtGenFaucet)
            {
                //std::cout<<"generating"<<std::endl;
                if (activeParticles_<SPHSettings::NParticles)
                {
                    for (int i=0;i<SPHSettings::LParticles;i++)
                    {
                        if (activeParticles_>=SPHSettings::NParticles) break;
                        cloud_[activeParticles_].position = InitialConditions::particleInitPos + glm::dvec2(0,SPHSettings::initDx*i);
                        cloud_[activeParticles_].velocity = InitialConditions::particleInitVel;
                        cloud_[activeParticles_].mass     = SPHSettings::particleMass;
                        activeParticles_++;
                    }
                }
                genTime = SimulationSettings::simTime;
            }
            break;
        }
        case (InitialConditions::DRIPPING):
        {
            if (SimulationSettings::simTime-genTime>=InitialConditions::particleGenTime)
            {
                //std::cout<<"generating"<<std::endl;
                if (activeParticles_<SPHSettings::NParticles)
                {
                    double     radius = 0.3*double(SPHSettings::LParticles)*SPHSettings::initDx; 
                    double     offset = 0.5*SPHSettings::initDx*double(SPHSettings::LParticles-1);
                    glm::dvec2 center = InitialConditions::particleInitPos + glm::dvec2(offset,offset);

                    for (int i=0;i<SPHSettings::LParticles;i++)
                    {
                        for (int j=0;j<SPHSettings::LParticles;j++)
                        {
                            if (activeParticles_>=SPHSettings::NParticles) break;
                            glm::dvec2 pos(InitialConditions::particleInitPos + glm::dvec2(SPHSettings::initDx*i,SPHSettings::initDx*j));
                            double dist = glm::length(pos-center);
                            if (dist>radius) continue;
                            cloud_[activeParticles_].position = pos;
                            cloud_[activeParticles_].velocity = InitialConditions::particleInitVel;
                            cloud_[activeParticles_].mass     = SPHSettings::particleMass;
                            activeParticles_++;
                        }
                    }
                }
                genTime = SimulationSettings::simTime;
            }
            break;
        }
        case (InitialConditions::ALLIN):
        {
            std::cerr<<"SPHSolver::generateParticles::ALLIN - Should be in here!!"<<std::endl;
            exit(1);
        }
        default: {}
    }
}

//********************************************************************************
void SPHSolver::calcDensity()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        Particle& iParticle = cloud_[iPart];

        double dens = 0.0;

        int Nnei = (int)iParticle.nei.size();
        for (int i = 0; i<Nnei;i++)
        {
            int neiPart = iParticle.nei[i];
            Particle& neiParticle = cloud_[neiPart];

            double Wij = Kernel::poly6::W(iParticle.position,neiParticle.position);
            dens += neiParticle.mass*Wij;
        }
        iParticle.density = dens;
        iParticle.ddensity = 1./dens;
    }
}

//********************************************************************************
void SPHSolver::calcDensityErr()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        cloud_[iPart].densityErr = cloud_[iPart].density - SPHSettings::particleDensity;
    }
}

//********************************************************************************
void SPHSolver::initPressure()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        cloud_[iPart].pressure = 0.0;
    }
}

//********************************************************************************
void SPHSolver::calcPressure()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        //p = k(rho-rho0)
        cloud_[iPart].pressure = fmax(SPHSettings::stiffness*(cloud_[iPart].density-SPHSettings::particleDensity),0.0);
      //cloud_[iPart].pressure =      SPHSettings::stiffness*(cloud_[iPart].density-SPHSettings::particleDensity);
    }
}

//********************************************************************************
void SPHSolver::updatePressure()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        cloud_[iPart].pressure = 1.*fmax(delta_*cloud_[iPart].densityErr,0.0);
    }
}


//********************************************************************************
void SPHSolver::calcOtherForces()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        //gravity
        cloud_[iPart].Fother = SPHSettings::grav;

        //boundary forces
        glm::dvec2 iPos = cloud_[iPart].position;
        glm::dvec2 bndPos(0.0,0.0);
        glm::dvec2 tmpiPos(0.0,0.0);
        double W0 = Kernel::poly6::W(tmpiPos,tmpiPos);
        double mult = 2.0;
        double scaledh  = mult*Kernel::SmoothingLength::h;
        double dScaledh = 1./scaledh;
        if (iPos.x-BoundaryConditions::bndBox.minX()<scaledh)
        {
            if (iPos.x>BoundaryConditions::bndBox.minX())
            {
                tmpiPos = glm::dvec2(iPos.x,0.0)*dScaledh;
                bndPos  = glm::dvec2(BoundaryConditions::bndBox.minX(),0.0)*dScaledh;
                cloud_[iPart].Fother.x+=BoundaryConditions::bndCoeff*Kernel::poly6::W(tmpiPos,bndPos);
            }
            else
            {
                //cloud_[iPart].velocity.x=0.;
                cloud_[iPart].Fother.x+=BoundaryConditions::bndCoeff*W0;
            }
        }
        else if (BoundaryConditions::bndBox.maxX()-iPos.x<scaledh)
        {
            if (iPos.x<BoundaryConditions::bndBox.maxX())
            {
                tmpiPos = glm::dvec2(iPos.x,0.0)*dScaledh;
                bndPos  = glm::dvec2(BoundaryConditions::bndBox.maxX(),0.0)*dScaledh;
                cloud_[iPart].Fother.x-=BoundaryConditions::bndCoeff*Kernel::poly6::W(tmpiPos,bndPos);
            }
            else
            {
                //cloud_[iPart].velocity.x=0.;
                cloud_[iPart].Fother.x-=BoundaryConditions::bndCoeff*W0;
            }
        }
        if (iPos.y-BoundaryConditions::bndBox.minY()<scaledh)
        {
            if (iPos.y>BoundaryConditions::bndBox.minY())
            {
                tmpiPos = glm::dvec2(0.0,iPos.y)*dScaledh;
                bndPos  = glm::dvec2(0.0,BoundaryConditions::bndBox.minY())*dScaledh;
                cloud_[iPart].Fother.y+=BoundaryConditions::bndCoeff*Kernel::poly6::W(tmpiPos,bndPos);
            }
            else
            {
                //cloud_[iPart].velocity.y=0.;
                cloud_[iPart].Fother.y+=BoundaryConditions::bndCoeff*W0;
            }
        }
        else if (BoundaryConditions::bndBox.maxY()-iPos.y<scaledh)
        {
            if (iPos.y>BoundaryConditions::bndBox.maxY())
            {
                tmpiPos = glm::dvec2(0.0,iPos.y)*dScaledh;
                bndPos  = glm::dvec2(0.0,BoundaryConditions::bndBox.maxY())*dScaledh;
                cloud_[iPart].Fother.y-=BoundaryConditions::bndCoeff*Kernel::poly6::W(tmpiPos,bndPos);
            }
            else
            {
                //cloud_[iPart].velocity.y=0.;
                cloud_[iPart].Fother.y-=BoundaryConditions::bndCoeff*W0;
            }
        }

        cloud_[iPart].Fother*= cloud_[iPart].density;
    }
}

//********************************************************************************
void SPHSolver::calcPressForces()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart = 0;iPart<activeParticles_;iPart++) 
    {
        Particle& iParticle = cloud_[iPart];

        glm::dvec2 Fp = glm::dvec2(0.0);

        int Nnei = (int)cloud_[iPart].nei.size();
        for (int i = 0; i<Nnei;i++)
        {
            int neiPart = cloud_[iPart].nei[i];
            if (iPart==neiPart) continue;
            Particle& neiParticle = cloud_[neiPart];

            Fp+=neiParticle.mass
               *neiParticle.ddensity 
               *(iParticle.pressure + neiParticle.pressure)
               *Kernel::spiky::gradW(iParticle.position,neiParticle.position);
        }

        Fp*=-0.5;
        
        cloud_[iPart].Fpress= Fp;
    }
}

//********************************************************************************
void SPHSolver::calcViscForces()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart = 0;iPart<activeParticles_;iPart++) 
    {
        glm::dvec2 Fv = glm::dvec2(0.0);

        Particle& iParticle = cloud_[iPart];

        int Nnei = (int)cloud_[iPart].nei.size();
        //std::cout<<"iPart= "<<iPart<<" nNei="<<Nnei;
        for (int i = 0; i<Nnei;i++)
        {
            int neiPart = cloud_[iPart].nei[i];
            //std::cout<<" "<<neiPart;
            if (iPart==neiPart) continue;
            Particle& neiParticle = cloud_[neiPart];

            Fv+= neiParticle.mass
               *neiParticle.ddensity
               *(neiParticle.velocity-iParticle.velocity)
               *Kernel::visc::laplW(iParticle.position,neiParticle.position);
        }
        Fv*=SPHSettings::viscosity;

        cloud_[iPart].Fvisc= Fv;
    }
}

//********************************************************************************
double SPHSolver::calcCFL()
//********************************************************************************
{
    double vMax = 0.0;

    for (int iPart = 0;iPart<activeParticles_;iPart++) 
    {
        double vLen2 = glm::length2(cloud_[iPart].velocity);
        if (vLen2 > vMax) vMax = vLen2;
    }

    double maxDt = 0.4*Kernel::SmoothingLength::h/(sqrt(vMax)+1.e-5);
    double CFL = SimulationSettings::dt/maxDt;
    return CFL;
}
