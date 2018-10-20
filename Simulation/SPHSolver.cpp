
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
#include "Settings.hpp"
#include "Statistics/Statistics.hpp"

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
                glm::dvec2 rij = cntParticle-offParticle;
                double dist = glm::length(rij);
                glm::dvec2 tmpGradWij = Kernel::poly6::gradW(rij,dist)
                                       *Kernel::poly6::gradW_coeff();
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
                if (activeParticles_>=SPHSettings::NParticles) break;
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

    neibhs_.clear();
    neibhs_.findNei(cloud_, activeParticles_);
    updateNei();
    calcDensity();
}

//********************************************************************************
void SPHSolver::WCSPHStep() 
//********************************************************************************
{
    static auto densCalcTimerID   = Statistics::createTimer("SPHSolver::Step::WCSPH::densCalcTimer  ");
    static auto pressCalcTimerID  = Statistics::createTimer("SPHSolver::Step::WCSPH::pressCalcTimer ");
    static auto colorCalcTimerID  = Statistics::createTimer("SPHSolver::Step::WCSPH::colorCalcTimer ");
    static auto pressForceTimerID = Statistics::createTimer("SPHSolver::Step::WCSPH::pressForceTimer");
    static auto viscForceTimerID  = Statistics::createTimer("SPHSolver::Step::WCSPH::viscForceTimer ");
    static auto surfForceTimerID  = Statistics::createTimer("SPHSolver::Step::WCSPH::surfForceTimer ");
    static auto otherForceTimerID = Statistics::createTimer("SPHSolver::Step::WCSPH::otherForceTimer");
    static auto updatePosTimerID  = Statistics::createTimer("SPHSolver::Step::WCSPH::updatePosTimer ");

    auto densCalculator       = Statistics::CountTime(densCalcTimerID,   [this](){ this->calcDensity();  });
    auto pressCalculator      = Statistics::CountTime(pressCalcTimerID,  [this](){ this->calcPressure(); });
    auto colorCalculator      = Statistics::CountTime(colorCalcTimerID,  [this](){ this->calcNormal();   });
    auto pressForceCalculator = Statistics::CountTime(pressForceTimerID, [this](){ this->calcPressForces(); });
    auto viscForceCalculator  = Statistics::CountTime(viscForceTimerID,  [this](){ this->calcViscForces();  });
    auto surfForceCalculator  = Statistics::CountTime(surfForceTimerID,  [this](){ this->calcSurfForces();  });
    auto otherForceCalculator = Statistics::CountTime(otherForceTimerID, [this](){ this->calcOtherForces(); });

    densCalculator();
    pressCalculator();
    colorCalculator();

    //calc forces
    pressForceCalculator();
    viscForceCalculator();
    surfForceCalculator();
    otherForceCalculator();

    //combine forces and update values
    {
        Statistics::TimerGuard g(updatePosTimerID);
        #pragma omp parallel for
        for (int iPart = 0;iPart<activeParticles_;iPart++) 
        {
            cloud_[iPart].Ftot = cloud_[iPart].Fpress
                               + cloud_[iPart].Fvisc
                               + cloud_[iPart].Fsurf
                               + cloud_[iPart].Fother;

            cloud_[iPart].velocity += SimulationSettings::dt*cloud_[iPart].ddensity*cloud_[iPart].Ftot;
            cloud_[iPart].position += SimulationSettings::dt*cloud_[iPart].velocity;
        }
    }
}

//********************************************************************************
void SPHSolver::PCISPHStep() 
//********************************************************************************
{
    static auto densCalcTimerID   = Statistics::createTimer("SPHSolver::Step::PCISPH::densCalcTimer  ");
    static auto colorCalcTimerID  = Statistics::createTimer("SPHSolver::Step::PCISPH::colorCalcTimer ");
    static auto viscForceTimerID  = Statistics::createTimer("SPHSolver::Step::PCISPH::viscForceTimer ");
    static auto surfForceTimerID  = Statistics::createTimer("SPHSolver::Step::PCISPH::surfForceTimer ");
    static auto otherForceTimerID = Statistics::createTimer("SPHSolver::Step::PCISPH::otherForceTimer");
    static auto pressForceTimerID = Statistics::createTimer("SPHSolver::Step::PCISPH::pressForceTimer");
    static auto updatePosTimerID  = Statistics::createTimer("SPHSolver::Step::PCISPH::updatePosTimer ");

    auto densCalculator       = Statistics::CountTime(densCalcTimerID,   [this](){ this->calcDensity(); });
    auto colorCalculator      = Statistics::CountTime(colorCalcTimerID,  [this](){ this->calcNormal();  });
    auto viscForceCalculator  = Statistics::CountTime(viscForceTimerID,  [this](){ this->calcViscForces();  });
    auto surfForceCalculator  = Statistics::CountTime(surfForceTimerID,  [this](){ this->calcSurfForces();  });
    auto otherForceCalculator = Statistics::CountTime(otherForceTimerID, [this](){ this->calcOtherForces(); });

    densCalculator();
    colorCalculator();
    
    //calc forces
    viscForceCalculator();
    surfForceCalculator();
    otherForceCalculator();

    {
        Statistics::TimerGuard g(updatePosTimerID);
        #pragma omp parallel for
        for (int iPart = 0;iPart<activeParticles_;iPart++) 
        {
            cloud_[iPart].Ftot = cloud_[iPart].Fvisc+cloud_[iPart].Fsurf+cloud_[iPart].Fother;

            cloud_[iPart].velocity += SimulationSettings::dt*cloud_[iPart].ddensity*cloud_[iPart].Ftot;
            cloud_[iPart].position += SimulationSettings::dt*cloud_[iPart].velocity;
        }
    }

    double densErr=SPHSettings::particleDensity;
    int iter=0;

    {
        //calculate pressure and pressure force
        Statistics::TimerGuard g(pressForceTimerID);
        initPressure();
        const double targetDensErr = SPHSettings::densityErr*SPHSettings::particleDensity;
        while(densErr>targetDensErr && iter<500)
        {
            updateNei();
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
                iParticle.Ftot+= iParticle.Fpress;

                glm::dvec2 update = SimulationSettings::dt*iParticle.ddensity*iParticle.Fpress;
                iParticle.velocity += update;
                iParticle.position += SimulationSettings::dt*update;
            }
        }
    }
}

//********************************************************************************
bool SPHSolver::step() 
//********************************************************************************
{
    static auto reorTimerID      = Statistics::createTimer("SPHSolver::Step::reorTime     ");
    static auto findNeiTimerID   = Statistics::createTimer("SPHSolver::Step::findNeiTime  ");
    static auto updateNeiTimerID = Statistics::createTimer("SPHSolver::Step::updateNeiTime");
    static auto stepTimerID      = Statistics::createTimer("SPHSolver::Step::SPHsolver    ");

    static int iReorder = 0;
    iReorder++;

    auto findNeighbours   = Statistics::CountTime(findNeiTimerID,   &HashTable::findNei);
    auto updateNeighbours = Statistics::CountTime(updateNeiTimerID, [this](){ this->updateNei();});

    #pragma omp parallel for
    for (int i=0;i<activeParticles_;i++)
    {
        if(cloud_[i].nei.size())
            cloud_[i].nei.clear();
    }

    if(iReorder%100==0)
    {
        Statistics::CountTime(reorTimerID, &HashTable::reorderCloud)(
          neibhs_, cloud_, activeParticles_
        );
    }

    generateParticles();

    neibhs_.clear();

    findNeighbours(neibhs_, cloud_, activeParticles_);

    updateNeighbours();

    {
        Statistics::TimerGuard g(stepTimerID);
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
    }

    return (!SimulationSettings::breakLoop());
}

//********************************************************************************
void SPHSolver::generateParticles()
//********************************************************************************
{
    static double genTime = 0;
    static double dtGenFaucet = 0.51*Kernel::SmoothingLength::h
                     /(glm::length(InitialConditions::particleInitVel)+1.e-8);

    if (activeParticles_>=SPHSettings::NParticles) return;
     
    switch(InitialConditions::particleGeneration)
    {
        case (InitialConditions::FAUCET):
        {
            if (SimulationSettings::simTime-genTime>=dtGenFaucet)
            {
                glm::dvec2 initVelNormalized = glm::normalize(InitialConditions::particleInitVel);
                glm::dvec2 generationDirection = glm::dvec2(-initVelNormalized.y,initVelNormalized.x);
                //std::cout<<"generating"<<std::endl;
                if (activeParticles_<SPHSettings::NParticles)
                {
                    for (int i=0;i<SPHSettings::LParticles;i++)
                    {
                        if (activeParticles_>=SPHSettings::NParticles) break;
                        cloud_[activeParticles_].position = InitialConditions::particleInitPos + i*SPHSettings::initDx*generationDirection;
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
                    //double     radius = 0.3*double(SPHSettings::LParticles)*SPHSettings::initDx; 
                    double     radius = 0.4*double(SPHSettings::LParticles)*SPHSettings::initDx; 
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
void SPHSolver::updateNei()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        Particle& iParticle = cloud_[iPart];
        int Nnei = (int)iParticle.nei.size();
        for (int i = 0; i<Nnei;i++)
        {
            int jPart = cloud_[iPart].nei[i].ID;
            if (iPart!=jPart)
            {
                Particle& jParticle = cloud_[jPart];

                iParticle.nei[i].dir  = iParticle.position-jParticle.position;
                iParticle.nei[i].dist = glm::length(iParticle.nei[i].dir);
            }
        }
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
            int jPart = iParticle.nei[i].ID;
            Particle& jParticle = cloud_[jPart];

            dens += jParticle.mass*Kernel::poly6::W(iParticle.nei[i].dist);
        }
        dens*=Kernel::poly6::W_coeff();

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
void SPHSolver::calcNormal()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        Particle& iParticle = cloud_[iPart];

        glm::dvec2 norm (0.0);

        int Nnei = (int)iParticle.nei.size();
        for (int i = 0; i<Nnei;i++)
        {
            int jPart = iParticle.nei[i].ID;
            Particle& jParticle = cloud_[jPart];

            norm += jParticle.mass*jParticle.ddensity
                  *Kernel::poly6::gradW(iParticle.nei[i].dir,iParticle.nei[i].dist);
        }
        norm *= Kernel::poly6::gradW_coeff()*Kernel::SmoothingLength::h;
        iParticle.normal = norm;
    }
}


//********************************************************************************
void SPHSolver::updatePressure()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        cloud_[iPart].pressure = fmax(delta_*cloud_[iPart].densityErr,0.0);
      //cloud_[iPart].pressure = delta_*cloud_[iPart].densityErr;
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
        double W0 = Kernel::poly6::W(0.0)*Kernel::poly6::W_coeff();
        double mult = 2.0;
        double scaledh  = mult*Kernel::SmoothingLength::h;
        double dScaledh = 1./scaledh;
        if (iPos.x-BoundaryConditions::bndBox.minX()<scaledh)
        {
            if (iPos.x>BoundaryConditions::bndBox.minX())
            {
                tmpiPos = glm::dvec2(iPos.x,0.0)*dScaledh;
                bndPos  = glm::dvec2(BoundaryConditions::bndBox.minX(),0.0)*dScaledh;
                double dist = glm::length(tmpiPos-bndPos);
                cloud_[iPart].Fother.x+=BoundaryConditions::bndCoeff
                                       *Kernel::poly6::W(dist)
                                       *Kernel::poly6::W_coeff();
            }
            else
            {
                cloud_[iPart].velocity.x=0.;
                cloud_[iPart].Fother.x+=BoundaryConditions::bndCoeff*W0;
            }
        }
        else if (BoundaryConditions::bndBox.maxX()-iPos.x<scaledh)
        {
            if (iPos.x<BoundaryConditions::bndBox.maxX())
            {
                tmpiPos = glm::dvec2(iPos.x,0.0)*dScaledh;
                bndPos  = glm::dvec2(BoundaryConditions::bndBox.maxX(),0.0)*dScaledh;
                double dist = glm::length(tmpiPos-bndPos);
                cloud_[iPart].Fother.x-=BoundaryConditions::bndCoeff
                                       *Kernel::poly6::W(dist)
                                       *Kernel::poly6::W_coeff();
            }
            else
            {
                cloud_[iPart].velocity.x=0.;
                cloud_[iPart].Fother.x-=BoundaryConditions::bndCoeff*W0;
            }
        }
        if (iPos.y-BoundaryConditions::bndBox.minY()<scaledh)
        {
            if (iPos.y>BoundaryConditions::bndBox.minY())
            {
                tmpiPos = glm::dvec2(0.0,iPos.y)*dScaledh;
                bndPos  = glm::dvec2(0.0,BoundaryConditions::bndBox.minY())*dScaledh;
                double dist = glm::length(tmpiPos-bndPos);
                cloud_[iPart].Fother.y+=BoundaryConditions::bndCoeff
                                       *Kernel::poly6::W(dist)
                                       *Kernel::poly6::W_coeff();
            }
            else
            {
                cloud_[iPart].velocity.y=0.;
                cloud_[iPart].Fother.y+=BoundaryConditions::bndCoeff*W0;
            }
        }
        else if (BoundaryConditions::bndBox.maxY()-iPos.y<scaledh)
        {
            if (iPos.y>BoundaryConditions::bndBox.maxY())
            {
                tmpiPos = glm::dvec2(0.0,iPos.y)*dScaledh;
                bndPos  = glm::dvec2(0.0,BoundaryConditions::bndBox.maxY())*dScaledh;
                double dist = glm::length(tmpiPos-bndPos);
                cloud_[iPart].Fother.y-=BoundaryConditions::bndCoeff
                                       *Kernel::poly6::W(dist)
                                       *Kernel::poly6::W_coeff();
            }
            else
            {
                cloud_[iPart].velocity.y=0.;
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

        int Nnei = (int)iParticle.nei.size();
        for (int i = 0; i<Nnei;i++)
        {
            int jPart = iParticle.nei[i].ID;
            if (iPart==jPart) continue;
            Particle& jParticle = cloud_[jPart];

            Fp+=jParticle.mass*jParticle.ddensity 
               *(iParticle.pressure + jParticle.pressure)
               *Kernel::spiky::gradW(iParticle.nei[i].dir,iParticle.nei[i].dist);
        }

        Fp*=-0.5*Kernel::spiky::gradW_coeff();
        
        iParticle.Fpress= Fp;
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

        int Nnei = (int)iParticle.nei.size();
        //std::cout<<"iPart= "<<iPart<<" nNei="<<Nnei;
        for (int i = 0; i<Nnei;i++)
        {
            int jPart = iParticle.nei[i].ID;
            //std::cout<<" "<<jPart;
            if (iPart==jPart) continue;
            Particle& jParticle = cloud_[jPart];

            Fv+= jParticle.mass
               *jParticle.ddensity
               *(jParticle.velocity-iParticle.velocity)
               *Kernel::visc::laplW(iParticle.nei[i].dist);
        }
        Fv*=SPHSettings::viscosity*Kernel::visc::laplW_coeff();

        iParticle.Fvisc= Fv;
    }
}

//********************************************************************************
void SPHSolver::calcSurfForces()
//********************************************************************************
{
    #pragma omp parallel for
    for (int iPart = 0;iPart<activeParticles_;iPart++) 
    {
        glm::dvec2 Fcohesion  = glm::dvec2(0.0);
        glm::dvec2 Fcurvature = glm::dvec2(0.0);
        double correction = 0.0;

        Particle& iParticle = cloud_[iPart];

        int Nnei = (int)iParticle.nei.size();
        //std::cout<<"iPart= "<<iPart<<" nNei="<<Nnei;
        for (int i = 0; i<Nnei;i++)
        {
            int jPart = iParticle.nei[i].ID;
            //std::cout<<" "<<jPart;
            if (iPart==jPart) continue;
            Particle& jParticle = cloud_[jPart];

            correction = 2.*SPHSettings::particleDensity/(iParticle.density+jParticle.density);

            Fcohesion+= correction
                   *iParticle.mass*jParticle.mass 
                   *Kernel::surface::C(iParticle.nei[i].dist)
                   *iParticle.nei[i].dir/iParticle.nei[i].dist;

            Fcurvature+= correction
                   *iParticle.mass
                   *(iParticle.normal - jParticle.normal);
        }

        Fcohesion *= Kernel::surface::C_coeff();

        iParticle.Fsurf = -SPHSettings::surfTension
                          * iParticle.density
                          * (Fcohesion+Fcurvature);
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
