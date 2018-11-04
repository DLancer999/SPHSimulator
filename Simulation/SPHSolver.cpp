
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
#include <numeric>

#include "SPHSolver.hpp"
#include "Settings.hpp"
#include "Kernels.hpp"
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

    calcDensity();
    calcPressure();
    calcNormal();

    //calc forces
    calcPressForces();
    calcViscForces();
    calcSurfForces();
    calcOtherForces();

    //combine forces and update values
    {
        static auto updatePosTimerID  = Statistics::createTimer("SPHSolver::WCSPH::updatePosTimer");
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
    calcDensity();
    calcNormal();
    
    //calc forces
    calcViscForces();
    calcSurfForces();
    calcOtherForces();

    {
        static auto updatePosTimerID  = Statistics::createTimer("SPHSolver::PCISPH::updatePosTimer");
        Statistics::TimerGuard g(updatePosTimerID);

        #pragma omp parallel for
        for (int iPart = 0;iPart<activeParticles_;iPart++) 
        {
            Particle& particle = cloud_[iPart];
            particle.Ftot = particle.Fvisc+particle.Fsurf+particle.Fother;

            particle.velocity += SimulationSettings::dt*particle.ddensity*particle.Ftot;
            particle.position += SimulationSettings::dt*particle.velocity;
        }
    }

    double densErr=SPHSettings::particleDensity;
    int iter=0;

    {
        static auto pressForceTimerID = Statistics::createTimer("SPHSolver::PCISPH::pressForceTimer");
        Statistics::TimerGuard g(pressForceTimerID);

        //calculate pressure and pressure force
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
    static int iReorder = 0;
    iReorder++;

    #pragma omp parallel for
    for (int i=0;i<activeParticles_;i++)
    {
        if(cloud_[i].nei.size())
            cloud_[i].nei.clear();
    }

    if(iReorder%100==0)
    {
        neibhs_.reorderCloud(cloud_,activeParticles_);
    }

    generateParticles();

    neibhs_.clear();

    neibhs_.findNei(cloud_, activeParticles_);

    updateNei();

    {
        static auto stepTimerID      = Statistics::createTimer("SPHSolver::Step::SPHsolver");
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
                        Particle& active = cloud_[activeParticles_];
                        active.position = InitialConditions::particleInitPos + i*SPHSettings::initDx*generationDirection;
                        active.velocity = InitialConditions::particleInitVel;
                        active.mass     = SPHSettings::particleMass;
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
                            Particle& active = cloud_[activeParticles_];
                            active.position = pos;
                            active.velocity = InitialConditions::particleInitVel;
                            active.mass     = SPHSettings::particleMass;
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
    static auto updateNeiTimerID = Statistics::createTimer("SPHSolver::updateNeiTime");
    Statistics::TimerGuard updateNeiTimerGuard(updateNeiTimerID);

    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        Particle& iParticle = cloud_[iPart];
        int Nnei = (int)iParticle.nei.size();
        for (int i = 0; i<Nnei;i++)
        {
            int jPart = cloud_[iPart].nei[i].ID;
            Neigbhor& iPartNeiI = iParticle.nei[i];
            if (iPart!=jPart)
            {
                Particle& jParticle = cloud_[jPart];

                iPartNeiI.dir  = iParticle.position-jParticle.position;
                iPartNeiI.dist = glm::length(iPartNeiI.dir);
            }
        }
    }
}

//********************************************************************************
void SPHSolver::calcDensity()
//********************************************************************************
{
    static auto densCalcTimerID   = Statistics::createTimer("SPHSolver::densCalcTimer");
    Statistics::TimerGuard densGuard(densCalcTimerID);

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
    static auto densErrCalcTimerID   = Statistics::createTimer("SPHSolver::densErrCalcTimer");
    Statistics::TimerGuard densErrGuard(densErrCalcTimerID);

    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        auto& particle = cloud_[iPart];
        particle.densityErr = particle.density - SPHSettings::particleDensity;
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
    static auto pressCalcTimerID   = Statistics::createTimer("SPHSolver::pressureCalcTimer");
    Statistics::TimerGuard pressGuard(pressCalcTimerID);

    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        auto& particle = cloud_[iPart];
        //p = k(rho-rho0)
        particle.pressure = fmax(SPHSettings::stiffness*(particle.density-SPHSettings::particleDensity),0.0);
    }
}

//********************************************************************************
void SPHSolver::calcNormal()
//********************************************************************************
{
    static auto normalCalcTimerID   = Statistics::createTimer("SPHSolver::normalCalcTimer");
    Statistics::TimerGuard normalGuard(normalCalcTimerID);

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
    static auto updatePressTimerID   = Statistics::createTimer("SPHSolver::updatePressure");
    Statistics::TimerGuard updatePressGuard(updatePressTimerID);

    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        auto& particle = cloud_[iPart];
        particle.pressure = fmax(delta_*particle.densityErr,0.0);
      //cloud_[iPart].pressure = delta_*cloud_[iPart].densityErr;
    }
}


//********************************************************************************
void SPHSolver::calcOtherForces()
//********************************************************************************
{
    static auto otherForceTimerID = Statistics::createTimer("SPHSolver::otherForceTimer");
    Statistics::TimerGuard otherForceGuard(otherForceTimerID);

    #pragma omp parallel for
    for (int iPart=0;iPart<activeParticles_;iPart++)
    {
        auto& particle = cloud_[iPart];
        //gravity
        particle.Fother = SPHSettings::grav;

        //boundary forces
        glm::dvec2 iPos = particle.position;
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
                particle.Fother.x+=BoundaryConditions::bndCoeff
                                  *Kernel::poly6::W(dist)
                                  *Kernel::poly6::W_coeff();
            }
            else
            {
                particle.velocity.x=0.;
                particle.Fother.x+=BoundaryConditions::bndCoeff*W0;
            }
        }
        else if (BoundaryConditions::bndBox.maxX()-iPos.x<scaledh)
        {
            if (iPos.x<BoundaryConditions::bndBox.maxX())
            {
                tmpiPos = glm::dvec2(iPos.x,0.0)*dScaledh;
                bndPos  = glm::dvec2(BoundaryConditions::bndBox.maxX(),0.0)*dScaledh;
                double dist = glm::length(tmpiPos-bndPos);
                particle.Fother.x-=BoundaryConditions::bndCoeff
                                  *Kernel::poly6::W(dist)
                                  *Kernel::poly6::W_coeff();
            }
            else
            {
                particle.velocity.x=0.;
                particle.Fother.x-=BoundaryConditions::bndCoeff*W0;
            }
        }
        if (iPos.y-BoundaryConditions::bndBox.minY()<scaledh)
        {
            if (iPos.y>BoundaryConditions::bndBox.minY())
            {
                tmpiPos = glm::dvec2(0.0,iPos.y)*dScaledh;
                bndPos  = glm::dvec2(0.0,BoundaryConditions::bndBox.minY())*dScaledh;
                double dist = glm::length(tmpiPos-bndPos);
                particle.Fother.y+=BoundaryConditions::bndCoeff
                                  *Kernel::poly6::W(dist)
                                  *Kernel::poly6::W_coeff();
            }
            else
            {
                particle.velocity.y=0.;
                particle.Fother.y+=BoundaryConditions::bndCoeff*W0;
            }
        }
        else if (BoundaryConditions::bndBox.maxY()-iPos.y<scaledh)
        {
            if (iPos.y>BoundaryConditions::bndBox.maxY())
            {
                tmpiPos = glm::dvec2(0.0,iPos.y)*dScaledh;
                bndPos  = glm::dvec2(0.0,BoundaryConditions::bndBox.maxY())*dScaledh;
                double dist = glm::length(tmpiPos-bndPos);
                particle.Fother.y-=BoundaryConditions::bndCoeff
                                  *Kernel::poly6::W(dist)
                                  *Kernel::poly6::W_coeff();
            }
            else
            {
                particle.velocity.y=0.;
                particle.Fother.y-=BoundaryConditions::bndCoeff*W0;
            }
        }

        particle.Fother*= particle.density;
    }
}

//********************************************************************************
void SPHSolver::calcPressForces()
//********************************************************************************
{
    static auto pressForceTimerID = Statistics::createTimer("SPHSolver::pressForceTimer");
    Statistics::TimerGuard pressForceGuard(pressForceTimerID);

    #pragma omp parallel for
    for (int iPart = 0;iPart<activeParticles_;iPart++) 
    {
        Particle& iParticle = cloud_[iPart];

        glm::dvec2 Fp = glm::dvec2(0.0);

        int Nnei = (int)iParticle.nei.size();
        for (int i = 0; i<Nnei;i++)
        {
            Neigbhor& iPartNeiI = iParticle.nei[i];
            int jPart = iPartNeiI.ID;
            if (iPart==jPart) continue;
            Particle& jParticle = cloud_[jPart];

            Fp+=jParticle.mass*jParticle.ddensity 
               *(iParticle.pressure + jParticle.pressure)
               *Kernel::spiky::gradW(iPartNeiI.dir,iPartNeiI.dist);
        }

        Fp*=-0.5*Kernel::spiky::gradW_coeff();
        
        iParticle.Fpress= Fp;
    }
}

//********************************************************************************
void SPHSolver::calcViscForces()
//********************************************************************************
{
    static auto viscForceTimerID  = Statistics::createTimer("SPHSolver::viscForceTimer");
    Statistics::TimerGuard viscForceGuard(viscForceTimerID);

    #pragma omp parallel for
    for (int iPart = 0;iPart<activeParticles_;iPart++) 
    {
        glm::dvec2 Fv = glm::dvec2(0.0);

        Particle& iParticle = cloud_[iPart];

        int Nnei = (int)iParticle.nei.size();
        //std::cout<<"iPart= "<<iPart<<" nNei="<<Nnei;
        for (int i = 0; i<Nnei;i++)
        {
            Neigbhor& iPartNeiI = iParticle.nei[i];
            int jPart = iPartNeiI.ID;
            //std::cout<<" "<<jPart;
            if (iPart==jPart) continue;
            Particle& jParticle = cloud_[jPart];

            Fv+= jParticle.mass
               *jParticle.ddensity
               *(jParticle.velocity-iParticle.velocity)
               *Kernel::visc::laplW(iPartNeiI.dist);
        }
        Fv*=SPHSettings::viscosity*Kernel::visc::laplW_coeff();

        iParticle.Fvisc= Fv;
    }
}

//********************************************************************************
void SPHSolver::calcSurfForces()
//********************************************************************************
{
    static auto surfForceTimerID  = Statistics::createTimer("SPHSolver::surfForceTimer");
    Statistics::TimerGuard surfForceGuard(surfForceTimerID);

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
            Neigbhor& iPartNeiI = iParticle.nei[i];
            int jPart = iPartNeiI.ID;
            //std::cout<<" "<<jPart;
            if (iPart==jPart) continue;
            Particle& jParticle = cloud_[jPart];

            correction = 2.*SPHSettings::particleDensity/(iParticle.density+jParticle.density);

            Fcohesion+= correction
                   *iParticle.mass*jParticle.mass 
                   *Kernel::surface::C(iPartNeiI.dist)
                   *iPartNeiI.dir/iPartNeiI.dist;

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
double SPHSolver::calcCFL() const
//********************************************************************************
{
    static auto clfTimerID  = Statistics::createTimer("SPHSolver::calcCFL");
    Statistics::TimerGuard clfGuard(clfTimerID);

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

//********************************************************************************
double SPHSolver::calcKineticEnergy() const
//********************************************************************************
{
    static auto kineticTimerID  = Statistics::createTimer("SPHSolver::calcKineticEnergy");
    Statistics::TimerGuard kinetickGuard(kineticTimerID);

    double kEnergy = 0.5*std::accumulate(cloud_.begin(), cloud_.end(), 0.,
        [](double k, const Particle& p) {
            return k + p.mass*glm::dot(p.velocity, p.velocity);
        }
    );

    return kEnergy;
}

