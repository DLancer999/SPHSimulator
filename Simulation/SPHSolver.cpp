
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
#include "Reorderer.hpp"
#include "Statistics/Statistics.hpp"

#include "Util/ParallelFor.hpp"

#include <boost/range/combine.hpp>
#include <boost/range/algorithm/for_each.hpp>

//********************************************************************************
void SPHSolver::init()
//********************************************************************************
{
  neibhs_.setHashTable(BoundaryConditions::bndBox.minPos(),BoundaryConditions::bndBox.delta());
  //neibhs_.writeGridRAW("neiGrid");

  cloud_.reserve(SPHSettings::NParticles);

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
    unsigned nPart = 0;
    unsigned i = 0;
    while (nPart<SPHSettings::NParticles)
    {
      for (unsigned j=0;j<SPHSettings::LParticles;j++)
      {
        if (cloud_.size()>=SPHSettings::NParticles) break;
        Particle p;
        p.get<Attr::ePosition>() = glm::dvec2
          (
           InitialConditions::particleInitPos.x+(SPHSettings::initDx*(double(j)+0.5)),
           InitialConditions::particleInitPos.y+SPHSettings::initDx*(double(i)+0.5)
          ); 
        p.get<Attr::eVelocity>() = InitialConditions::particleInitVel; 
        p.get<Attr::eDensErr>()  = SPHSettings::particleDensity;
        p.get<Attr::eDDensity>() = SPHSettings::dParticleDensity;
        p.get<Attr::eMass>()     = SPHSettings::particleMass;
        nPart++;
        cloud_.push_back(p);
      }
      i++;
    }
  }
  else
  {
    //particles will be generated in generateParticles()
  }

  neibhs_.clear();
  neibhs_.findNei(cloud_);
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

    auto& particlePos = cloud_.get<Attr::ePosition>();
    auto& particleVel = cloud_.get<Attr::eVelocity>();
    const auto& particleFPress = cloud_.get<Attr::ePressForce>();
    const auto& particleFVisc  = cloud_.get<Attr::eViscForce >();
    const auto& particleFSurf  = cloud_.get<Attr::eSurfForce >();
    const auto& particleFOther = cloud_.get<Attr::eOtherForce>();
    const auto& particleDDens  = cloud_.get<Attr::eDDensity  >();

    auto& particleFTotal = cloud_.get<Attr::eTotalForce>();

    Parallel::For (cloud_.size(), [ &particlePos,   &particleVel,    &particleFPress, &particleFVisc,
                                    &particleFSurf, &particleFOther, &particleDDens,  &particleFTotal
                                  ] (size_t iPart)
    {
      particleFTotal[iPart] = particleFPress[iPart] + particleFVisc [iPart]
                            + particleFSurf [iPart] + particleFOther[iPart];

      particleVel[iPart] += SimulationSettings::dt*particleDDens[iPart]*particleFTotal[iPart];
      particlePos[iPart] += SimulationSettings::dt*particleVel[iPart];
    });
  }
}

//********************************************************************************
void SPHSolver::PCISPHStep() 
//********************************************************************************
{
  auto& particlePos = cloud_.get<Attr::ePosition>();
  auto& particleVel = cloud_.get<Attr::eVelocity>();
  auto& particleFTot = cloud_.get<Attr::eTotalForce>();
  const auto& particleFPress  = cloud_.get<Attr::ePressForce>();
  const auto& particleDDens   = cloud_.get<Attr::eDDensity>();
  const auto& particleDensErr = cloud_.get<Attr::eDensErr>();

  calcDensity();
  calcNormal();
  
  //calc forces
  calcViscForces();
  calcSurfForces();
  calcOtherForces();

  const size_t nPart = cloud_.size();
  {
    static auto updatePosTimerID  = Statistics::createTimer("SPHSolver::PCISPH::updatePosTimer");
    Statistics::TimerGuard g(updatePosTimerID);

    const auto& particleFVisc  = cloud_.get<Attr::eViscForce >();
    const auto& particleFSurf  = cloud_.get<Attr::eSurfForce >();
    const auto& particleFOther = cloud_.get<Attr::eOtherForce>();

    auto& particleFTotal = cloud_.get<Attr::eTotalForce>();

    Parallel::For (nPart, [ &particleFTotal, &particleVel,    &particlePos, &particleFVisc,
                            &particleFSurf,  &particleFOther, &particleDDens
                          ](size_t iPart){
      particleFTotal[iPart] = particleFVisc [iPart]
                            + particleFSurf [iPart]
                            + particleFOther[iPart];

      particleVel[iPart] += SimulationSettings::dt*particleDDens[iPart]*particleFTotal[iPart];
      particlePos[iPart] += SimulationSettings::dt*particleVel[iPart];
    });
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

      auto itMax = std::max_element(particleDensErr.begin(), particleDensErr.end());
      densErr = (itMax != particleDensErr.end()) ? (*itMax) : 0.;

      if (iter==0)
      iter++;

      calcPressForces();

      const double dt = SimulationSettings::dt;

      Parallel::For (nPart, [ &particleFTot, &particleFPress, &particleVel,
                              &particlePos,  &particleDDens, dt
                            ] (size_t iPart){
        const glm::dvec2 iPress = particleFPress[iPart];
        particleFTot[iPart] += iPress;

        glm::dvec2 update = dt*particleDDens[iPart]*iPress;
        particleVel[iPart] += update;
        particlePos[iPart] += dt*update;
      });
    }
  }
}

//********************************************************************************
bool SPHSolver::step() 
//********************************************************************************
{
  static int iReorder = 0;
  iReorder++;

  auto& particleNei = cloud_.get<Attr::eNei>();
  boost::for_each(particleNei, [](auto& v){ v.clear();});

  if(iReorder == 100)
  {
    Reorderer::reorderCloud(cloud_);
    iReorder = 0;
  }

  generateParticles();

  neibhs_.clear();

  neibhs_.findNei(cloud_);

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

  if (cloud_.size()>=SPHSettings::NParticles) return;
   
  switch(InitialConditions::particleGeneration)
  {
    case (InitialConditions::FAUCET):
    {
      if (SimulationSettings::simTime-genTime>=dtGenFaucet)
      {
        glm::dvec2 initVelNormalized = glm::normalize(InitialConditions::particleInitVel);
        glm::dvec2 generationDirection = glm::dvec2(-initVelNormalized.y,initVelNormalized.x);
        //std::cout<<"generating"<<std::endl;
        if (cloud_.size()<SPHSettings::NParticles)
        {
          for (unsigned i=0;i<SPHSettings::LParticles;i++)
          {
            if (cloud_.size()>=SPHSettings::NParticles) break;
            Particle active;
            active.get<Attr::ePosition>() = InitialConditions::particleInitPos + i*SPHSettings::initDx*generationDirection;
            active.get<Attr::eVelocity>() = InitialConditions::particleInitVel;
            active.get<Attr::eMass>()     = SPHSettings::particleMass;
            cloud_.push_back(active);
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
        if (cloud_.size()<SPHSettings::NParticles)
        {
          //double     radius = 0.3*double(SPHSettings::LParticles)*SPHSettings::initDx; 
          double     radius = 0.4*double(SPHSettings::LParticles)*SPHSettings::initDx; 
          double     offset = 0.5*SPHSettings::initDx*double(SPHSettings::LParticles-1);
          glm::dvec2 center = InitialConditions::particleInitPos + glm::dvec2(offset,offset);

          for (unsigned i=0;i<SPHSettings::LParticles;i++)
          {
            for (unsigned j=0;j<SPHSettings::LParticles;j++)
            {
              if (cloud_.size()>=SPHSettings::NParticles) break;
              glm::dvec2 pos(InitialConditions::particleInitPos + glm::dvec2(SPHSettings::initDx*i,SPHSettings::initDx*j));
              double dist = glm::length(pos-center);
              if (dist>radius) continue;
              Particle active;
              active.get<Attr::ePosition>() = pos;
              active.get<Attr::eVelocity>() = InitialConditions::particleInitVel;
              active.get<Attr::eMass    >() = SPHSettings::particleMass;
              cloud_.push_back(active);
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

  const auto& particlePos = cloud_.get<Attr::ePosition>();
  auto& particleNei = cloud_.get<Attr::eNei>();

  Parallel::For (cloud_.size(), [ &particlePos, &particleNei ] (size_t iPart) {
    const glm::dvec2 iPos = particlePos[iPart];
    for ( Neigbhor& iPartNeiI :  particleNei[iPart] )
    {
        const glm::dvec2 dir = iPos-particlePos[iPartNeiI.ID];
        iPartNeiI.dir  = dir;
        iPartNeiI.dist = glm::length(dir);
    }
  });
}

//********************************************************************************
void SPHSolver::calcDensity()
//********************************************************************************
{
  static auto densCalcTimerID   = Statistics::createTimer("SPHSolver::densCalcTimer");
  Statistics::TimerGuard densGuard(densCalcTimerID);

  const auto& particleMass = cloud_.get<Attr::eMass>();
  const auto& particleNei  = cloud_.get<Attr::eNei>();

  auto& particleDens  = cloud_.get<Attr::eDensity>();
  auto& particleDDens = cloud_.get<Attr::eDDensity>();

  const double Wcoeff = Kernel::poly6::W_coeff();

  Parallel::For (cloud_.size(), [ &particleMass, &particleNei,
                                  &particleDens, &particleDDens, Wcoeff
                                ] (size_t iPart) {
    const auto& iNei = particleNei[iPart];

    double dens = 0.0;

    const unsigned Nnei = unsigned(iNei.size());
    for (unsigned i = 0; i<Nnei;i++)
    {
        unsigned jPart = iNei[i].ID;

        dens += particleMass[jPart]*Kernel::poly6::W(iNei[i].dist);
    }
    dens*=Wcoeff;

    particleDens[iPart]  = dens;
    particleDDens[iPart] = 1./dens;
  });
}

//********************************************************************************
void SPHSolver::calcDensityErr()
//********************************************************************************
{
  static auto densErrCalcTimerID   = Statistics::createTimer("SPHSolver::densErrCalcTimer");
  Statistics::TimerGuard densErrGuard(densErrCalcTimerID);

  const auto& particleDens = cloud_.get<Attr::eDensity>();
  auto& particleDensErr = cloud_.get<Attr::eDensErr>();

  const double pDens = SPHSettings::particleDensity;

  Parallel::For (cloud_.size(), [ &particleDensErr, &particleDens, pDens ] (size_t iPart) {
    particleDensErr[iPart] = particleDens[iPart] - pDens;
  });
}

//********************************************************************************
void SPHSolver::initPressure()
//********************************************************************************
{
  auto& particlePress = cloud_.get<Attr::ePressure>();
  Parallel::For (cloud_.size(), [ &particlePress ] (size_t iPart) {
    particlePress[iPart]=0.0;
  });
}

//********************************************************************************
void SPHSolver::calcPressure()
//********************************************************************************
{
  static auto pressCalcTimerID   = Statistics::createTimer("SPHSolver::pressureCalcTimer");
  Statistics::TimerGuard pressGuard(pressCalcTimerID);

  const auto& particleDens = cloud_.get<Attr::eDensity>();
  auto& particlePress = cloud_.get<Attr::ePressure>();

  const double stiff = SPHSettings::stiffness;
  const double pDens = SPHSettings::particleDensity;

  Parallel::For (cloud_.size(), [ &particlePress, &particleDens, stiff, pDens ] (size_t iPart) {
    //p = k(rho-rho0)
    particlePress[iPart] = fmax(stiff*(particleDens[iPart] - pDens), 0.0);
  });
}

//********************************************************************************
void SPHSolver::calcNormal()
//********************************************************************************
{
  static auto normalCalcTimerID   = Statistics::createTimer("SPHSolver::normalCalcTimer");
  Statistics::TimerGuard normalGuard(normalCalcTimerID);

  const auto& particleMass  = cloud_.get<Attr::eMass>();
  const auto& particleDDens = cloud_.get<Attr::eDDensity>();
  const auto& particleNei  = cloud_.get<Attr::eNei>();
  auto& particleNormal = cloud_.get<Attr::eNormal>();

  Parallel::For (cloud_.size(), [ &particleMass, &particleDDens, &particleNei, &particleNormal] (size_t iPart) {
    const auto& iNei = particleNei[iPart];

    glm::dvec2 norm (0.0);

    unsigned Nnei = unsigned(iNei.size());
    for (unsigned i = 0; i<Nnei;i++)
    {
        unsigned jPart = iNei[i].ID;

        norm += particleMass[jPart]*particleDDens[jPart]
              *Kernel::poly6::gradW(iNei[i].dir,iNei[i].dist);
    }
    norm *= Kernel::poly6::gradW_coeff()*Kernel::SmoothingLength::h;
    particleNormal[iPart] = norm;
  });
}


//********************************************************************************
void SPHSolver::updatePressure()
//********************************************************************************
{
  static auto updatePressTimerID   = Statistics::createTimer("SPHSolver::updatePressure");
  Statistics::TimerGuard updatePressGuard(updatePressTimerID);

  const auto& particleDensErr = cloud_.get<Attr::eDensErr>();
  auto& particlePress = cloud_.get<Attr::ePressure>();

  const double delta = delta_;

  Parallel::For (cloud_.size(), [ &particlePress, &particleDensErr, delta] (size_t iPart) {
    particlePress[iPart] = fmax(delta*particleDensErr[iPart],0.0);
    //particlePress[iPart] = delta*particleDensErr[iPart];
  });
}


//********************************************************************************
void SPHSolver::calcOtherForces()
//********************************************************************************
{
  static auto otherForceTimerID = Statistics::createTimer("SPHSolver::otherForceTimer");
  Statistics::TimerGuard otherForceGuard(otherForceTimerID);

  const auto& particlePos = cloud_.get<Attr::ePosition>();
  const auto& particleDens = cloud_.get<Attr::eDensity>();
  auto& particleVel = cloud_.get<Attr::eVelocity>();
  auto& particleFOther = cloud_.get<Attr::eOtherForce>();

  const glm::dvec2 grav = SPHSettings::grav;

  Parallel::For (cloud_.size(), [ &particlePos, &particleDens, &particleVel, &particleFOther, grav] (size_t iPart) {
    //gravity
    particleFOther[iPart] = grav;

    //boundary forces
    glm::dvec2 iPos = particlePos[iPart];
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
        particleFOther[iPart].x+=BoundaryConditions::bndCoeff
                                *Kernel::poly6::W(dist)
                                *Kernel::poly6::W_coeff();
      }
      else
      {
        particleVel[iPart].x=0.;
        particleFOther[iPart].x+=BoundaryConditions::bndCoeff*W0;
      }
    }
    else if (BoundaryConditions::bndBox.maxX()-iPos.x<scaledh)
    {
      if (iPos.x<BoundaryConditions::bndBox.maxX())
      {
        tmpiPos = glm::dvec2(iPos.x,0.0)*dScaledh;
        bndPos  = glm::dvec2(BoundaryConditions::bndBox.maxX(),0.0)*dScaledh;
        double dist = glm::length(tmpiPos-bndPos);
        particleFOther[iPart].x-=BoundaryConditions::bndCoeff
                                *Kernel::poly6::W(dist)
                                *Kernel::poly6::W_coeff();
      }
      else
      {
        particleVel[iPart].x=0.;
        particleFOther[iPart].x-=BoundaryConditions::bndCoeff*W0;
      }
    }
    if (iPos.y-BoundaryConditions::bndBox.minY()<scaledh)
    {
      if (iPos.y>BoundaryConditions::bndBox.minY())
      {
        tmpiPos = glm::dvec2(0.0,iPos.y)*dScaledh;
        bndPos  = glm::dvec2(0.0,BoundaryConditions::bndBox.minY())*dScaledh;
        double dist = glm::length(tmpiPos-bndPos);
        particleFOther[iPart].y+=BoundaryConditions::bndCoeff
                               *Kernel::poly6::W(dist)
                               *Kernel::poly6::W_coeff();
      }
      else
      {
        particleVel[iPart].y=0.;
        particleFOther[iPart].y+=BoundaryConditions::bndCoeff*W0;
      }
    }
    else if (BoundaryConditions::bndBox.maxY()-iPos.y<scaledh)
    {
      if (iPos.y>BoundaryConditions::bndBox.maxY())
      {
        tmpiPos = glm::dvec2(0.0,iPos.y)*dScaledh;
        bndPos  = glm::dvec2(0.0,BoundaryConditions::bndBox.maxY())*dScaledh;
        double dist = glm::length(tmpiPos-bndPos);
        particleFOther[iPart].y-=BoundaryConditions::bndCoeff
                                *Kernel::poly6::W(dist)
                                *Kernel::poly6::W_coeff();
      }
      else
      {
        particleVel[iPart].y=0.;
        particleFOther[iPart].y-=BoundaryConditions::bndCoeff*W0;
      }
    }

    particleFOther[iPart] *= particleDens[iPart];
  });
}

//********************************************************************************
void SPHSolver::calcPressForces()
//********************************************************************************
{
  static auto pressForceTimerID = Statistics::createTimer("SPHSolver::pressForceTimer");
  Statistics::TimerGuard pressForceGuard(pressForceTimerID);

  const auto& particleMass  = cloud_.get<Attr::eMass>();
  const auto& particleDDens = cloud_.get<Attr::eDDensity>();
  const auto& particlePress = cloud_.get<Attr::ePressure>();
  const auto& particleNei   = cloud_.get<Attr::eNei>();
  auto& particleFPress = cloud_.get<Attr::ePressForce>();

  Parallel::For (cloud_.size(), [ &particleMass,  &particleDDens, &particlePress,
                                  &particleNei, &particleFPress
                                ] (size_t iPart) {
    const auto& iNei = particleNei[iPart];

    glm::dvec2 Fp = glm::dvec2(0.0);
    const double iPress = particlePress[iPart];

    const unsigned Nnei = unsigned(iNei.size());
    for (unsigned i = 0; i<Nnei;i++)
    {
      const Neigbhor& iPartNeiI = iNei[i];
      unsigned jPart = iPartNeiI.ID;
      if (iPart==jPart) continue;

      Fp+=particleMass[jPart]*particleDDens[jPart]
         *(iPress + particlePress[jPart])
         *Kernel::spiky::gradW(iPartNeiI.dir,iPartNeiI.dist);
    }

    Fp*=-0.5*Kernel::spiky::gradW_coeff();
    
    particleFPress[iPart] = Fp;
  });
}

//********************************************************************************
void SPHSolver::calcViscForces()
//********************************************************************************
{
  static auto viscForceTimerID  = Statistics::createTimer("SPHSolver::viscForceTimer");
  Statistics::TimerGuard viscForceGuard(viscForceTimerID);
  
  const auto& particleVel   = cloud_.get<Attr::eVelocity>();
  const auto& particleMass  = cloud_.get<Attr::eMass>();
  const auto& particleDDens = cloud_.get<Attr::eDDensity>();
  const auto& particleNei   = cloud_.get<Attr::eNei>();
  auto& particleFVisc = cloud_.get<Attr::eViscForce>();
  
  Parallel::For (cloud_.size(), [ &particleVel, &particleMass, &particleDDens,
                                  &particleNei, &particleFVisc
                                ] (size_t iPart) {
    glm::dvec2 Fv = glm::dvec2(0.0);
  
    glm::dvec2 iVel = particleVel[iPart];
    const auto& iNei = particleNei[iPart];
  
    const unsigned Nnei = unsigned(iNei.size());
    //std::cout<<"iPart= "<<iPart<<" nNei="<<Nnei;
    for (unsigned i = 0; i<Nnei;i++)
    {
      const Neigbhor& iPartNeiI = iNei[i];
      unsigned jPart = iPartNeiI.ID;
      if (iPart==jPart) continue;
  
      Fv+= particleMass[jPart]
         *particleDDens[jPart]
         *(particleVel[jPart]-iVel)
         *Kernel::visc::laplW(iPartNeiI.dist);
    }
    Fv*=SPHSettings::viscosity*Kernel::visc::laplW_coeff();
  
    particleFVisc[iPart] = Fv;
  });
}

//********************************************************************************
void SPHSolver::calcSurfForces()
//********************************************************************************
{
  static auto surfForceTimerID  = Statistics::createTimer("SPHSolver::surfForceTimer");
  Statistics::TimerGuard surfForceGuard(surfForceTimerID);

  const auto& particleNormal = cloud_.get<Attr::eNormal>();
  const auto& particleMass = cloud_.get<Attr::eMass>();
  const auto& particleDens = cloud_.get<Attr::eDensity>();
  const auto& particleNei  = cloud_.get<Attr::eNei>();
  auto& particleFSurf = cloud_.get<Attr::eSurfForce>();

  Parallel::For (cloud_.size(), [ &particleNormal, &particleMass, &particleDens,
                                  &particleNei, &particleFSurf
                                ] (size_t iPart) {
    glm::dvec2 Fcohesion  = glm::dvec2(0.0);
    glm::dvec2 Fcurvature = glm::dvec2(0.0);
    double correction = 0.0;

    const glm::dvec2 iNorm = particleNormal[iPart];
    const double iMass = particleMass[iPart];
    const double iDens = particleDens[iPart];
    const auto&  iNei  = particleNei[iPart];

    const unsigned Nnei = unsigned(iNei.size());
    for (unsigned i = 0; i<Nnei;i++)
    {
      const Neigbhor& iPartNeiI = iNei[i];
      const unsigned jPart = iPartNeiI.ID;
      //std::cout<<" "<<jPart;
      if (iPart==jPart) continue;

      correction = 2.*SPHSettings::particleDensity/(iDens + particleDens[jPart]);

      Fcohesion+= correction
             *iMass*particleMass[jPart]
             *Kernel::surface::C(iPartNeiI.dist)
             *iPartNeiI.dir/iPartNeiI.dist;

      Fcurvature+= correction
             *iMass
             *(iNorm - particleNormal[jPart]);
    }

    Fcohesion *= Kernel::surface::C_coeff();

    particleFSurf[iPart] = -SPHSettings::surfTension
                         * iDens
                         * (Fcohesion+Fcurvature);
  });
}

//********************************************************************************
double SPHSolver::calcCFL() const
//********************************************************************************
{
  static auto clfTimerID  = Statistics::createTimer("SPHSolver::calcCFL");
  Statistics::TimerGuard clfGuard(clfTimerID);

  const auto& particleVel = cloud_.get<Attr::eVelocity>();

  double vMax = 0.0;

  for (size_t iPart = 0, nPart = cloud_.size(); iPart<nPart; ++iPart) 
  {
    double vLen2 = glm::length2(particleVel[iPart]);
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

  const auto& particleVel  = cloud_.get<Attr::eVelocity>();
  const auto& particleMass = cloud_.get<Attr::eMass>();
  auto accumRange = boost::combine(particleMass, particleVel);

  double kEnergy = 0.5*std::accumulate(accumRange.begin(), accumRange.end(), 0.,
    [](double k, const auto& elem) {
      const double mass = boost::get<0>(elem);
      const glm::dvec2& vel = boost::get<1>(elem);
      return k + mass*glm::dot(vel, vel);
    }
  );

  return kEnergy;
}

