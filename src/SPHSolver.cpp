
#include "SPHSolver.hpp"
#include "WriteFunctions.hpp"

#include <chrono>

void SPHSolver::init()
{
    neibhs_.setHashTable(BoundaryConditions::bndBox.minPos(),BoundaryConditions::bndBox.delta());
    //neibhs_.writeGridGNU("neiGrid");

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
        for (int i=0;i<SPHSettings::NParticles;i++)
        {
            cloud_[i].active = false;
        }
    }

    neibhs_.findNei(cloud_, SPHSettings::NParticles);
    calcDensity();

    if (RenderSettings::fileRender==RenderSettings::GNUPLOT)
    {
        writeGNUfile("step", cloud_);
    }
    else
    {
        renderImage("render", cloud_, neibhs_);
    }
}

void SPHSolver::WCSPHStep() 
{
    calcDensity();
    calcPressure();
    
    //calc forces
    calcPressForces();
    calcViscForces();
    calcOtherForces();

    //combine forces and update values
    #pragma omp parallel for
    for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
    {
        if (!cloud_[iPart].active) continue;
        cloud_[iPart].Ftot = cloud_[iPart].Fpress+cloud_[iPart].Fvisc+cloud_[iPart].Fother;

        cloud_[iPart].velocity += SimulationSettings::dt*cloud_[iPart].ddensity*cloud_[iPart].Ftot;
        cloud_[iPart].position += SimulationSettings::dt*cloud_[iPart].velocity;
    }
}

void SPHSolver::PCISPHStep() 
{
    calcDensity();
    //calcPressure();
    
    //calc forces
    calcViscForces();
    calcOtherForces();
    //calcPressForces();

    #pragma omp parallel for
    for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
    {
        if (!cloud_[iPart].active) continue;
        cloud_[iPart].Ftot = cloud_[iPart].Fvisc+cloud_[iPart].Fother;

        cloud_[iPart].velocity += SimulationSettings::dt*cloud_[iPart].ddensity*cloud_[iPart].Ftot;
        cloud_[iPart].position += SimulationSettings::dt*cloud_[iPart].velocity;
    }

    double densErr=1000001.;
    int iter=0;

    initPressure();
    const double targetDensErr = SPHSettings::densityErr*SPHSettings::particleDensity;
    while(densErr>targetDensErr && iter<500)
    {
        calcDensity();
        calcDensityErr();
        updatePressure();

        densErr = 0.;
        for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
        {
            if (!cloud_[iPart].active) continue;
            if (densErr<cloud_[iPart].densityErr)
                densErr = cloud_[iPart].densityErr;
        }

        if (iter==0)
        //std::cout<<"densErrInit="<<densErr;
        iter++;

        calcPressForces();
          //glm::dvec2 Ftot = Fpress[iPart]+Fvisc[iPart]+Fother[iPart];
          //cloud_[iPart].position += SimulationSettings::dt*cloud_[iPart].velocity;
        //#pragma omp parallel for
        for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
        {
            if (!cloud_[iPart].active) continue;
            cloud_[iPart].Ftot = cloud_[iPart].Fpress;

            glm::dvec2 update = SimulationSettings::dt*cloud_[iPart].ddensity*cloud_[iPart].Ftot;
            cloud_[iPart].velocity += update;
            cloud_[iPart].position += SimulationSettings::dt*update;
        }
    }
    //std::cout<<" densErrFinal="<<densErr<<" iter="<<iter<<std::endl;
}

void SPHSolver::step(GLFWwindow* window) 
{
    generateParticles();
    neibhs_.findNei(cloud_, SPHSettings::NParticles);
    switch (SPHSettings::SPHstep)
    {
        case SPHSettings::Solver::WCSPH:
        {
            WCSPHStep();
            break;
        }
        case SPHSettings::Solver::PCISPH:
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

    double vMax = 0.0;
    double tMax = 0.0;

    for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
    {
        if (!cloud_[iPart].active) continue;
        double vLen2 = glm::length2(cloud_[iPart].velocity);
        if (vLen2 > vMax) vMax = vLen2;
        tMax+=vLen2;
    }

    
    double maxDt = 0.4*Kernel::SmoothingLength::h/(sqrt(vMax)+1.e-5);
    double CFL = SimulationSettings::dt/maxDt;

    static int iStep=0;
    static int printStep=0;
    static std::chrono::time_point<std::chrono::system_clock> codeStart = std::chrono::system_clock::now();
    static std::chrono::time_point<std::chrono::system_clock> start     = std::chrono::system_clock::now();
    static std::chrono::time_point<std::chrono::system_clock> now       = std::chrono::system_clock::now();
    static std::chrono::duration<double> deltaTime =std::chrono::duration<double>(0.0);

    SimulationSettings::updateSimTime();
    iStep++;
    printStep++;
    now = std::chrono::system_clock::now();
    deltaTime = now - start;

    if (SimulationSettings::breakLoop()) 
    {
	    if (!glfwWindowShouldClose(window))
        {
            now = std::chrono::system_clock::now();
            deltaTime = now - codeStart;
            double dTime =  deltaTime.count();
            std::cout<<"SimulationTime="<<SimulationSettings::simTime<<"\n"
                     <<"compTime      ="<<dTime<<"\n"
                     <<"nSteps        ="<<printStep<<"\n"
                     <<"dt            ="<<SimulationSettings::dt<<"\n"
                     <<"timePerStep   ="<<dTime/printStep<<"\n"
                     <<"simTimePerSec ="<<SimulationSettings::simTime/dTime<<std::endl;
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }
    
    if (deltaTime.count()>1.)
    {
        std::string winName = "time="+std::to_string(SimulationSettings::simTime)+" CFL="+std::to_string(CFL)+" nStep="+std::to_string(iStep)+" tVel="+std::to_string(tMax);
        //glutSetWindowTitle(winName.c_str());
        glfwSetWindowTitle(window, winName.c_str());
        start = std::chrono::system_clock::now();
        iStep = 0;
    }

    if (printStep%RenderSettings::printEvr==0)
    {
        if (RenderSettings::fileRender==RenderSettings::GNUPLOT)
        {
            writeGNUfile("step", cloud_);
        }
        else
        {
            renderImage("render", cloud_, neibhs_);
        }
    }
}

void SPHSolver::generateParticles()
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
                        cloud_[activeParticles_].active = true;
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
                            cloud_[activeParticles_].active = true;
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

void SPHSolver::calcDensity()
{
  //for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
#pragma omp parallel for
    for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
    {
        if (!cloud_[iPart].active) continue;
        double dens = 0.0;

        //double tmpValue = cloud_[iPart].pressure/(cloud_[iPart].density*cloud_[iPart].density);
        int Nnei = (int)cloud_[iPart].nei.size();
        for (int i = 0; i<Nnei;i++)
        {
            int neiPart = cloud_[iPart].nei[i];
            Particle& neiParticle = cloud_[neiPart];

            double Wij = Kernel::poly6::W(cloud_[iPart].position,neiParticle.position);
            dens += neiParticle.mass*Wij;
        }
        cloud_[iPart].density = dens;
        cloud_[iPart].ddensity = 1./dens;
    }
}

void SPHSolver::calcDensityErr()
{
  //for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
#pragma omp parallel for
    for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
    {
        if (!cloud_[iPart].active) continue;
        cloud_[iPart].densityErr = cloud_[iPart].density - SPHSettings::particleDensity;
    }
}

void SPHSolver::initPressure()
{
  //for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
#pragma omp parallel for
    for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
    {
        if (!cloud_[iPart].active) continue;
        cloud_[iPart].pressure = 0.0;
    }
}

void SPHSolver::calcPressure()
{
#pragma omp parallel for
    for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
    {
        if (!cloud_[iPart].active) continue;
        //p = k(rho-rho0)
        cloud_[iPart].pressure = fmax(SPHSettings::stiffness*(cloud_[iPart].density-SPHSettings::particleDensity),0.0);
      //cloud_[iPart].pressure = SPHSettings::stiffness*(cloud_[iPart].density-SPHSettings::particleDensity);
        //p = K((rho/rho0)-1)
      //double tmp = cloud_[iPart].density*SPHSettings::dParticleDensity;
      //double tmp2 = tmp*tmp*tmp*tmp*tmp*tmp*tmp;
      //cloud_[iPart].pressure = fmax(1./7.*SPHSettings::stiffness*SPHSettings::particleDensity*(tmp2-1.),0.);
      //cloud_[iPart].pressure = 1./7.*SPHSettings::stiffness*SPHSettings::particleDensity*(tmp2-1.);
    }
}

void SPHSolver::updatePressure()
{
#pragma omp parallel for
    for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
    {
        if (!cloud_[iPart].active) continue;
        cloud_[iPart].pressure = 1.*fmax(delta_*cloud_[iPart].densityErr,0.0);
    }
}


void SPHSolver::calcOtherForces()
{
//#pragma omp parallel for
    for (int iPart=0;iPart<SPHSettings::NParticles;iPart++)
    {
        if (!cloud_[iPart].active) continue;
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
      //cloud_[iPart].Fother+=BoundaryConditions::bndCoeff*Kernel::poly6::W(tmpiPos,bndPos)*glm::normalize(tmpiPos-bndPos);

        cloud_[iPart].Fother*= cloud_[iPart].density;
    }
}

void SPHSolver::calcPressForces()
{
    for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
    {
        if (!cloud_[iPart].active) continue;
        glm::dvec2 Fp = glm::dvec2(0.0);

        Particle& iParticle = cloud_[iPart];
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

void SPHSolver::calcViscForces()
{
    for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
    {
        glm::dvec2 Fv = glm::dvec2(0.0);

        Particle& iParticle = cloud_[iPart];
        //double tmpValue = cloud_[iPart].pressure/(cloud_[iPart].density*cloud_[iPart].density);
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
