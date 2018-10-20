
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Description
    Main source file of SPHSimulator program. 
    2D SPH solver with WCSPH and PCISPH capabilities.

\************************************************************************/
#include <iostream>            

#include "Kernels.hpp"        
#include "SPHSolver.hpp"      
#include "WriteFunctions.hpp" 
#include "Statistics.hpp"     
#include "Settings.hpp"       

#define EnableGLRender 1

#if EnableGLRender
#include "DisplayView.hpp"    
#endif

void reportProgress(int iStep, SPHSolver& s)
{
    static auto wrtngTimerID = Statistics::createTimer("writingTime  ");

    if (iStep%SimulationSettings::showProgressEvery==0)
    {
        std::cout<<"Time="<<SimulationSettings::simTime<<" "
                 <<"Step="<< iStep<<" "
                 <<"CFL="<<s.calcCFL()<<" "
                 <<std::endl;
    }
    if (iStep%RenderSettings::printEvr==0)
    {
        Statistics::TimerGuard wrtTimer(wrtngTimerID);
        if (RenderSettings::fileRender==RenderSettings::RAWDATA)
        {
            writeRAWfile("step", s.cloud());
        }
        else
        {
            renderImage("render", s.cloud(), s.neibhs());
        }
    }
}

#if EnableGLRender
bool updateRender(DisplayView::WindowManager& w, SPHSolver& s, int iStep) {

    static auto rendeTimerID = Statistics::createTimer("renderingTime");
    static int  stepsPerSec = 0;
    ++stepsPerSec;

    auto windowRenderer = Statistics::CountTime(rendeTimerID, &DisplayView::WindowManager::renderParticles);

    bool ret = true;
    // Render to screen - 60 fps max
    if (w.shouldRender()) 
    {
        ret = windowRenderer(
            w, s.cloud(), s.NParticles()
        );
    }
    
    // Rename window - once per sec
    if (w.shouldRename()) 
    {
        std::string newName ="Time="+std::to_string(SimulationSettings::simTime)+" "
                            +"Step="+std::to_string(iStep)+" "
                            +"Steps/Sec="+std::to_string(stepsPerSec)+" "
                            +"CFL=" +std::to_string(s.calcCFL());
        w.renameWindow(newName);
        stepsPerSec=0;
    }
    return ret;
}
#endif


int main()
{
    //init kernel - smoothing lengths
    Kernel::SmoothingLength::setSmoothingLength(SPHSettings::initDx*2.);

    SPHSolver SPHsolver;
    SPHsolver.init();


    auto totalTimerID = Statistics::createTimer("totalTime    ");
    auto compuTimerID = Statistics::createTimer("SPHSolver    ");

    auto SPHSolverStep  = Statistics::CountTime(compuTimerID, &SPHSolver::step);

#if EnableGLRender
    DisplayView::WindowManager windowManager;
    windowManager.init(SPHsolver.cloud());

#endif

    bool continueLoop = true;
    int  iStep = 0;

    std::cout<<"Starting Simulation"<<std::endl;
    {
        Statistics::TimerGuard g(totalTimerID);
        while (continueLoop)
        {
            iStep++;
        
            continueLoop &= SPHSolverStep(SPHsolver);
        
            reportProgress(iStep, SPHsolver);
        
#if EnableGLRender
            continueLoop &= updateRender(windowManager, SPHsolver, iStep);
#endif
        }
    }

    Statistics::printStatistics();
}

