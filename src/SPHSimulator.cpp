
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
#include "DisplayView.hpp"    
#include "Statistics.hpp"     
#include "Settings.hpp"       

int main()
{
    //init kernel - smoothing lengths
    Kernel::SmoothingLength::setSmoothingLength(SPHSettings::initDx*2.);

    SPHSolver SPHsolver;
    SPHsolver.init();

    DisplayView::WindowManager windowManager;
    windowManager.init(SPHsolver.cloud());

    int totalTimerID = Statistics::createTimer("totalTime    ");
    int compuTimerID = Statistics::createTimer("SPHSolver    ");
    int wrtngTimerID = Statistics::createTimer("writingTime  ");
    int rendeTimerID = Statistics::createTimer("renderingTime");

    bool continueLoop = true;
    int  iStep = 0;
    int  stepsPerSec = 0;

    std::cout<<"Starting Simulation"<<std::endl;
	while (continueLoop)
	{
        iStep++;
        stepsPerSec++;

        //solve
        Statistics::timers[compuTimerID].start();
        continueLoop &= SPHsolver.step();
        Statistics::timers[compuTimerID].end();
        Statistics::timers[compuTimerID].addTime();

        //print progresss
        if (iStep%SimulationSettings::showProgressEvery==0)
        {
            std::cout<<"Time="<<SimulationSettings::simTime<<" "
                     <<"Step="<< iStep<<" "
                     <<"CFL="<<SPHsolver.calcCFL()<<" "
                     <<std::endl;
        }

        //write to file
        if (iStep%RenderSettings::printEvr==0)
        {
            Statistics::timers[wrtngTimerID].start();
            if (RenderSettings::fileRender==RenderSettings::RAWDATA)
            {
                writeRAWfile("step", SPHsolver.cloud());
            }
            else
            {
                renderImage("render", SPHsolver.cloud(), SPHsolver.neibhs());
            }
            Statistics::timers[wrtngTimerID].end();
            Statistics::timers[wrtngTimerID].addTime();
        }

	    // Render to screen - 60 fps max
        if (windowManager.shouldRender()) 
        {
            Statistics::timers[rendeTimerID].start();
            continueLoop &= windowManager.renderParticles(SPHsolver.cloud(), SPHsolver.NParticles());
            Statistics::timers[rendeTimerID].end();
            Statistics::timers[rendeTimerID].addTime();
        }

	    // rename window - once per sec
        if (windowManager.shouldRename()) 
        {
            std::string newName ="Time="+std::to_string(SimulationSettings::simTime)+" "
                                +"Step="+std::to_string(iStep)+" "
                                +"Steps/Sec="+std::to_string(stepsPerSec)+" "
                                +"CFL=" +std::to_string(SPHsolver.calcCFL());
            windowManager.renameWindow(newName);
            stepsPerSec=0;
        }
	}
    Statistics::timers[totalTimerID].end();
    Statistics::timers[totalTimerID].addTime();

    Statistics::printStatistics();
}

