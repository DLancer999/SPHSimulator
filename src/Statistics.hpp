
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Statistics
 
Description
    Class to measure performance/timers

SourceFiles
    Statistics.cpp

\************************************************************************/

#ifndef STATISTICSS_H
#define STATISTICSS_H

#include <chrono>
#include <iostream>
#include <vector>

#define COUNT_TIME(FUNCTION,timerID)\
Statistics::timers[timerID].start();\
FUNCTION;\
Statistics::timers[timerID].end();\
Statistics::timers[timerID].addTime();

typedef std::chrono::time_point<std::chrono::system_clock> Clock;
typedef std::chrono::duration<double> Duration;

class Timer
{
private:
    Clock start_;
    Clock end_;
    Duration duration_;
    std::string timerName_;
public:
    Timer():
    start_(std::chrono::system_clock::now()),
    end_(std::chrono::system_clock::now()),
    duration_(end_-start_),
    timerName_("")
    {}

    Timer(const std::string& name):
    start_(std::chrono::system_clock::now()),
    end_(std::chrono::system_clock::now()),
    duration_(end_-start_),
    timerName_(name)
    {}

    void start()
    {
        start_ = std::chrono::system_clock::now();
    }
    void end()
    {
        end_ = std::chrono::system_clock::now();
    }
    void addTime()
    {
        duration_+=end_-start_;
    }
    double lapTime()
    {
        return Duration(end_-start_).count();
    }
    void printDuration()
    {
        std::cout<<timerName_<<":: "<<duration_.count()<<" sec"<<std::endl;
    }
};

class Statistics
{
public:
    static std::vector<Timer> timers;
          
    static int createTimer(const std::string& name)
    {
        int timerID = (int)timers.size();;
        timers.push_back(Timer(name));
        return timerID;
    }

    static void printStatistics()
    {
        std::cout<<std::endl <<"----------------Timers Start------------------" <<std::endl;
        int nTimers = (int)timers.size();;
        for (int i=0;i<nTimers;i++)
        {
            timers[i].printDuration();
        }
        std::cout<<"----------------Timers End------------------" <<std::endl<<std::endl;
    }
};

#endif
