
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
#include <functional>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <vector>
#include <algorithm>

class Timer
{
public:
    using Clock      = std::chrono::high_resolution_clock;
    using TimePoint  = std::chrono::time_point<Clock>;
    using Duration   = std::chrono::duration<double>;

private:
    TimePoint start_;
    TimePoint end_;
    Duration duration_;
    std::string timerName_;
public:
    Timer()
    : start_(Clock::now())
    , end_(start_)
    , duration_(end_-start_)
    , timerName_("")
    {}

    Timer(const std::string& name)
    : start_(Clock::now())
    , end_(start_)
    , duration_(end_-start_)
    , timerName_(name)
    {}

    Timer(const Timer& other) = default;

    void start()
    {
        start_ = Clock::now();
    }
    void end()
    {
        end_ = Clock::now();
    }
    void addTime()
    {
        duration_+=end_-start_;
    }
    double lapTime() const
    {
        return Duration(end_-start_).count();
    }
    double totalTime() const
    {
        return duration_.count();
    }
    const std::string& name() const
    {
        return timerName_;
    }
    void printDuration() const
    {
        std::cout<<timerName_<<":: "<<duration_.count()<<" sec"<<std::endl;
    }
};

class Statistics
{
public:
    struct TimerID
    {
        size_t index;
    };

    class TimerGuard
    {
    public:
        TimerGuard(TimerID timerID)
        : _timerID(timerID)
        {
          timers[_timerID.index].start();
        }

        ~TimerGuard() {
          timers[_timerID.index].end();
          timers[_timerID.index].addTime();
        }
    private:
        TimerID _timerID;
    };

    template <typename Functor>
    class CountTime {
    public:
      TimerID _id;
      Functor _f;

      CountTime(TimerID id, Functor f):_id(id),_f(f){}

      template <typename ... Args>
      std::invoke_result_t<Functor, Args...>
      operator()(Args&&... args)
      {
          TimerGuard g(_id);
          return std::invoke(_f, std::forward<Args>(args)...);
      }
    };

          
    static TimerID createTimer(const std::string& name)
    {
        TimerID timerID{timers.size()};
        timers.push_back(Timer(name));
        return timerID;
    }

    static void printStatistics(std::ostream& os = std::cout)
    {
        if (timers.empty())
          return;

        std::vector<Timer> sortedTimers = timers;
        std::sort(sortedTimers.begin(), sortedTimers.end(),
            [](const Timer& a, const Timer& b){ return a.totalTime() > b.totalTime(); }
        );

        auto it = std::max_element( sortedTimers.begin(), sortedTimers.end(),
            [](const Timer& t1, const Timer& t2){ return t1.name().size() < t2.name().size(); }
        );
        size_t namesWidth = (*it).name().size();

        for (const Timer& t : sortedTimers) {
            os<<std::left<<std::setw(int(namesWidth))<<t.name()<<" ::";
            os<<std::right<<std::fixed<<std::setprecision(6)
              <<std::setw(11)<<t.totalTime()<<" sec"<<'\n';
        }
    }

private:
    static std::vector<Timer> timers;
};

#endif
