
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Description
    Utility functions for ranges

SourceFiles
    -

\************************************************************************/

#ifndef ParallelOMPFor_H
#define ParallelOMPFor_H

#include <cstdio>
#include <fstream>
#include <thread>
#include <iostream>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task_scheduler_init.h>

#include "ParallelImpl.hpp"
#include <boost/program_options.hpp>

struct TBBFor : public ParallelImpl<TBBFor>
{
  template <typename Functor>
  static void For(size_t rangeSize, Functor f) {
    auto blocked_f = [f](const tbb::blocked_range<size_t>& rng) {
      size_t iEnd = rng.end();
      for (size_t i = rng.begin(); i<iEnd; ++i)
        f(i);
    };
    const size_t grainSize = std::max(rangeSize/64,size_t(128));
    tbb::parallel_for(
      tbb::blocked_range<size_t>(size_t{0}, rangeSize, grainSize),
      blocked_f,
      tbb::simple_partitioner{}
    );
  }

  static void readSettings() {
    int nThreads = 0;
    std::ifstream Config_File("ParallelConfig.ini");

    using namespace boost::program_options;
    options_description SimSet("Settings");
    SimSet.add_options()
        ("ParallelSettings.nThreads", value<int>(&nThreads));
    
    variables_map vm;
    store(parse_config_file(Config_File, SimSet), vm);
    notify(vm);

    if (!nThreads) {
      nThreads = static_cast<int>(std::thread::hardware_concurrency());
      if (!nThreads) {
        std::cerr<<"Hardware_concurrency method failed!!"<<std::endl;
        nThreads = 1;
      }
    }

    scheduler.initialize(nThreads);
  }

private:
  static tbb::task_scheduler_init scheduler;
};

#endif
