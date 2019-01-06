
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
#include <string>
#include <omp.h>

#include "ParallelImpl.hpp"
#include <boost/program_options.hpp>

struct OMPFor : private ParallelImpl<OMPFor>
{
  template <typename Functor>
  static void For(size_t rangeSize, Functor f) {
    #pragma omp parallel for
    for (size_t i=0; i<rangeSize; ++i) {
      f(i);
    }
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

    omp_set_num_threads(nThreads);
  }
};

#endif
