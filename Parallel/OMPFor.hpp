
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

#include "ParallelImpl.hpp"
#include <boost/program_options.hpp>

struct OMPFor : public ParallelImpl<OMPFor>
{
  template <typename Functor>
  static void For(size_t rangeSize, Functor f) {
    #pragma omp parallel for
    for (size_t i=0; i<rangeSize; ++i) {
      f(i);
    }
  }

  //static void readSettings() {
  //  size_t nThreads;
  //  std::ifstream Config_File("ParallelConfig.ini");

  //  using namespace boost::program_options;
  //  options_description SimSet("Settings");
  //  SimSet.add_options()
  //      ("ParallelSettings.nThreads", value<size_t>(&nThreads));
  //  
  //  variables_map vm;
  //  store(parse_config_file(Config_File, SimSet), vm);
  //  notify(vm);
  //}
};

#endif
