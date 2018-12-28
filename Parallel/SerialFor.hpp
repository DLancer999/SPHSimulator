
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

#include "ParallelImpl.hpp"

struct SerialFor : public ParallelImpl<SerialFor>
{
  template <typename Functor>
  static void For(size_t rangeSize, Functor f) {
    for (size_t i=0; i<rangeSize; ++i) {
      f(i);
    }
  }

  static void readSettings() {
    //nothing to setup
  }
};

#endif
