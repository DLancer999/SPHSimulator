
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

#ifndef ParallelFor_H
#define ParallelFor_H

#include <cstdio>

#define PARAMIMPL 1

#define SERIAL_VAL 0
#define OPENMP_VAL 1
#define TBB_VAL 2

namespace Parallel
{
  template <typename Functor>
  void For(size_t rangeSize, Functor f) {
#if PARAMIMPL == SERIAL_VAL
    for (size_t i=0; i<rangeSize; ++i) {
      f(i);
    }
#elif PARAMIMPL == OPENMP_VAL
    #pragma omp parallel for
    for (size_t i=0; i<rangeSize; ++i) {
      f(i);
    }
#elif PARAMIMPL == TBB_VAL
    static_assert(false, "TBB version not implemented"); 
#endif
  }
}

#endif
