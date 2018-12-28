
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

#include "Definitions.hpp"

#if PARAMIMPL == TBB_VAL
  #include "tbb/blocked_range.h"
  #include "tbb/parallel_for.h"
  #include "tbb/partitioner.h"
#endif

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
#endif
  }
}

#endif
