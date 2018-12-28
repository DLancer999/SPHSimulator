
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

#ifndef ParallelImpl_H
#define ParallelImpl_H

#include <cstdio>
#include <utility>

template <typename Base>
struct ParallelImpl
{
  template <typename Functor>
  static void For(size_t rangeSize, Functor&& f) {
    Base::For(rangeSize, std::forward<Functor>(f));
  }

  //static void readSettings() {
  //  Base::readSettings();
  //}
};

#endif
