
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

#ifndef Parallel_H
#define Parallel_H

#include "Definitions.hpp"

#if PARAMIMPL == SERIAL_VAL
  #include "SerialFor.hpp"
#elif PARAMIMPL == OPENMP_VAL
  #include "OMPFor.hpp"
#elif PARAMIMPL == TBB_VAL
  #include "TBBFor.hpp"
#endif

#if PARAMIMPL == SERIAL_VAL
using ForType = SerialFor;
#elif PARAMIMPL == OPENMP_VAL
using ForType = OMPFor;
#elif PARAMIMPL == TBB_VAL
using ForType = TBBFor;
#endif

using Parallel = ParallelImpl<ForType>;
//This static assert forces ForType to inherit privatly from ParallelImpl.
//This is needed to ensure that ParallelImpl interface functions are properly overriden by ForType,
//otherwise a "function is private" compile time error is issued.
static_assert(Parallel::checkInheritance(),"ForType does not inherit privately from ParallelImpl");

#endif
