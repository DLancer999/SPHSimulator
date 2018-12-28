
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include "TBBFor.hpp"

tbb::task_scheduler_init TBBFor::scheduler{tbb::task_scheduler_init::deferred};
