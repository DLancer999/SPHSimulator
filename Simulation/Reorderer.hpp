
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Reorderer
 
Description
    Class for reordering particle cloud

SourceFiles
    -

\************************************************************************/

#ifndef REORDERER_H
#define REORDERER_H

#include "ParticleCloud.hpp"
#include "Util/Matrix.hpp"

class Reorderer
{
public:
    static void reorderCloud(ParticleCloud& cloud);
};

#endif

