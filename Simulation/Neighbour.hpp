
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Neighbour
 
Description
    All info needed to describe the relation of two close particles

SourceFiles
    -

\************************************************************************/

#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <glm/glm.hpp>

struct Neigbhor
{
  glm::dvec2 dir;
  double     dist;
  unsigned   ID;

  Neigbhor():
  dir(0.0),
  dist(0.0),
  ID(0)
  {}

  Neigbhor(const unsigned neiID):
  dir(0.0),
  dist(0.0),
  ID(neiID)
  {}
};

#endif
