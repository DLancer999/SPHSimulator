
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  
 
Description
    global functions to write results either in RAW format or 
    to render to png

SourceFiles
    WriteFunctions.cpp

\************************************************************************/

#ifndef WRITEFUNCTIONS_H
#define WRITEFUNCTIONS_H

#include <iostream>
#include <vector>
#include <glm/glm.hpp>

#include "Particle.hpp"
#include "HashTable.hpp"

void writeRAWfile(std::string fileName, std::vector<Particle>& cloud);
void renderImage(std::string fileName, std::vector<Particle>& cloud, HashTable& neibhs);
void writePNG(std::string fileName, std::vector<glm::vec3>& color);
glm::ivec3 floatToIntColor(glm::vec3 fCol, int res);

#endif
