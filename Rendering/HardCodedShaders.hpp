
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Description
    Hardcoded shaders

SourceFiles
    HardCodedShaders.cpp

\************************************************************************/

#ifndef HARDCODEDSHADERS_H
#define HARDCODEDSHADERS_H

#include <string>

class HardCodedShader
{
public:
  static const std::string pointShader_vs;
  static const std::string pointShader_fs;
  static const std::string forceShader_vs;
  static const std::string forceShader_gs;
  static const std::string forceShader_fs;
};

#endif
