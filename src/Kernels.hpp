
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

namespace
    Kernel
 
Description
    definition of smoothing length and different kernels

SourceFiles
    Kernels.cpp

\************************************************************************/

#ifndef KERNELS_H
#define KERNELS_H

#include <glm/glm.hpp>

namespace Kernel
{
    class SmoothingLength
    {
    public:
        static double h  ; //smoothing length
        static double h2 ; //h*h
        static double dh ; //1./h
        static double dh2; //1./(pow(h,2))
        static double dh4; //1./(pow(h,4))
        static double dh6; //1./(pow(h,6))
        static double dh8; //1./(pow(h,8))

        static void setSmoothingLength(const double);
    };

    namespace poly6
    {
        double         W(glm::dvec2& xi, glm::dvec2& xj);
        double         W_coeff();
        glm::dvec2 gradW(glm::dvec2& xi, glm::dvec2& xj);
        double     gradW_coeff();
        double     laplW(glm::dvec2& xi, glm::dvec2& xj);
        double     laplW_coeff();
    }
    namespace spiky
    {
        glm::dvec2 gradW(glm::dvec2& xi, glm::dvec2& xj);
        double     gradW_coeff();
    }
    namespace visc
    {
        double laplW(glm::dvec2& xi, glm::dvec2& xj);
        double laplW_coeff();
    }
    namespace surface
    {
        double C(const double& r);
        double C_coeff();
    }
}

#endif
