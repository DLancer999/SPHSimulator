
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include "Kernels.hpp"
#include <glm/gtx/norm.hpp>

//----------------Kernel Static Variables-----------------//
double Kernel::SmoothingLength::h   = 0.;
double Kernel::SmoothingLength::h2  = 0.;
double Kernel::SmoothingLength::dh  = 0.;
double Kernel::SmoothingLength::dh4 = 0.;
double Kernel::SmoothingLength::dh8 = 0.;
//to be initialized by setSmoothingLength() call

//********************************************************************************
void Kernel::SmoothingLength::setSmoothingLength(const double smoothingLength)
//********************************************************************************
{
    h   = smoothingLength;
    h2  = h*h;
    dh  = 1.0/h;
    dh4 = dh*dh*dh*dh;
    dh8 = dh4*dh4;
}

//********************************************************************************
double Kernel::poly6::W(glm::dvec2& xi, glm::dvec2& xj)
//********************************************************************************
{
    glm::dvec2 rij = xi - xj;
    double len2 = glm::length2(rij);
    if (len2>=SmoothingLength::h2 || len2<0.0) return 0.0;
    else 
    {
        double tmp = SmoothingLength::h2-len2;
        return 4.*glm::one_over_pi<double>()
                 *SmoothingLength::dh8*tmp*tmp*tmp;
    }
}

//********************************************************************************
glm::dvec2 Kernel::poly6::gradW(glm::dvec2& xi, glm::dvec2& xj)
//********************************************************************************
{
    glm::dvec2 rij = xi - xj;
    double len2 = glm::length2(rij);
    if (len2>=SmoothingLength::h2 || len2<0.0) return glm::dvec2(0.0);
    else 
    {
        double tmp = SmoothingLength::h2-len2;
        return -24.*glm::one_over_pi<double>()
               *SmoothingLength::dh8*tmp*tmp*rij;
    }
}

//********************************************************************************
double Kernel::poly6::laplW(glm::dvec2& xi, glm::dvec2& xj)
//********************************************************************************
{
    glm::dvec2 rij = xi - xj;
    double len2 = glm::length2(rij);
    if (len2>=SmoothingLength::h2 || len2<0.0) return 0.0;
    else 
    {
        return -48.*glm::one_over_pi<double>()
                   *SmoothingLength::dh8
                   *(SmoothingLength::h2-len2)
                   *(SmoothingLength::h2-3.0*len2);
    }
}

//********************************************************************************
glm::dvec2 Kernel::spiky::gradW(glm::dvec2& xi, glm::dvec2& xj)
//********************************************************************************
{
    glm::dvec2 rij = xi - xj;
    double q = glm::length(rij)*SmoothingLength::dh;
    if (q>=1. || q<=0.0) return glm::dvec2(0.0);
    else 
    {
        double tmp = 1.-q;
        return -30.*glm::one_over_pi<double>()
                   *SmoothingLength::dh4*tmp*tmp/(q+1.e-6)*rij;
    }
}

//********************************************************************************
double Kernel::visc::laplW(glm::dvec2& xi, glm::dvec2& xj)
//********************************************************************************
{
    glm::dvec2 rij = xi - xj;
    double q = glm::length(rij)*SmoothingLength::dh;
    if (q>=1. || q<0.0) return 0.0;
    else 
    {
        double tmp = 1.-q;
        return 40.*glm::one_over_pi<double>()*SmoothingLength::dh4*tmp;
    }
}
