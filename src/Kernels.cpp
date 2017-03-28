
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
double Kernel::SmoothingLength::dh2 = 0.;
double Kernel::SmoothingLength::dh4 = 0.;
double Kernel::SmoothingLength::dh6 = 0.;
double Kernel::SmoothingLength::dh8 = 0.;
//to be initialized by setSmoothingLength() call

//********************************************************************************
void Kernel::SmoothingLength::setSmoothingLength(const double smoothingLength)
//********************************************************************************
{
    h   = smoothingLength;
    h2  = h*h;
    dh  = 1.0/h;
    dh2 = dh*dh;
    dh4 = dh2*dh2;
    dh6 = dh4*dh2;
    dh8 = dh4*dh4;
}

//********************************************************************************
double Kernel::poly6::W(double mag_rij)
//********************************************************************************
{
    double mag2 = mag_rij*mag_rij;
    if (mag2>=SmoothingLength::h2 || mag2<0.0) return 0.0;
    else 
    {
        double tmp = SmoothingLength::h2-mag2;
        return tmp*tmp*tmp;
    }
}

//********************************************************************************
double Kernel::poly6::W_coeff()
//********************************************************************************
{
    static double coeff = 4.*glm::one_over_pi<double>()*SmoothingLength::dh8;

    return coeff;
}

//********************************************************************************
glm::dvec2 Kernel::poly6::gradW(glm::dvec2& rij, double mag_rij)
//********************************************************************************
{
    double mag2 = mag_rij*mag_rij;
    if (mag2>=SmoothingLength::h2 || mag2<0.0) return glm::dvec2(0.0);
    else 
    {
        double tmp = SmoothingLength::h2-mag2;
        return tmp*tmp*rij;
    }
}

//********************************************************************************
double Kernel::poly6::gradW_coeff()
//********************************************************************************
{
    static double coeff = -24.*glm::one_over_pi<double>()*SmoothingLength::dh8;

    return coeff;
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
        return (SmoothingLength::h2-len2)*(SmoothingLength::h2-3.0*len2);
    }
}

//********************************************************************************
double Kernel::poly6::laplW_coeff()
//********************************************************************************
{
    static double coeff = -48.*glm::one_over_pi<double>()*SmoothingLength::dh8;

    return coeff;
}

//********************************************************************************
glm::dvec2 Kernel::spiky::gradW(glm::dvec2& rij, double mag_rij)
//********************************************************************************
{
    double q = mag_rij*SmoothingLength::dh;
    
    double tmp = 1.-q;
    return tmp*tmp/(q+1.e-6)*rij;
}

//********************************************************************************
double Kernel::spiky::gradW_coeff()
//********************************************************************************
{
    static const double coeff = -30.*glm::one_over_pi<double>()*SmoothingLength::dh4;

    return coeff;
}

//********************************************************************************
double Kernel::visc::laplW(double mag_rij)
//********************************************************************************
{
    double q = mag_rij*SmoothingLength::dh;

    return (1.-q);
}

//********************************************************************************
double Kernel::visc::laplW_coeff()
//********************************************************************************
{
    static double coeff = 40.*glm::one_over_pi<double>()*SmoothingLength::dh4;

    return coeff;
}

//********************************************************************************
double Kernel::surface::C(const double& r)
//********************************************************************************
{
    double q = r*SmoothingLength::dh;
    double tmp = 1.-q;
    double f1 = tmp*tmp*tmp*q*q*q;
    if (q<=0.5)
    {
        f1*=2.;
        f1-=1./64.;
    }
    return f1;
}

//********************************************************************************
double Kernel::surface::C_coeff()
//********************************************************************************
{
    static const double coeff = 32.*glm::one_over_pi<double>()*SmoothingLength::dh2;

    return coeff;
}
