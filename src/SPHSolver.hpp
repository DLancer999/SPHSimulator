#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "Particle.hpp"
#include "HashTable.hpp"
#include "Settings.hpp"

class SPHSolver
{
protected:
    std::vector<Particle> cloud_;
    HashTable             neibhs_;

    double delta_; //used for PCISPH
    int    activeParticles_;

    void generateParticles();
    void calcDensity();
    void calcDensityErr();
    void initPressure();
    void calcPressure();
    void calcPressForces();
    void calcViscForces();
    void calcOtherForces();
    void updatePressure();

    void WCSPHStep();
    void PCISPHStep();
public:
    SPHSolver():
    cloud_(),
    neibhs_(),
    delta_(0.0),
    activeParticles_(0)
    {}

    void init();
    void step(GLFWwindow* window);

    std::vector<Particle>& cloud(){return cloud_;}
    HashTable& neibhs(){return neibhs_;}
};
#endif
