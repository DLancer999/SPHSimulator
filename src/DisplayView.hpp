

/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

namespace
    DisplayView
 
Description
    objects regarding display rendering 

SourceFiles
    DisplayView.cpp

\************************************************************************/

#ifndef DISPLAYVIEW_H
#define DISPLAYVIEW_H

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "Particle.hpp"
#include "Shader.hpp"

namespace DisplayView
{

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mode);

class WindowManager
{
protected:
	GLFWwindow* window_;

    Shader particleShader_;
    Shader forceShader_;
    GLuint VAO_,VBO_p_,VBO_c_,VBO_f_;
    glm::mat4 cameraView_;

    std::vector<glm::vec2> sPosition_;
    std::vector<glm::vec3> sColor_;
    std::vector<glm::vec2> sForce_;

public:
    WindowManager():
	window_(),
    particleShader_(),
    forceShader_(),
    VAO_(0),
    VBO_p_(0),
    VBO_c_(0),
    VBO_f_(0),
    cameraView_(),
    sPosition_(),
    sColor_(),
    sForce_()
    {}

    //disable copy constructor
    WindowManager(const WindowManager& copy);

    ~WindowManager()
    {
        glDeleteBuffers(1, &VBO_p_ );
        glDeleteBuffers(1, &VBO_c_ );
        glDeleteBuffers(1, &VBO_f_ );
        glDeleteVertexArrays(1, &VAO_ );
        
        std::cout<< "destroying window" <<std::endl;
        glfwDestroyWindow(window_);
        
        std::cout<< "terminating glfw" <<std::endl;
        glfwTerminate();
    }
    
    GLFWwindow* window(){return window_;}

    glm::dvec2 calcMaxForce(std::vector<glm::vec2>& force);
    float calcScale(glm::vec2 Fmax);
    void init(std::vector<Particle>& cloud);
    bool renderParticles(std::vector<Particle>& cloud, const int NParticles);

    //disable operator = 
    WindowManager& operator=(const WindowManager& op);
};


}//DisplayView namespace end

#endif
