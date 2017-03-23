
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of mandelbrotExplorer.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Shader
 
Description
    Class to obscure shader compilation and linking
    Can compile shaders of type
        --vertex
        --fragment
        --compute
 
SourceFiles
    Shader.cpp

\************************************************************************/

#ifndef SHADER_H
#define SHADER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

// GLEW
#include <GL/glew.h>

class Shader
{
protected:
    //protected members
    GLuint program_;
	GLuint vertex_;
	GLuint geometry_;
	GLuint fragment_;
	GLuint compute_;

public:
    //constructor
    Shader();
    Shader(const GLchar* computeSourcePath);
    Shader(const GLchar* vertexSourcePath, const GLchar* fragmentSourcePath);
    Shader(const GLchar* vertexSourcePath, const GLchar* geometrySourcePath, const GLchar* fragmentSourcePath);

    //public member functions
	void compileVertexShader(const GLchar* vertexSourcePath);
	void compileGeometryShader(const GLchar* geometrySourcePath);
	void compileFragmentShader(const GLchar* fragmentSourcePath);
	void compileComputeShader(const GLchar* computeSourcePath);
	void linkShaders();
    void use();
    GLuint program(){return program_;}
};

#endif
