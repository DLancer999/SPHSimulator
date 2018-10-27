
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Shader
 
Description
    Class to obscure shader compilation and linking

SourceFiles
    Shader.cpp

\************************************************************************/

#ifndef SHADER_H
#define SHADER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>

// GLEW
#include <GL/glew.h>

class Shader
{
public:
    static std::string stringShaderType (const GLenum& shdrType);

    struct ShaderComponent
    {
        GLuint componentID;
        GLenum componentType;

        ShaderComponent():
        componentID(0),
        componentType(0)
        {}

        ShaderComponent(const GLenum& shdrType):
        componentID(0),
        componentType(shdrType)
        {}

        ~ShaderComponent()
        {
            if (componentID) glDeleteShader(componentID);
        }
    };

    using pShaderComponent = std::unique_ptr<ShaderComponent>;
protected:
    //protected members
    GLuint program_;
    std::vector<pShaderComponent> parts_;

public:
    //constructor
    Shader():
    program_(0), parts_()
    { }

    Shader(const GLchar* computeSourcePath):
    program_(0), parts_()
    {
        compileShaderPart(computeSourcePath, GL_COMPUTE_SHADER);
        linkProgram();
    }

    Shader(const GLchar* vertexSourcePath, const GLchar* fragmentSourcePath):
    program_(0), parts_()
    {
        compileShaderPart(vertexSourcePath,   GL_VERTEX_SHADER);
        compileShaderPart(fragmentSourcePath, GL_FRAGMENT_SHADER);
        linkProgram();
    }

    Shader(const GLchar* vertexSourcePath, const GLchar* geometrySourcePath, const GLchar* fragmentSourcePath):
    program_(0), parts_()
    {
        compileShaderPart(vertexSourcePath,   GL_VERTEX_SHADER);
        compileShaderPart(geometrySourcePath, GL_GEOMETRY_SHADER);
        compileShaderPart(fragmentSourcePath, GL_FRAGMENT_SHADER);
        linkProgram();
    }

    //destructor
    ~Shader()
    {
        if (program_) glDeleteProgram(program_);
    }

    //public member functions
    void compileShaderPart (const std::string& shaderCode, const GLenum& shdrType); //compile each shader
    void linkProgram(); //link shaders and program

    void use() { glUseProgram(program_); }

    GLuint program(){return program_;}

    static std::string readFromFile(const GLchar* sourcePath); //read shader file
};

#endif
