
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of mandelbrotExplorer.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include "Shader.hpp"
 
//********************************************************************************
Shader::Shader():
//********************************************************************************
program_(0),
vertex_(0),
geometry_(0),
fragment_(0),
compute_(0)
{
}

//********************************************************************************
Shader::Shader(const GLchar* computeSourcePath):
//********************************************************************************
program_(0),
vertex_(0),
geometry_(0),
fragment_(0),
compute_(0)
{
	compileComputeShader(computeSourcePath);
	linkShaders();
}

//********************************************************************************
Shader::Shader(const GLchar* vertexSourcePath, const GLchar* fragmentSourcePath):
//********************************************************************************
program_(0),
vertex_(0),
geometry_(0),
fragment_(0),
compute_(0)
{
	compileVertexShader(vertexSourcePath);
	compileFragmentShader(fragmentSourcePath);
	linkShaders();
}

//********************************************************************************
Shader::Shader(const GLchar* vertexSourcePath, const GLchar* geometrySourcePath, const GLchar* fragmentSourcePath):
//********************************************************************************
program_(0),
vertex_(0),
geometry_(0),
fragment_(0),
compute_(0)
{
	compileVertexShader(vertexSourcePath);
	compileGeometryShader(geometrySourcePath);
	compileFragmentShader(fragmentSourcePath);
	linkShaders();
}

//********************************************************************************
void Shader::compileVertexShader(const GLchar* vertexSourcePath)
//********************************************************************************
{
    std::cout<<"Compiling vertex shader from file::"<<vertexSourcePath<<std::endl;
    // 1. Retrieve the vertex source code from filePath
    std::string vertexCode;
    std::ifstream vShaderFile;
    // ensures ifstream objects can throw exceptions:
    vShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
    try
    {
        // Open files
        vShaderFile.open(vertexSourcePath);
        std::stringstream vShaderStream;
        // Read file’s buffer contents into streams
        vShaderStream << vShaderFile.rdbuf();
        // close file handlers
        vShaderFile.close();
        // Convert stream into GLchar array
        vertexCode = vShaderStream.str();
    }
    catch(std::ifstream::failure e)
    {
        std::cerr << "ERROR::vSHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
        exit(1);
    }
    const GLchar* vShaderCode = vertexCode.c_str();

    // 2. Compile shaders
    GLint success;
    GLchar infoLog[512];
    // Vertex Shader
    vertex_ = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_, 1, &vShaderCode, NULL);
    glCompileShader(vertex_);
    // Print compile errors if any
    glGetShaderiv(vertex_, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(vertex_, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
        exit(1);
    };
}

//********************************************************************************
void Shader::compileGeometryShader(const GLchar* geometrySourcePath)
//********************************************************************************
{
    std::cout<<"Compiling geometry shader from file::"<<geometrySourcePath<<std::endl;
    // 1. Retrieve the fragment source code from filePath
    std::string geometryCode;
    std::ifstream gShaderFile;
    // ensures ifstream objects can throw exceptions:
    gShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit) ;
    try
    {
        // Open files
        gShaderFile.open(geometrySourcePath);
        std::stringstream gShaderStream;
        // Read file’s buffer contents into streams
        gShaderStream << gShaderFile.rdbuf();
        // close file handlers
        gShaderFile.close();
        // Convert stream into GLchar array
        geometryCode = gShaderStream.str();
    }
    catch(std::ifstream::failure e)
    {
        std::cerr << "ERROR::gSHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
        exit(1);
    }
    const GLchar* gShaderCode = geometryCode.c_str();

    // 2. Compile shaders
    GLint success;
    GLchar infoLog[512];
    // geometry Shader
    geometry_ = glCreateShader(GL_GEOMETRY_SHADER);
    glShaderSource(geometry_, 1, &gShaderCode, NULL);
    glCompileShader(geometry_);
    // Print compile errors if any
    glGetShaderiv(geometry_, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(geometry_, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::GEOMETRY::COMPILATION_FAILED\n" << infoLog << std::endl;
        exit(1);
    };
}

//********************************************************************************
void Shader::compileFragmentShader(const GLchar* fragmentSourcePath)
//********************************************************************************
{
    std::cout<<"Compiling fragment shader from file::"<<fragmentSourcePath<<std::endl;
    // 1. Retrieve the fragment source code from filePath
    std::string fragmentCode;
    std::ifstream fShaderFile;
    // ensures ifstream objects can throw exceptions:
    fShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit) ;
    try
    {
        // Open files
        fShaderFile.open(fragmentSourcePath);
        std::stringstream fShaderStream;
        // Read file’s buffer contents into streams
        fShaderStream << fShaderFile.rdbuf();
        // close file handlers
        fShaderFile.close();
        // Convert stream into GLchar array
        fragmentCode = fShaderStream.str();
    }
    catch(std::ifstream::failure e)
    {
        std::cerr << "ERROR::fSHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
        exit(1);
    }
    const GLchar* fShaderCode = fragmentCode.c_str();

    // 2. Compile shaders
    GLint success;
    GLchar infoLog[512];
    // fragment Shader
    fragment_ = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_, 1, &fShaderCode, NULL);
    glCompileShader(fragment_);
    // Print compile errors if any
    glGetShaderiv(fragment_, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(fragment_, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
        exit(1);
    };
}

//********************************************************************************
void Shader::compileComputeShader(const GLchar* computeSourcePath)
//********************************************************************************
{
    std::cout<<"Compiling compute shader from file::"<<computeSourcePath<<std::endl;
    // 1. Retrieve the compute source code from filePath
    std::string computeCode;
    std::ifstream cShaderFile;
    // ensures ifstream objects can throw exceptions:
    cShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit) ;
    try
    {
        // Open files
        cShaderFile.open(computeSourcePath);
        std::stringstream cShaderStream;
        // Read file’s buffer contents into streams
        cShaderStream << cShaderFile.rdbuf();
        // close file handlers
        cShaderFile.close();
        // Convert stream into GLchar array
        computeCode = cShaderStream.str();
    }
    catch(std::ifstream::failure e)
    {
        std::cerr << "ERROR::cSHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
        exit(1);
    }
    const GLchar* cShaderCode = computeCode.c_str();

    // 2. Compile shaders
    GLint success;
    GLchar infoLog[512];
    // fragment Shader
    compute_ = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(compute_, 1, &cShaderCode, NULL);
    glCompileShader(compute_);
    // Print compile errors if any
    glGetShaderiv(compute_, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(compute_, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::COMPUTE::COMPILATION_FAILED\n" << infoLog << std::endl;
        exit(1);
    };
}

//********************************************************************************
void Shader::linkShaders()
//********************************************************************************
{
    GLint success;
    GLchar infoLog[512];
    // Shader program_
    this->program_ = glCreateProgram();
    if (vertex_)   
    {
        std::cout<<"linking vertex shader"<<std::endl;
        glAttachShader(this->program_, vertex_);
    }
    if (geometry_)   
    {
        std::cout<<"linking geometry shader"<<std::endl;
        glAttachShader(this->program_, geometry_);
    }
    if (fragment_) 
    {
        std::cout<<"linking fragment shader"<<std::endl;
        glAttachShader(this->program_, fragment_);
    }
    if (compute_)  
    {
        std::cout<<"linking compute shader"<<std::endl;
        glAttachShader(this->program_, compute_);
    }
    std::cout<<"linking program"<<std::endl;
    glLinkProgram(this->program_);
    // Print linking errors if any
    glGetProgramiv(this->program_, GL_LINK_STATUS, &success);
    if (!success)
    {
        glGetProgramInfoLog(this->program_, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
        exit(1);
    }
    // Delete the shaders as they’re linked into our program now and no longer necessery
    if (vertex_)   glDeleteShader(vertex_);
    if (geometry_)   glDeleteShader(geometry_);
    if (fragment_) glDeleteShader(fragment_);
    if (compute_)  glDeleteShader(compute_);
}

//********************************************************************************
void Shader::use()
//********************************************************************************
{
    glUseProgram(this->program_);
}
