
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include "Shader.hpp"

#define ShaderLogging 0
 
//********************************************************************************
std::string Shader::stringShaderType(const GLenum& shdrType)
//********************************************************************************
{
    std::string type;
    switch (shdrType)
    {
        case (GL_VERTEX_SHADER         ): { type = "vertex shader";          break; }
        case (GL_TESS_CONTROL_SHADER   ): { type = "tess control shader";    break; }
        case (GL_TESS_EVALUATION_SHADER): { type = "tess evaluation shader"; break; }
        case (GL_GEOMETRY_SHADER       ): { type = "geometry shader";        break; }
        case (GL_FRAGMENT_SHADER       ): { type = "fragment shader";        break; }
        case (GL_COMPUTE_SHADER        ): { type = "compute shader";         break; }
        default:
        {
            std::cerr<<"Shader::invalid shader type";  
            exit(1); 
        }
    }
    return type;
}

//********************************************************************************
void Shader::compileShaderPart(const GLchar* sourcePath, const GLenum& shdrType)
//********************************************************************************
{
    //stringShaderType also checks for shdrType... throws error for invalid values
#if ShaderLogging
    std::cout<<"Compiling "<<stringShaderType(shdrType)<<" from file \""<<sourcePath<<"\"";
#endif

    pShaderComponent shaderPart(new ShaderComponent(shdrType));

    // 1. Retrieve the shader source code from filePath
    std::string shaderCode;
    std::ifstream shaderFile;
    // ensures ifstream objects can throw exceptions:
    shaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit) ;
    try
    {
        // Open files
        shaderFile.open(sourcePath);
        std::stringstream shaderStream;
        // Read file’s buffer contents into streams
        shaderStream << shaderFile.rdbuf();
        // close file handlers
        shaderFile.close();
        // Convert stream into GLchar array
        shaderCode = shaderStream.str();
    }
    catch(std::ifstream::failure e)
    {
#if ShaderLogging
        std::cout<<" - Failed"<<std::endl;
#endif
        std::cerr << "Shader::"<<stringShaderType(shdrType)<<"::FILE_NOT_SUCCESFULLY_READ" << std::endl;
        exit(1);
    }
    const GLchar* shaderCodeCharPtr = shaderCode.c_str();

    // 2. Compile shaders
    GLint success;
    GLchar infoLog[512];
    // create Shader
    shaderPart->componentID=glCreateShader(shdrType);
    // read source
    glShaderSource(shaderPart->componentID, 1, &shaderCodeCharPtr, NULL);
    // compile source
    glCompileShader(shaderPart->componentID);
    // Print compile errors if any
    glGetShaderiv(shaderPart->componentID, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(shaderPart->componentID, 512, NULL, infoLog);
#if ShaderLogging
        std::cout<<" - Failed"<<std::endl;
#endif
        std::cerr << "Shader::"<<stringShaderType(shdrType)<<"::COMPILATION_FAILED\n" 
                  << infoLog << std::endl;
        exit(1);
    };
#if ShaderLogging
    std::cout<<" - Successful"<<std::endl;
#endif

    // add shader part to list for linking
    parts_.push_back(shaderPart);
}


//********************************************************************************
void Shader::linkProgram()
//********************************************************************************
{
    GLint success;
    GLchar infoLog[512];
    // Shader program
    program_ = glCreateProgram();

    // Link each part
    GLuint nParts = (GLuint)parts_.size();
    for (GLuint iPart=0; iPart<nParts; iPart++)
    {
#if ShaderLogging
        std::cout<<"linking "<<stringShaderType(parts_[iPart]->componentType);
#endif
        glAttachShader(program_, parts_[iPart]->componentID);
#if ShaderLogging
        std::cout<<" - Successful"<<std::endl;
#endif
    }

    // Link to shader program
#if ShaderLogging
    std::cout<<"linking program";
#endif
    glLinkProgram(program_);

    // Print linking errors if any
    glGetProgramiv(program_, GL_LINK_STATUS, &success);
    if (!success)
    {
        glGetProgramInfoLog(program_, 512, NULL, infoLog);
#if ShaderLogging
        std::cout<<" - Failed"<<std::endl;
#endif
        std::cerr << "Shader::program::LINKING_FAILED\n" << infoLog << std::endl;
        exit(1);
    }
#if ShaderLogging
    std::cout<<" - Successful"<<std::endl;
#endif

    // Delete the shader parts -  no longer necessery
    parts_.clear();
}
