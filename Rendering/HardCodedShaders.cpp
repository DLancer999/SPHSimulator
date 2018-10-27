
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include "HardCodedShaders.hpp"

const std::string HardCodedShader::pointShader_vs = "\
#version 330 core                                                     \n\
layout (location = 0) in vec2 position;                               \n\
layout (location = 1) in vec3 color;                                  \n\
out vec3 Color;                                                       \n\
uniform mat4 cameraView;                                              \n\
uniform float pointSize;                                              \n\
void main()                                                           \n\
{                                                                     \n\
    gl_Position  = cameraView*vec4(position.x, position.y, 0.0, 1.0); \n\
    //gl_Position  = vec4(position.x, position.y, 0.0, 1.0);          \n\
    gl_PointSize = pointSize;                                         \n\
    Color = color;                                                    \n\
}";

const std::string HardCodedShader::pointShader_fs = "\
#version 330 core                               \n\
in  vec3 Color;                                 \n\
out vec4 color;                                 \n\
void main()                                     \n\
{                                               \n\
    vec2 circCoord = 2.0 * gl_PointCoord - 1.0; \n\
    if (dot(circCoord,circCoord) > 0.98)        \n\
    {                                           \n\
        discard;                                \n\
    }                                           \n\
    color = vec4(Color, 1.0f);                  \n\
}";

const std::string HardCodedShader::forceShader_vs = "\
#version 330 core                                                     \n\
layout (location = 0) in vec2 position;                               \n\
layout (location = 2) in vec2 force;                                  \n\
uniform mat4 cameraView;                                              \n\
out VS_OUT {                                                          \n\
    vec2 Force;                                                       \n\
} vs_out;                                                             \n\
void main()                                                           \n\
{                                                                     \n\
    gl_Position  = cameraView*vec4(position.x, position.y, 0.0, 1.0); \n\
    vs_out.Force = force;                                                           \n\
}";

const std::string HardCodedShader::forceShader_gs = "\
#version 330 core                                                                \n\
layout (points) in;                                                              \n\
layout (line_strip, max_vertices = 2) out;                                       \n\
in VS_OUT {                                                                      \n\
    vec2 Force;                                                                  \n\
} gs_in[];                                                                       \n\
uniform float scale;                                                             \n\
void main()                                                                      \n\
{                                                                                \n\
    gl_Position = gl_in[0].gl_Position;                                          \n\
    EmitVertex();                                                                \n\
    gl_Position = gl_in[0].gl_Position + vec4(gs_in[0].Force, 0.0, 0.0) * scale; \n\
    EmitVertex();                                                                \n\
    EndPrimitive();                                                              \n\
}";

const std::string HardCodedShader::forceShader_fs = "\
#version 330 core              \n\
uniform vec3 Color;            \n\
out vec4 color;                \n\
void main()                    \n\
{                              \n\
    color = vec4(Color, 1.0f); \n\
}";

