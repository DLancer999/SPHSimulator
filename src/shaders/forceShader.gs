#version 330 core
layout (points) in;
layout (line_strip, max_vertices = 2) out;

in VS_OUT {
    vec2 Force;
} gs_in[];

uniform float scale;

void main()
{
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();
    gl_Position = gl_in[0].gl_Position + vec4(gs_in[0].Force, 0.0, 0.0) * scale;
    EmitVertex();
    EndPrimitive();
}  
