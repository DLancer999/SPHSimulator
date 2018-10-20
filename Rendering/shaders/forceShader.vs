#version 330 core
layout (location = 0) in vec2 position;
layout (location = 2) in vec2 force;

uniform mat4 cameraView;

out VS_OUT {
    vec2 Force;
} vs_out;

void main()
{
    gl_Position  = cameraView*vec4(position.x, position.y, 0.0, 1.0);
    vs_out.Force = force;
}
