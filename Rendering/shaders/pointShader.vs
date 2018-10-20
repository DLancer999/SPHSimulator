#version 330 core
layout (location = 0) in vec2 position;
layout (location = 1) in vec3 color;
out vec3 Color;
uniform mat4 cameraView;
uniform float pointSize;
void main()
{
    gl_Position  = cameraView*vec4(position.x, position.y, 0.0, 1.0);
    //gl_Position  = vec4(position.x, position.y, 0.0, 1.0);
    gl_PointSize = pointSize;
    Color = color;
}
