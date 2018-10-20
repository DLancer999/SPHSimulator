#version 330 core
in  vec3 Color;
out vec4 color;
void main()
{
    vec2 circCoord = 2.0 * gl_PointCoord - 1.0;
    if (dot(circCoord,circCoord) > 0.98)
    {
        discard;
    }
    color = vec4(Color, 1.0f);
}
