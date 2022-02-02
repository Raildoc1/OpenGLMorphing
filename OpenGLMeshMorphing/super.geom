#version 330 core
layout ( triangles ) in;
layout ( triangle_strip, max_vertices = 3 ) out;

out vec3 crntNormal;

uniform mat4 camMatrix;

vec3 GetNormal()
{
   mat4 inv = inverse(camMatrix);
   vec3 a = vec3(inv * gl_in[0].gl_Position) - vec3(inv * gl_in[1].gl_Position);
   vec3 b = vec3(inv * gl_in[2].gl_Position) - vec3(inv * gl_in[1].gl_Position);
   return normalize(cross(a, b));
}  

void main(void)
{
    crntNormal = GetNormal();
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();
    gl_Position = gl_in[1].gl_Position;
    EmitVertex();
    gl_Position = gl_in[2].gl_Position;
    EmitVertex();
    EndPrimitive();
}