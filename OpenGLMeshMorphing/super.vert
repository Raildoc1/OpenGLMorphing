#version 330 core

layout (location = 0) in vec3 srcPos;
layout (location = 1) in vec3 tarPos;
layout (location = 2) in vec3 srcNorm;
layout (location = 3) in vec3 tarNorm;

out vec3 crntPos;
out vec3 crntNormal;

uniform mat4 camMatrix;
uniform mat4 model;
uniform mat4 translation;
uniform mat4 rotation;
uniform mat4 scale;
uniform float t;

void main()
{
	vec3 intPos = t * tarPos + (1.0f - t) * srcPos;
	vec3 intNorm = t * tarNorm + (1.0f - t) * srcNorm;

	crntPos = vec3(model * translation * -rotation * scale * vec4(intPos, 1.0f));
	crntNormal = intNorm;

	gl_Position = camMatrix * vec4(crntPos, 1.0);
}