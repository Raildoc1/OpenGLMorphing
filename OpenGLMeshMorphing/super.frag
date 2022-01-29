#version 330 core

out vec4 FragColor;

in vec3 crntPos;
in vec3 crntNormal;

uniform vec4 lightColor;
uniform vec3 lightPos;
uniform vec3 camPos;

vec4 direcLight()
{
	float ambient = 0.20f;

	vec3 normal = normalize(Normal);
	vec3 lightDirection = normalize(vec3(-3.0f, -5.0f, -2.0f));
	float diffuse = max(dot(normal, lightDirection), 0.0f);

	return vec4(diffuse + ambient) * lightColor;
}

void main()
{
	FragColor = direcLight();
}