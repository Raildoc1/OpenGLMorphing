#version 330 core

out vec4 FragColor;

in vec3 crntNormal;

uniform vec4 lightColor;
uniform vec3 lightPos;
uniform vec3 camPos;

vec4 direcLight()
{
	float ambient = 0.20f;

	vec3 normal = normalize(crntNormal);
	vec3 lightDirection = normalize(-lightPos);
	float diffuse = clamp(dot(normal, lightDirection), 0.0f, 1.0f);

	return vec4(diffuse + ambient) * lightColor;
}

void main()
{
	FragColor = direcLight();
}