#pragma once
#include<glm/glm.hpp>
#include<glad/glad.h>
#include<vector>

struct SuperVertex
{
	glm::vec3 position1;
	glm::vec3 position2;
	glm::vec3 normal1;
	glm::vec3 normal2;
};

class SuperVBO
{
public:
	GLuint ID;
	SuperVBO(std::vector<SuperVertex>& vertices);

	void Bind();
	void Unbind();
	void Delete();
};

