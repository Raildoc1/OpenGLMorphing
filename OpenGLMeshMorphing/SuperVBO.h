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

	SuperVertex(glm::vec3 srcPos, glm::vec3 tarPos, glm::vec3 srcNorm, glm::vec3 tarNorm) :
		position1(srcPos),
		position2(tarPos),
		normal1(srcNorm),
		normal2(tarNorm)
	{ }
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

