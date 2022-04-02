#pragma once
#include<glm/glm.hpp>
#include<glad/glad.h>
#include<vector>

struct SuperVertex
{
	glm::vec3 position1;
	glm::vec3 position2;

	SuperVertex(glm::vec3 srcPos, glm::vec3 tarPos) :
		position1(srcPos),
		position2(tarPos)
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

