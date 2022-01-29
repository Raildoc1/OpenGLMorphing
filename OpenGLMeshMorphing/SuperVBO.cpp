#include "SuperVBO.h"

SuperVBO::SuperVBO(std::vector<SuperVertex>& vertices)
{
	glGenBuffers(1, &ID);
	glBindBuffer(GL_ARRAY_BUFFER, ID);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(SuperVertex), vertices.data(), GL_STATIC_DRAW);
}

void SuperVBO::Bind()
{
	glBindBuffer(GL_ARRAY_BUFFER, ID);
}

void SuperVBO::Unbind()
{
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void SuperVBO::Delete()
{
	glDeleteBuffers(1, &ID);
}