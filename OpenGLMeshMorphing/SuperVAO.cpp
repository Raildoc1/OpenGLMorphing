#include"SuperVAO.h"

SuperVAO::SuperVAO()
{
	glGenVertexArrays(1, &ID);
}

void SuperVAO::LinkAttrib(SuperVBO& VBO, GLuint layout, GLuint numComponents, GLenum type, GLsizeiptr stride, void* offset)
{
	VBO.Bind();
	glVertexAttribPointer(layout, numComponents, type, GL_FALSE, stride, offset);
	glEnableVertexAttribArray(layout);
	VBO.Unbind();
}

void SuperVAO::Bind()
{
	glBindVertexArray(ID);
}

void SuperVAO::Unbind()
{
	glBindVertexArray(0);
}

void SuperVAO::Delete()
{
	glDeleteVertexArrays(1, &ID);
}