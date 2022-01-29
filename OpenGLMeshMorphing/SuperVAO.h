#include<glad/glad.h>
#include"SuperVBO.h"

class SuperVAO
{
public:
	GLuint ID;
	SuperVAO();

	void LinkAttrib(SuperVBO& VBO, GLuint layout, GLuint numComponents, GLenum type, GLsizeiptr stride, void* offset);
	void Bind();
	void Unbind();
	void Delete();
};

