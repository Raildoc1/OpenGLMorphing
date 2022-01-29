#include"SuperVAO.h"
#include"SuperVBO.h"
#include"Camera.h"
#include"Texture.h"
#include"EBO.h"

class SuperMesh
{
public:
	std::vector <SuperVertex> vertices;
	std::vector <GLuint> indices;

	SuperVAO VAO;

	SuperMesh() {}
	SuperMesh(std::vector <SuperVertex>& vertices, std::vector <GLuint>& indices);

	void Draw
	(
		Shader& shader,
		Camera& camera,
		glm::mat4 matrix = glm::mat4(1.0f),
		glm::vec3 translation = glm::vec3(0.0f, 0.0f, 0.0f),
		glm::quat rotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f),
		glm::vec3 scale = glm::vec3(1.0f, 1.0f, 1.0f)
	);
};

