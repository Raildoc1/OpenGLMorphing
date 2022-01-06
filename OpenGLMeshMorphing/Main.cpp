//------- Ignore this ----------
#include<filesystem>
namespace fs = std::filesystem;
//------------------------------

#define GLM_ENABLE_EXPERIMENTAL
#include"Model.h"
#include "MeshData.h"
#include <glm/gtx/string_cast.hpp>


const unsigned int width = 800;
const unsigned int height = 800;


int main()
{
	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(width, height, "Morphing", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	gladLoadGL();
	glViewport(0, 0, width, height);

	Shader shaderProgram("default.vert", "default.frag");

	glm::vec4 lightColor = glm::vec4(0.43f, 0.91f, 0.85f, 1.0f);
	glm::vec3 lightPos = glm::vec3(0.5f, 0.5f, 0.5f);
	glm::mat4 lightModel = glm::mat4(1.0f);
	lightModel = glm::translate(lightModel, lightPos);

	shaderProgram.Activate();
	glUniform4f(glGetUniformLocation(shaderProgram.ID, "lightColor"), lightColor.x, lightColor.y, lightColor.z, lightColor.w);
	glUniform3f(glGetUniformLocation(shaderProgram.ID, "lightPos"), lightPos.x, lightPos.y, lightPos.z);

	glEnable(GL_DEPTH_TEST);

	Camera camera(width, height, glm::vec3(0.0f, 0.0f, 2.0f));

	std::string parentDir = (fs::current_path().fs::path::parent_path()).string();
	std::string modelPath = "/Resources/models/tree/tree.gltf";

	Model model((parentDir + modelPath).c_str());

	MeshData data = MeshData(model.GetMesh());
	data.init();

	/*for (size_t i = 0; i < data.getVertexCount(); i++)
	{
		std::cout << data.vertices[i].index << glm::to_string(data.vertices[i].vertex.position) << " -> " << data.vertices[i].eqClass << std::endl;
	}*/

	std::cout << "size = " << data.uniqueEdges.size() << std::endl;

	for (size_t i = 0; i < data.uniqueEdges.size(); i++)
	{
		std::cout << "(" << data.uniqueEdges[i].v1 << ", " << data.uniqueEdges[i].v2 << ")" << std::endl;
	}

	std::cout << "border length = " << data.border.size() << std::endl;

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		camera.Inputs(window);
		camera.updateMatrix(45.0f, 0.1f, 100.0f);

		model.Draw(shaderProgram, camera);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	shaderProgram.Delete();
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}