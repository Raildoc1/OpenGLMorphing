//------- Ignore this ----------
#include<filesystem>
namespace fs = std::filesystem;
//------------------------------

#define GLM_ENABLE_EXPERIMENTAL

#include "Model.h"
#include "HarmonicMapper.h"
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
	glm::vec3 lightPos = glm::vec3(0.0f, 0.0f, -5.0f);
	glm::mat4 lightModel = glm::mat4(1.0f);
	lightModel = glm::translate(lightModel, lightPos);

	shaderProgram.Activate();
	glUniform4f(glGetUniformLocation(shaderProgram.ID, "lightColor"), lightColor.x, lightColor.y, lightColor.z, lightColor.w);
	glUniform3f(glGetUniformLocation(shaderProgram.ID, "lightPos"), lightPos.x, lightPos.y, lightPos.z);

	glEnable(GL_DEPTH_TEST);

	Camera camera(width, height, glm::vec3(0.0f, 0.0f, 2.0f));

	std::string parentDir = (fs::current_path().fs::path::parent_path()).string();

	std::string sourceModelPath = "/Resources/models/tree/tree.gltf";
	std::string targetModelPath = "/Resources/models/tree2/tree.gltf";

	Model sourceModel((parentDir + sourceModelPath).c_str());
	Model targetModel((parentDir + targetModelPath).c_str());

	MeshData sourceData = MeshData(sourceModel.GetMesh());
	MeshData targetData = MeshData(targetModel.GetMesh());

	sourceData.init();
	targetData.init();

	Shader debugShader("dot.vert", "dot.frag");
	debugShader.Activate();

	/*for (size_t i = 0; i < data.getVertexCount(); i++)
	{
		std::cout << data.vertices[i].index << glm::to_string(data.vertices[i].vertex.position) << " -> " << data.vertices[i].eqClass << std::endl;
	}*/

	std::cout << "size = " << sourceData.uniqueEdges.size() << std::endl;

	for (size_t i = 0; i < sourceData.uniqueEdges.size(); i++)
	{
		std::cout << "(" << sourceData.uniqueEdges[i].v1 << ", " << sourceData.uniqueEdges[i].v2 << ")" << std::endl;
	}

	std::cout << "border length = " << sourceData.border.size() << std::endl;

	for (size_t i = 0; i < sourceData.border.size(); i++)
	{
		std::cout << "(" << sourceData.border[i].v1.eqClass << ", " << sourceData.border[i].v2.eqClass << ")" << std::endl;
	}

	HarmonicMapper mapper(sourceData, targetData);
	mapper.init();

	int vertexAmount = mapper.sourceMap.size();

	GLfloat* map = (GLfloat*)calloc(vertexAmount * 2, sizeof(GLfloat*));

	for (size_t i = 0; i < vertexAmount * 2; i += 2)
	{
		map[i + 0] = mapper.sourceMap[i / 2].image.x;
		map[i + 1] = mapper.sourceMap[i / 2].image.y;

		std::cout << "(" << map[i + 0] << ", " << map[i + 1] << ")" << std::endl;
	}

	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	GLuint vbo;
	glGenBuffers(1, &vbo);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertexAmount * 2, map, GL_STATIC_DRAW);

	GLint position_attribute = glGetAttribLocation(debugShader.ID, "position");
	
	glVertexAttribPointer(position_attribute, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_attribute);

	glEnable(GL_PROGRAM_POINT_SIZE);

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
		//glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		camera.Inputs(window);
		camera.updateMatrix(45.0f, 0.1f, 100.0f);

		//sourceModel.Draw(shaderProgram, camera);

		glBindVertexArray(vao);
		glDrawArrays(GL_POINTS, 0, vertexAmount);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	shaderProgram.Delete();
	debugShader.Delete();
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}