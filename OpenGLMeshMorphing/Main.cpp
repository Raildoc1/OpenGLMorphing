//------- Ignore this ----------
#include<filesystem>
namespace fs = std::filesystem;
//------------------------------

#define GLM_ENABLE_EXPERIMENTAL

#include "Model.h"
#include "HarmonicMapper.h"
#include <glm/gtx/string_cast.hpp>

enum class ViewMode { Super, Map };

const unsigned int width = 1000;
const unsigned int height = 1000;

float scale = 0.0f;

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	scale += yoffset * 0.1f;
}

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

	Shader shaderProgram("default.vert", "default.frag", false);

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

	//std::string sourceModelPath = "/Resources/models/pyramide1/pyramide1.gltf";
	//std::string targetModelPath = "/Resources/models/pyramide2/pyramide2.gltf";
	//std::string sourceModelPath = "/Resources/models/pyramide3/pyramide3.gltf";
	//std::string targetModelPath = "/Resources/models/pyramide4/pyramide4.gltf";
	std::string sourceModelPath = "/Resources/models/tree/tree.gltf";
	std::string targetModelPath = "/Resources/models/tree2/tree.gltf";
	//std::string sourceModelPath = "/Resources/models/Sphere/sphere.gltf";
	//std::string targetModelPath = "/Resources/models/Suzanne/suzanne.gltf";

	Model sourceModel((parentDir + sourceModelPath).c_str());
	Model targetModel((parentDir + targetModelPath).c_str());

	MeshData sourceData = MeshData(sourceModel.GetMesh());
	MeshData targetData = MeshData(targetModel.GetMesh());

	sourceData.init();
	targetData.init();

	//sourceData.harmonizeMap();
	//targetData.harmonizeMap();

	Shader debugShader("dot.vert", "dot.frag", false);
	debugShader.Activate();

	/*for (size_t i = 0; i < data.getVertexCount(); i++)
	{
		std::cout << data.vertices[i].index << glm::to_string(data.vertices[i].vertex.position) << " -> " << data.vertices[i].eqClass << std::endl;
	}*/

	/*std::cout << "size = " << sourceData.uniqueEdges.size() << std::endl;

	for (size_t i = 0; i < sourceData.uniqueEdges.size(); i++)
	{
		std::cout << "(" << sourceData.uniqueEdges[i].v1 << ", " << sourceData.uniqueEdges[i].v2 << ")" << std::endl;
	}

	std::cout << "border length = " << sourceData.border.size() << std::endl;

	for (size_t i = 0; i < sourceData.border.size(); i++)
	{
		std::cout << "(" << sourceData.border[i].v1.eqClass << ", " << sourceData.border[i].v2.eqClass << ")" << std::endl;
	}*/

	HarmonicMapper mapper(sourceData, targetData);
	mapper.init();

	int sourceEdgesAmount = sourceData.uniqueEdges.size();
	int targetEdgesAmount = targetData.uniqueEdges.size();
	int superEdgesAmount = mapper.uniqueEdges.size();
	int totalEdgesAmount = sourceEdgesAmount + targetEdgesAmount + superEdgesAmount;

	GLfloat* map = (GLfloat*)calloc(totalEdgesAmount * 4, sizeof(GLfloat*));

	int j = 0;

	for (size_t i = 0; i < sourceEdgesAmount; i++, j += 4)
	{
		map[j + 0] = sourceData.map[sourceData.uniqueEdges[i].v1].image.x * 0.5f + 0.5f;
		map[j + 1] = sourceData.map[sourceData.uniqueEdges[i].v1].image.y * 0.5f + 0.5f;
		map[j + 2] = sourceData.map[sourceData.uniqueEdges[i].v2].image.x * 0.5f + 0.5f;
		map[j + 3] = sourceData.map[sourceData.uniqueEdges[i].v2].image.y * 0.5f + 0.5f;
	}

	for (size_t i = 0; i < targetEdgesAmount; i++, j += 4)
	{
		map[j + 0] = targetData.map[targetData.uniqueEdges[i].v1].image.x * 0.5f - 0.5f;
		map[j + 1] = targetData.map[targetData.uniqueEdges[i].v1].image.y * 0.5f - 0.5f;
		map[j + 2] = targetData.map[targetData.uniqueEdges[i].v2].image.x * 0.5f - 0.5f;
		map[j + 3] = targetData.map[targetData.uniqueEdges[i].v2].image.y * 0.5f - 0.5f;
	}

	for (size_t i = 0; i < superEdgesAmount; i++, j += 4)
	{
		map[j + 0] = mapper.map[mapper.uniqueEdges[i].v1].image.x * 0.5f - 0.5f;
		map[j + 1] = mapper.map[mapper.uniqueEdges[i].v1].image.y * 0.5f + 0.5f;
		map[j + 2] = mapper.map[mapper.uniqueEdges[i].v2].image.x * 0.5f - 0.5f;
		map[j + 3] = mapper.map[mapper.uniqueEdges[i].v2].image.y * 0.5f + 0.5f;
	}

	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	GLuint vbo;
	glGenBuffers(1, &vbo);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * totalEdgesAmount * 4, map, GL_STATIC_DRAW);

	GLint position_attribute = glGetAttribLocation(debugShader.ID, "position");
	
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);

	glEnable(GL_PROGRAM_POINT_SIZE);

	glm::vec2 offset = glm::vec2(0.0f, 0.0f);

	glfwSetScrollCallback(window, scroll_callback);

	// Super Mesh

	ViewMode mode = ViewMode::Super;

	Shader superShader("super.vert", "super.frag", true);
	SuperMesh* superMesh = mapper.generateSuperMesh();;

	if (mode == ViewMode::Super) {
		superShader.Activate();

		std::cout << "Super mesh vertices amount = " << superMesh->vertices.size() << std::endl;
		std::cout << "Super mesh indices amount = " << superMesh->indices.size() << std::endl;

		glUniform4f(glGetUniformLocation(superShader.ID, "lightColor"), lightColor.x, lightColor.y, lightColor.z, lightColor.w);
		glUniform3f(glGetUniformLocation(superShader.ID, "lightPos"), lightPos.x, lightPos.y, lightPos.z);

		glEnable(GL_DEPTH_TEST);
	}

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
		//glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		camera.Inputs(window);
		camera.updateMatrix(45.0f, 0.1f, 100.0f);

		//std::cout << to_string(offset) << " : " << camera.mouseDeltaX << ", " << camera.mouseDeltaY << std::endl;
		//sourceModel.Draw(shaderProgram, camera);

		if (mode == ViewMode::Map) {
			float w = glm::exp(-scale);
			offset.x += camera.mouseDeltaX * 0.001f * (1.0f / w);
			offset.y -= camera.mouseDeltaY * 0.001f * (1.0f / w);

			glBindVertexArray(vao);
			glUniform4f(glGetUniformLocation(debugShader.ID, "offset"), offset.x, offset.y, 0.0f, 0.0f);
			glUniform1f(glGetUniformLocation(debugShader.ID, "scale"), w);
			glUniform4f(glGetUniformLocation(debugShader.ID, "baseColor"), 1.0f, 1.0f, 1.0f, 1.0f);
			glDrawArrays(GL_LINES, 0, totalEdgesAmount * 2);
			glUniform4f(glGetUniformLocation(debugShader.ID, "baseColor"), 1.0f, 0.0f, 0.0f, 1.0f);
			glDrawArrays(GL_POINTS, 0, totalEdgesAmount * 2);
		} else if (mode == ViewMode::Super) {
			glUniform1f(glGetUniformLocation(superShader.ID, "t"), camera.t);
			(*superMesh).Draw(superShader, camera);
		}

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	shaderProgram.Delete();
	debugShader.Delete();
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}