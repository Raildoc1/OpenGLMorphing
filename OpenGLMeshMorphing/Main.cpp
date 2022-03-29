#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

//------- Ignore this ----------
#include<filesystem>
namespace fs = std::filesystem;
//------------------------------

#define GLM_ENABLE_EXPERIMENTAL

#include "Model.h"
#include "HarmonicMapper.h"
#include <ctime>
#include <glm/gtx/string_cast.hpp>

#include <Eigen/Dense>

enum class ViewMode { Super, Map };

const unsigned int width = 1000;
const unsigned int height = 1000;

float scale = 0.0f;

#define mode           ViewMode::Super
#define draw_src_map   true
#define draw_tar_map   true
#define draw_super_map true

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	scale += yoffset * 0.1f;
}

GLFWwindow* create_window();
void draw_map(MeshData* src, MeshData* tar, HarmonicMapper* mapper, GLFWwindow* window, Camera& camera);
void draw_super_mesh(HarmonicMapper* mapper, GLFWwindow* window, Camera& camera);

int main() {
	const clock_t begin_time = clock();

	glfwInit();

	std::cout << "glfwInit finished in " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds." << std::endl;

	GLFWwindow* window = create_window();

	if (window == nullptr) {
		return -1;
	}

	std::string parentDir = (fs::current_path().fs::path::parent_path()).string();
	//std::string sourceModelPath = "/Resources/models/pyramide1/pyramide1.gltf";
	//std::string targetModelPath = "/Resources/models/pyramide2/pyramide2.gltf";
	//std::string sourceModelPath = "/Resources/models/pyramide3/pyramide3.gltf";
	//std::string targetModelPath = "/Resources/models/pyramide4/pyramide4.gltf";
	//std::string sourceModelPath = "/Resources/models/tree/tree.gltf";
	//std::string targetModelPath = "/Resources/models/tree2/tree.gltf";

	std::string sourceModelPath = "/Resources/models/Sphere/sphere.gltf";
	std::string targetModelPath = "/Resources/models/Suzanne/suzanne.gltf";

	//std::string sourceModelPath = "/Resources/models/Suzanne1/suzanne1.gltf";
	//std::string targetModelPath = "/Resources/models/Icosphere1/icosphere1.gltf";

	//std::string sourceModelPath = "/Resources/models/suzanne_head/suzanne_head.gltf";
	//std::string targetModelPath = "/Resources/models/human_head/human_head.gltf";

	//std::string sourceModelPath = "/Resources/models/suzanne_head/suzanne_head.gltf";
	//std::string targetModelPath = "/Resources/models/Alberd/alberd_low.gltf";

	Model sourceModel((parentDir + sourceModelPath).c_str());
	Model targetModel((parentDir + targetModelPath).c_str());

	MeshData sourceData = MeshData(sourceModel.GetMesh(), 0.0f, false);
	MeshData targetData = MeshData(targetModel.GetMesh(), 3 * glm::pi<float>() / 4.0f, true);
	//MeshData sourceData = MeshData(sourceModel.GetMesh(), 0.0f, false);
	//MeshData targetData = MeshData(targetModel.GetMesh(), 0.0f, false);

	sourceData.init();
	targetData.init();

	HarmonicMapper mapper(sourceData, targetData);
	Camera camera(width, height, glm::vec3(0.0f, 0.0f, 2.0f));

	if (mode == ViewMode::Super || draw_super_map) {
		mapper.init();
	}

	glfwSetScrollCallback(window, scroll_callback);

	std::cout << "morphing finished in " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds." << std::endl;

	if (mode == ViewMode::Map) {
		draw_map(&sourceData, &targetData, &mapper, window, camera);
	}
	else if (mode == ViewMode::Super) {
		draw_super_mesh(&mapper, window, camera);
	}

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

GLFWwindow* create_window()
{
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(width, height, "Morphing", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return nullptr;
	}
	glfwMakeContextCurrent(window);

	gladLoadGL();
	glViewport(0, 0, width, height);

	return window;
}

void draw_map(MeshData* src, MeshData* tar, HarmonicMapper* mapper, GLFWwindow* window, Camera& camera) {
	Shader debugShader("dot.vert", "dot.frag", false);
	debugShader.Activate();

	int total_edges_amount = 0;
	int src_edges_amount = 0;
	int tar_edges_amount = 0;
	int super_edges_amount = 0;

	if (draw_src_map) {
		src_edges_amount = src->uniqueEdges.size();
		total_edges_amount += src_edges_amount;
	}

	if (draw_tar_map) {
		tar_edges_amount = tar->uniqueEdges.size();
		total_edges_amount += tar_edges_amount;
	}

	if (draw_super_map) {
		super_edges_amount = mapper->uniqueEdges.size();
		total_edges_amount += super_edges_amount;
	}

	GLfloat* map = (GLfloat*)calloc(total_edges_amount * 4, sizeof(GLfloat*));

	int j = 0;

	if (draw_src_map) {
		for (size_t i = 0; i < src_edges_amount; i++, j += 4)
		{
			map[j + 0] = src->map[src->uniqueEdges[i].v1].image.x * 0.5f + 0.5f;
			map[j + 1] = src->map[src->uniqueEdges[i].v1].image.y * 0.5f + 0.5f;
			map[j + 2] = src->map[src->uniqueEdges[i].v2].image.x * 0.5f + 0.5f;
			map[j + 3] = src->map[src->uniqueEdges[i].v2].image.y * 0.5f + 0.5f;
		}
	}

	if (draw_tar_map) {
		for (size_t i = 0; i < tar_edges_amount; i++, j += 4)
		{
			map[j + 0] = tar->map[tar->uniqueEdges[i].v1].image.x * 0.5f - 0.5f;
			map[j + 1] = tar->map[tar->uniqueEdges[i].v1].image.y * 0.5f - 0.5f;
			map[j + 2] = tar->map[tar->uniqueEdges[i].v2].image.x * 0.5f - 0.5f;
			map[j + 3] = tar->map[tar->uniqueEdges[i].v2].image.y * 0.5f - 0.5f;
		}
	}

	if (draw_super_map) {
		for (size_t i = 0; i < super_edges_amount; i++, j += 4)
		{
			map[j + 0] = mapper->map[mapper->uniqueEdges[i].v1].image.x * 0.5f - 0.5f;
			map[j + 1] = mapper->map[mapper->uniqueEdges[i].v1].image.y * 0.5f + 0.5f;
			map[j + 2] = mapper->map[mapper->uniqueEdges[i].v2].image.x * 0.5f - 0.5f;
			map[j + 3] = mapper->map[mapper->uniqueEdges[i].v2].image.y * 0.5f + 0.5f;
		}
	}

	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	GLuint vbo;
	glGenBuffers(1, &vbo);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * total_edges_amount * 4, map, GL_STATIC_DRAW);

	GLint position_attribute = glGetAttribLocation(debugShader.ID, "position");

	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);

	glEnable(GL_PROGRAM_POINT_SIZE);

	glm::vec2 offset = glm::vec2(0.0f, 0.0f);

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		camera.Inputs(window);
		camera.updateMatrix(45.0f, 0.1f, 100.0f);

		float w = glm::exp(-scale);
		offset.x += camera.mouseDeltaX * 0.001f * (1.0f / w);
		offset.y -= camera.mouseDeltaY * 0.001f * (1.0f / w);

		glBindVertexArray(vao);
		glUniform4f(glGetUniformLocation(debugShader.ID, "offset"), offset.x, offset.y, 0.0f, 0.0f);
		glUniform1f(glGetUniformLocation(debugShader.ID, "scale"), w);
		glUniform4f(glGetUniformLocation(debugShader.ID, "baseColor"), 1.0f, 1.0f, 1.0f, 1.0f);
		glDrawArrays(GL_LINES, 0, total_edges_amount * 2);
		glUniform4f(glGetUniformLocation(debugShader.ID, "baseColor"), 1.0f, 0.0f, 0.0f, 1.0f);
		glDrawArrays(GL_POINTS, 0, total_edges_amount * 2);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	debugShader.Delete();
}

void draw_super_mesh(HarmonicMapper* mapper, GLFWwindow* window, Camera& camera) {
	Shader meshShader("default.vert", "default.frag", false);

	glm::vec4 lightColor = glm::vec4(0.43f, 0.91f, 0.85f, 1.0f);
	glm::vec3 lightPos = glm::vec3(0.0f, 0.0f, -5.0f);
	glm::mat4 lightModel = glm::mat4(1.0f);
	lightModel = glm::translate(lightModel, lightPos);

	meshShader.Activate();
	glUniform4f(glGetUniformLocation(meshShader.ID, "lightColor"), lightColor.x, lightColor.y, lightColor.z, lightColor.w);
	glUniform3f(glGetUniformLocation(meshShader.ID, "lightPos"), lightPos.x, lightPos.y, lightPos.z);

	glEnable(GL_DEPTH_TEST);

	Shader superShader("super.vert", "super.frag", true);
	SuperMesh* superMesh = mapper->generateSuperMesh();

	superShader.Activate();

	std::cout << "Super mesh vertices amount = " << superMesh->vertices.size() << std::endl;
	std::cout << "Super mesh indices amount = " << superMesh->indices.size() << std::endl;

	glUniform4f(glGetUniformLocation(superShader.ID, "lightColor"), lightColor.x, lightColor.y, lightColor.z, lightColor.w);
	glUniform3f(glGetUniformLocation(superShader.ID, "lightPos"), lightPos.x, lightPos.y, lightPos.z);

	glEnable(GL_DEPTH_TEST);

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		camera.Inputs(window);
		camera.updateMatrix(45.0f, 0.1f, 100.0f);

		glUniform1f(glGetUniformLocation(superShader.ID, "t"), camera.t);
		(*superMesh).Draw(superShader, camera);
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	meshShader.Delete();
}

//void draw_wired_mesh(HarmonicMapper* mapper, GLFWwindow* window, Camera& camera) {
//	Shader debugShader("dot.vert", "dot.frag", false);
//	debugShader.Activate();
//
//	glUniform4f(glGetUniformLocation(meshShader.ID, "lightColor"), lightColor.x, lightColor.y, lightColor.z, lightColor.w);
//	glUniform3f(glGetUniformLocation(meshShader.ID, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
//
//	glEnable(GL_DEPTH_TEST);
//
//	Shader superShader("super.vert", "super.frag", true);
//	SuperMesh* superMesh = mapper->generateSuperMesh();
//
//	superShader.Activate();
//
//	std::cout << "Super mesh vertices amount = " << superMesh->vertices.size() << std::endl;
//	std::cout << "Super mesh indices amount = " << superMesh->indices.size() << std::endl;
//
//	glUniform4f(glGetUniformLocation(superShader.ID, "lightColor"), lightColor.x, lightColor.y, lightColor.z, lightColor.w);
//	glUniform3f(glGetUniformLocation(superShader.ID, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
//
//	glEnable(GL_DEPTH_TEST);
//
//	while (!glfwWindowShouldClose(window))
//	{
//		glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
//		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//
//		camera.Inputs(window);
//		camera.updateMatrix(45.0f, 0.1f, 100.0f);
//
//		glUniform1f(glGetUniformLocation(superShader.ID, "t"), camera.t);
//		(*superMesh).Draw(superShader, camera);
//		glfwSwapBuffers(window);
//		glfwPollEvents();
//	}
//
//	meshShader.Delete();
//}