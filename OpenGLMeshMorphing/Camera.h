#ifndef CAMERA_CLASS_H
#define CAMERA_CLASS_H

#include<glad/glad.h>
#include<GLFW/glfw3.h>
#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtc/type_ptr.hpp>
#include<glm/gtx/rotate_vector.hpp>
#include<glm/gtx/vector_angle.hpp>

#include"shaderClass.h"

class Camera
{
public:
	// Stores the main vectors of the camera
	glm::vec3 Position;
	glm::vec3 Orientation = glm::vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 Up = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::mat4 cameraMatrix = glm::mat4(1.0f);

	double lastMouseX = 0.0f;
	double lastMouseY = 0.0f;

	double mouseDeltaX = 0.0f;
	double mouseDeltaY = 0.0f;

	float scale = 1.0f;

	float t = 0.0f;
	
	float phi = 0.0f;
	float teta = 0.0f;

	bool enterEverPressed = false;
	bool drawSourceMesh = true;

	// Prevents the camera from jumping around when first clicking left click
	bool firstClick = true;
	bool firstRClick = true;

	// Stores the width and height of the window
	int width;
	int height;

	// Adjust the speed of the camera and it's sensitivity when looking around
	float speed = 0.1f;
	float sensitivity = 100.0f;

	// Camera constructor to set up initial values
	Camera(int width, int height, glm::vec3 position);

	// Updates the camera matrix to the Vertex Shader
	void updateMatrix(float FOVdeg, float nearPlane, float farPlane);
	// Exports the camera matrix to a shader
	void Matrix(Shader& shader, const char* uniform);
	// Handles camera inputs
	void Inputs(GLFWwindow* window);

	glm::vec3 getLightPosition() { return glm::vec3(sin(teta) * sin(phi), cos(teta), sin(teta) * cos(phi)); }
};
#endif