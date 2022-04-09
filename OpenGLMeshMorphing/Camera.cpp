#include"Camera.h"

Camera::Camera(int width, int height, glm::vec3 position)
{
	Camera::width = width;
	Camera::height = height;
	Position = position;
}

void Camera::updateMatrix(float FOVdeg, float nearPlane, float farPlane)
{
	// Initializes matrices since otherwise they will be the null matrix
	glm::mat4 view = glm::mat4(1.0f);
	glm::mat4 projection = glm::mat4(1.0f);

	// Makes camera look in the right direction from the right position
	view = glm::lookAt(Position, Position + Orientation, Up);
	// Adds perspective to the scene
	projection = glm::perspective(glm::radians(FOVdeg), (float)width / height, nearPlane, farPlane);

	// Sets new camera matrix
	cameraMatrix = projection * view;
}

void Camera::Matrix(Shader& shader, const char* uniform)
{
	// Exports camera matrix
	glUniformMatrix4fv(glGetUniformLocation(shader.ID, uniform), 1, GL_FALSE, glm::value_ptr(cameraMatrix));
}

void Camera::Inputs(GLFWwindow* window)
{
	// Handles key inputs
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
	{
		Position += speed * Orientation;
	}
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
	{
		Position += speed * -glm::normalize(glm::cross(Orientation, Up));
	}
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
	{
		Position += speed * -Orientation;
	}
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
	{
		Position += speed * glm::normalize(glm::cross(Orientation, Up));
	}
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
	{
		Position += speed * Up;
	}
	if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
	{
		Position += speed * -Up;
	}
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
	{
		speed = 0.2f;
	}
	else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_RELEASE)
	{
		speed = 0.03f;
	}
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
	{
		t += 0.01f;
		if (t > 1.0f) {
			t = 1.0f;
		}
	}
	if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
	{
		t -= 0.01f;
		if (t < 0.0f) {
			t = 0.0f;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_LEFT))
	{
		phi -= 0.05f;
		if (phi < 0.0f) {
			phi = 6.28f;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_RIGHT))
	{
		phi += 0.05f;
		if (phi > 6.28f) {
			phi = 0.0f;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_UP))
	{
		teta += 0.05f;
		if (teta > 3.14f) {
			teta = 3.14f;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_DOWN))
	{
		teta -= 0.05f;
		if (teta < 0.0f) {
			teta = 0.0f;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_ENTER)) {
		enterEverPressed = true;
	}

	if (glfwGetKey(window, GLFW_KEY_1)) {
		t = 0.0f;
		drawSourceMesh = true;
	}

	if (glfwGetKey(window, GLFW_KEY_2)) {
		t = 0.2f;
		drawSourceMesh = false;
	}

	if (glfwGetKey(window, GLFW_KEY_3)) {
		t = 0.4f;
	}

	if (glfwGetKey(window, GLFW_KEY_4)) {
		t = 0.6f;
	}

	if (glfwGetKey(window, GLFW_KEY_5)) {
		t = 0.8f;
	}

	if (glfwGetKey(window, GLFW_KEY_6)) {
		t = 1.0f;
	}

	// Handles mouse inputs
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
	{
		// Hides mouse cursor
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

		// Prevents camera from jumping on the first click
		if (firstClick)
		{
			glfwSetCursorPos(window, (width / 2), (height / 2));
			glfwGetCursorPos(window, &lastMouseX, &lastMouseY);
			firstClick = false;
		}

		// Stores the coordinates of the cursor
		double mouseX;
		double mouseY;
		// Fetches the coordinates of the cursor
		glfwGetCursorPos(window, &mouseX, &mouseY);

		// Normalizes and shifts the coordinates of the cursor such that they begin in the middle of the screen
		// and then "transforms" them into degrees 
		float rotX = sensitivity * (float)(mouseY - (height / 2)) / height;
		float rotY = sensitivity * (float)(mouseX - (width / 2)) / width;

		// Calculates upcoming vertical change in the Orientation
		glm::vec3 newOrientation = glm::rotate(Orientation, glm::radians(-rotX), glm::normalize(glm::cross(Orientation, Up)));

		// Decides whether or not the next vertical Orientation is legal or not
		if (abs(glm::angle(newOrientation, Up) - glm::radians(90.0f)) <= glm::radians(85.0f))
		{
			Orientation = newOrientation;
		}

		// Rotates the Orientation left and right
		Orientation = glm::rotate(Orientation, glm::radians(-rotY), Up);

		// Sets mouse cursor to the middle of the screen so that it doesn't end up roaming around
		glfwSetCursorPos(window, (width / 2), (height / 2));
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE)
	{
		// Unhides cursor since camera is not looking around anymore
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
		// Makes sure the next time the camera looks around it doesn't jump
		firstClick = true;
	}
	
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
		if (firstRClick) {
			firstRClick = false;
			glfwGetCursorPos(window, &lastMouseX, &lastMouseY);
		}
		else {
			double currentMouseX;
			double currentMouseY;

			glfwGetCursorPos(window, &currentMouseX, &currentMouseY);

			mouseDeltaX = currentMouseX - lastMouseX;
			mouseDeltaY = currentMouseY - lastMouseY;

			lastMouseX = currentMouseX;
			lastMouseY = currentMouseY;
		}
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE)
	{
		firstRClick = true;

		mouseDeltaX = 0.0f;
		mouseDeltaY = 0.0f;
	}
}
