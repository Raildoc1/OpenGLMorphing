#version 150

in vec4 position;

void main() {
	gl_PointSize = 7.5;
	gl_Position = position;
}