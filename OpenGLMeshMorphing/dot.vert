#version 150

in vec4 position;

void main() {
	gl_PointSize = 10.0;
	gl_Position = position;
}