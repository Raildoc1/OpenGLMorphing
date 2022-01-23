#version 150

in vec4 position;
uniform vec4 offset;
uniform float scale;

void main() {
	gl_PointSize = 7.5;
	vec4 pos = position + offset;
	gl_Position = vec4(pos.x, pos.y, pos.z, 1f / scale);
}