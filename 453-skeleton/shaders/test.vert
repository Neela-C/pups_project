#version 330 core
layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 color;
layout (location = 2) in vec3 normal;


out vec3 fragPos;
out vec3 fragColor;
out vec3 n;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;


void main() {
	fragPos = pos;
	fragColor = color;
	mat4 M =  model * view * projection;
	n = normalize(vec4( normal , 0.f) *M).xyz;
	gl_Position = vec4(pos, 1.0)* M;
	
}
