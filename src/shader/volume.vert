#version 330 core

layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 texturePosition;

out vec3 texCoords;

uniform mat4 MVP;

void main(void) {
  gl_Position = MVP *  vec4(vertexPosition,1.0);
  texCoords = texturePosition;
}