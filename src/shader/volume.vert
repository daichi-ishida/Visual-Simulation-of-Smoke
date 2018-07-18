#version 330 core

layout(location = 0) in vec3 vertexPosition;

out vec3 texPos;

uniform mat4 MVP;
uniform vec3 ratio;

void main(void) {
  vec3 vPos = 2.0 * vertexPosition - vec3(1.0);
  vPos *= ratio;
  gl_Position = MVP *  vec4(vPos, 1.0);
  texPos = vertexPosition;
}