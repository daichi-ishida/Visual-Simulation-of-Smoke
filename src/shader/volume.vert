#version 330 core

layout(location = 0) in vec3 texturePosition;

out vec3 vertPos;
out vec3 texPos;

uniform mat4 MVP;
uniform vec3 ratio;

void main(void) {
  vec3 pos = 2.0 * texturePosition - vec3(1.0);
  pos = ratio.xyz * pos.xyz;
  gl_Position = MVP *  vec4(pos, 1.0);
  vertPos = pos;
  texPos = texturePosition;
}