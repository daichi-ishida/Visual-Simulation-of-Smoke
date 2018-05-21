layout(location = 0) in vec3 vertexPosition;
uniform mat4 MVP;
 
void main(void) {
  gl_Position = MVP * vertexPosition;
}