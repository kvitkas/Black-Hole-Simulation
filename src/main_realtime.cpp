#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// Global variables
int width = 800;
int height = 600;
GLuint shaderProgram;
float time_val = 0.0f;

// Function to load shader source
std::string loadShaderSource(const std::string &filename) {
  std::ifstream file(filename);
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

// Compile shader
GLuint compileShader(GLenum type, const std::string &source) {
  GLuint shader = glCreateShader(type);
  const char *src = source.c_str();
  glShaderSource(shader, 1, &src, nullptr);
  glCompileShader(shader);

  GLint success;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
  if (!success) {
    char infoLog[512];
    glGetShaderInfoLog(shader, 512, nullptr, infoLog);
    std::cerr << "Shader Compilation Failed:\n" << infoLog << std::endl;
  }
  return shader;
}

void initShaders() {
  std::string vertSource = "#version 120\n"
                           "attribute vec2 position;\n"
                           "void main() {\n"
                           "   gl_Position = vec4(position, 0.0, 1.0);\n"
                           "}\n";

  std::string fragSource = loadShaderSource("src/shaders/blackhole.frag");

  GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertSource);
  GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragSource);

  shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glLinkProgram(shaderProgram);

  GLint success;
  glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
  if (!success) {
    char infoLog[512];
    glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
    std::cerr << "Shader Linking Failed:\n" << infoLog << std::endl;
  }

  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);
}

void display() {
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(shaderProgram);

  // Update uniforms
  glUniform2f(glGetUniformLocation(shaderProgram, "u_resolution"), (float)width,
              (float)height);
  glUniform1f(glGetUniformLocation(shaderProgram, "u_time"), time_val);

  // Draw full screen quad
  glBegin(GL_QUADS);
  glVertex2f(-1.0f, -1.0f);
  glVertex2f(1.0f, -1.0f);
  glVertex2f(1.0f, 1.0f);
  glVertex2f(-1.0f, 1.0f);
  glEnd();

  glutSwapBuffers();

  time_val += 0.016f; // Approx 60 FPS increment
}

void idle() { glutPostRedisplay(); }

void reshape(int w, int h) {
  width = w;
  height = h;
  glViewport(0, 0, w, h);
}

int main(int argc, char **argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(width, height);
  glutCreateWindow("Real-Time Black Hole Simulation");

  initShaders();

  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glutReshapeFunc(reshape);

  glutMainLoop();
  return 0;
}
