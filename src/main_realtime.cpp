#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

int width = 800;
int height = 600;
GLuint shaderProgram;
float time_val = 0.0f;
float mouseX = 0.0f;
float mouseY = 0.0f;

struct Particle {
  double t, r, theta, phi;
  double dt, dr, dtheta, dphi;
  bool active;
};

std::vector<Particle> particles;
const double M = 1.0;

void integrate_particle(Particle &p, double h) {
  auto get_derivs = [&](const Particle &s, Particle &out) {
    double r = s.r;
    double st = sin(s.theta);
    double ct = cos(s.theta);
    double st_2 = st * st;

    double G_t_tr = M / (r * (r - 2.0 * M));
    double G_r_tt = M * (r - 2.0 * M) / (r * r * r);
    double G_r_rr = -M / (r * (r - 2.0 * M));
    double G_r_thth = -(r - 2.0 * M);
    double G_r_phph = -(r - 2.0 * M) * st_2;
    double G_th_rth = 1.0 / r;
    double G_th_phph = -st * ct;
    double G_ph_rph = 1.0 / r;
    double G_ph_thph = ct / st;

    out.t = s.dt;
    out.r = s.dr;
    out.theta = s.dtheta;
    out.phi = s.dphi;

    out.dt = -2.0 * G_t_tr * s.dt * s.dr;
    out.dr = -(G_r_tt * s.dt * s.dt + G_r_rr * s.dr * s.dr +
               G_r_thth * s.dtheta * s.dtheta + G_r_phph * s.dphi * s.dphi);
    out.dtheta =
        -(2.0 * G_th_rth * s.dr * s.dtheta + G_th_phph * s.dphi * s.dphi);
    out.dphi =
        -(2.0 * G_ph_rph * s.dr * s.dphi + 2.0 * G_ph_thph * s.dtheta * s.dphi);
  };

  Particle k1, k2, k3, k4, temp;

  get_derivs(p, k1);

  temp = p;
  temp.t += k1.t * h * 0.5;
  temp.r += k1.r * h * 0.5;
  temp.theta += k1.theta * h * 0.5;
  temp.phi += k1.phi * h * 0.5;
  temp.dt += k1.dt * h * 0.5;
  temp.dr += k1.dr * h * 0.5;
  temp.dtheta += k1.dtheta * h * 0.5;
  temp.dphi += k1.dphi * h * 0.5;
  get_derivs(temp, k2);

  temp = p;
  temp.t += k2.t * h * 0.5;
  temp.r += k2.r * h * 0.5;
  temp.theta += k2.theta * h * 0.5;
  temp.phi += k2.phi * h * 0.5;
  temp.dt += k2.dt * h * 0.5;
  temp.dr += k2.dr * h * 0.5;
  temp.dtheta += k2.dtheta * h * 0.5;
  temp.dphi += k2.dphi * h * 0.5;
  get_derivs(temp, k3);

  temp = p;
  temp.t += k3.t * h;
  temp.r += k3.r * h;
  temp.theta += k3.theta * h;
  temp.phi += k3.phi * h;
  temp.dt += k3.dt * h;
  temp.dr += k3.dr * h;
  temp.dtheta += k3.dtheta * h;
  temp.dphi += k3.dphi * h;
  get_derivs(temp, k4);

  p.t += h * (k1.t + 2 * k2.t + 2 * k3.t + k4.t) / 6.0;
  p.r += h * (k1.r + 2 * k2.r + 2 * k3.r + k4.r) / 6.0;
  p.theta += h * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta) / 6.0;
  p.phi += h * (k1.phi + 2 * k2.phi + 2 * k3.phi + k4.phi) / 6.0;

  p.dt += h * (k1.dt + 2 * k2.dt + 2 * k3.dt + k4.dt) / 6.0;
  p.dr += h * (k1.dr + 2 * k2.dr + 2 * k3.dr + k4.dr) / 6.0;
  p.dtheta += h * (k1.dtheta + 2 * k2.dtheta + 2 * k3.dtheta + k4.dtheta) / 6.0;
  p.dphi += h * (k1.dphi + 2 * k2.dphi + 2 * k3.dphi + k4.dphi) / 6.0;
}

void spawn_particle() {
  Particle p;
  p.t = 0;
  p.r = 15.0;
  p.theta = 1.57;
  p.phi = 0;

  double v_phi = sqrt(M / p.r) / p.r;

  p.dt = 1.0;
  p.dr = -0.05;
  p.dtheta = 0.01;
  p.dphi = v_phi;
  p.active = true;

  particles.push_back(p);
}

std::string loadShaderSource(const std::string &filename) {
  std::ifstream file(filename);
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

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

  glUniform2f(glGetUniformLocation(shaderProgram, "u_resolution"), (float)width,
              (float)height);
  glUniform1f(glGetUniformLocation(shaderProgram, "u_time"), time_val);
  glUniform2f(glGetUniformLocation(shaderProgram, "u_mouse"), mouseX, mouseY);

  glBegin(GL_QUADS);
  glVertex2f(-1.0f, -1.0f);
  glVertex2f(1.0f, -1.0f);
  glVertex2f(1.0f, 1.0f);
  glVertex2f(-1.0f, 1.0f);
  glEnd();

  glUseProgram(0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, (double)width / height, 0.1, 100.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  double camX = 0, camY = 1.0, camZ = -15.0;
  double rotY = (mouseX / width - 0.5) * 6.0;
  double c = cos(rotY);
  double s = sin(rotY);
  double rx = c * camX + s * camZ;
  double rz = -s * camX + c * camZ;

  gluLookAt(rx, camY, rz, 0, 0, 0, 0, 1, 0);

  glPointSize(5.0);
  glBegin(GL_POINTS);
  glColor3f(1.0, 0.2, 0.2);
  for (auto &p : particles) {
    if (!p.active)
      continue;

    double px = p.r * sin(p.theta) * cos(p.phi);
    double pz = p.r * sin(p.theta) * sin(p.phi);
    double py = p.r * cos(p.theta);

    glVertex3f(px, py, pz);
  }
  glEnd();

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, width, 0, height);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glColor3f(0.0, 0.8, 0.0);
  glBegin(GL_QUADS);
  glVertex2f(20, height - 60);
  glVertex2f(140, height - 60);
  glVertex2f(140, height - 20);
  glVertex2f(20, height - 20);
  glEnd();

  glutSwapBuffers();

  for (auto &p : particles) {
    if (p.active) {
      integrate_particle(p, 0.1);
      if (p.r < 2.0 * M)
        p.active = false;
      if (p.r > 50.0)
        p.active = false;
    }
  }

  time_val += 0.016f;
}

void idle() { glutPostRedisplay(); }

void reshape(int w, int h) {
  width = w;
  height = h;
  glViewport(0, 0, w, h);
}

void mouseMotion(int x, int y) {
  mouseX = (float)x;
  mouseY = (float)y;
}

void mouseClick(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    if (x >= 20 && x <= 140 && y >= 20 && y <= 60) {
      spawn_particle();
      std::cout << "Spawned Particle!" << std::endl;
    }
  }
}

int main(int argc, char **argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(width, height);
  glutCreateWindow("Real-Time Black Hole Simulation");

  initShaders();

  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glutReshapeFunc(reshape);
  glutPassiveMotionFunc(mouseMotion);
  glutMouseFunc(mouseClick);

  glutMainLoop();
  return 0;
}
