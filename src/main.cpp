#include "RayTracer.h"
#include <iostream>

int main() {
  int width = 800;
  int height = 600;

  std::cout << "Initializing Black Hole Ray Tracer..." << std::endl;
  std::cout << "Image size: " << width << "x" << height << std::endl;

  RayTracer tracer(width, height);
  tracer.render("output/render.ppm");

  return 0;
}
