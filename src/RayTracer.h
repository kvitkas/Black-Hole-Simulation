#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "Vec3.h"
#include <cmath>
#include <vector>

struct State {
  double t, r, theta, phi;
  double dt, dr, dtheta, dphi;
};

class RayTracer {
public:
  RayTracer(int width, int height);
  void render(const std::string &filename);

private:
  int width, height;
  double M = 1.0; // Mass of the black hole
  double step_size = 0.05;
  double max_dist = 50.0;
  double horizon = 2.0 * M;

  Vec3 get_color(double theta, double phi);
  void integrate(State &s, double h);
  State initial_state(const Vec3 &origin, const Vec3 &direction);

  // Derivatives for RK4
  void derivatives(const State &s, State &out);
};

#endif
