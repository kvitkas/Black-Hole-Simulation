#include "RayTracer.h"
#include "ImageWriter.h"
#include <cmath>
#include <iostream>

RayTracer::RayTracer(int width, int height) : width(width), height(height) {}

void RayTracer::render(const std::string &filename) {
  std::vector<unsigned char> image(width * height * 3);

  // Camera setup
  Vec3 cam_pos(0, 0, -15); // Start far away on Z axis
  Vec3 cam_lookat(0, 0, 0);
  Vec3 cam_up(0, 1, 0);

  Vec3 forward = unit_vector(cam_lookat - cam_pos);
  Vec3 right = unit_vector(cross(forward, cam_up));
  Vec3 up = cross(right, forward);

  double fov = 60.0;
  double scale = tan((fov * 0.5 * M_PI / 180.0));
  double aspect_ratio = (double)width / height;

  std::cout << "Rendering..." << std::endl;

#pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      // Normalized device coordinates
      double x = (2 * (i + 0.5) / (double)width - 1) * aspect_ratio * scale;
      double y = (1 - 2 * (j + 0.5) / (double)height) * scale;

      Vec3 dir = unit_vector(forward + x * right + y * up);

      State s = initial_state(cam_pos, dir);

      bool hit_horizon = false;
      // Integration loop
      for (int step = 0; step < 5000; ++step) {
        integrate(s, step_size);

        if (s.r < horizon * 1.01) {
          hit_horizon = true;
          break;
        }
        if (s.r > max_dist) {
          break;
        }
      }

      Vec3 color(0, 0, 0);
      if (hit_horizon) {
        color = Vec3(0, 0, 0); // Black hole shadow
      } else {
        // Map final direction to background
        // We need to convert final spherical coords back to direction vector to
        // sample a skybox Or just use theta/phi directly for a procedural
        // pattern
        color = get_color(s.theta, s.phi);
      }

      int idx = (j * width + i) * 3;
      image[idx] = static_cast<unsigned char>(std::min(255.0, color.x * 255));
      image[idx + 1] =
          static_cast<unsigned char>(std::min(255.0, color.y * 255));
      image[idx + 2] =
          static_cast<unsigned char>(std::min(255.0, color.z * 255));
    }
    if (j % 10 == 0)
      std::cout << "Row " << j << " / " << height << "\r" << std::flush;
  }
  std::cout << std::endl;

  ImageWriter::writePPM(filename, width, height, image);
}

Vec3 RayTracer::get_color(double theta, double phi) {
  // Simple procedural background: Checkerboard or Grid
  // Map theta/phi to UV
  double u = phi / (2 * M_PI);
  double v = theta / M_PI;

  // Galaxy band
  double galaxy = std::exp(-10.0 * (v - 0.5) * (v - 0.5));

  // Grid lines
  bool grid = (int(u * 20) % 2 == 0) ^ (int(v * 10) % 2 == 0);

  Vec3 base = grid ? Vec3(0.1, 0.1, 0.2) : Vec3(0.05, 0.05, 0.1);

  // Add galaxy band
  base += Vec3(0.8, 0.7, 0.5) * galaxy;

  return base;
}

State RayTracer::initial_state(const Vec3 &origin, const Vec3 &direction) {
  State s;

  // Convert Cartesian (x,y,z) to Spherical (r, theta, phi)
  // Assuming origin is in standard Cartesian where Z is axis?
  // Actually, let's align physics coordinates with standard math:
  // x = r sin theta cos phi
  // y = r sin theta sin phi
  // z = r cos theta

  double r = origin.length();
  double theta = acos(origin.z / r);
  double phi = atan2(origin.y, origin.x);

  s.t = 0;
  s.r = r;
  s.theta = theta;
  s.phi = phi;

  // Convert velocity vector (direction) to spherical basis
  // We need to project 'direction' onto e_r, e_theta, e_phi

  // Basis vectors at position (r, theta, phi)
  double st = sin(theta);
  double ct = cos(theta);
  double sp = sin(phi);
  double cp = cos(phi);

  Vec3 e_r(st * cp, st * sp, ct);
  Vec3 e_theta(ct * cp, ct * sp, -st);
  Vec3 e_phi(-sp, cp, 0);

  // Local derivatives
  s.dr = dot(direction, e_r);
  s.dtheta = dot(direction, e_theta) / r;
  s.dphi = dot(direction, e_phi) / (r * st);

  // Initial dt/dlambda for a photon
  // Null geodesic condition: ds^2 = 0
  // -(1-2M/r)dt^2 + (1-2M/r)^-1 dr^2 + r^2 dtheta^2 + r^2 sin^2theta dphi^2 = 0
  // (1-2M/r)dt^2 = (1-2M/r)^-1 dr^2 + r^2 (dtheta^2 + sin^2theta dphi^2)
  // dt = sqrt( (1-2M/r)^-2 dr^2 + (1-2M/r)^-1 r^2 (...) )

  double metric_rr = 1.0 / (1.0 - 2.0 * M / r);
  double metric_ang = r * r * (s.dtheta * s.dtheta + st * st * s.dphi * s.dphi);
  double rhs = metric_rr * s.dr * s.dr + metric_ang;

  // (1-2M/r) dt^2 = rhs
  // dt^2 = rhs / (1-2M/r)

  double factor = (1.0 - 2.0 * M / r);
  s.dt = sqrt(rhs / factor) / sqrt(factor); // Wait, check algebra

  // (1-2M/r) dt^2 = (1-2M/r)^-2 dr^2 + (1-2M/r)^-1 r^2 dOmega^2
  // dt^2 = (1-2M/r)^-2 dr^2 + (1-2M/r)^-1 r^2 dOmega^2
  // dt = sqrt(...)

  s.dt = sqrt((s.dr * s.dr) / pow(factor, 2) +
              (r * r * (s.dtheta * s.dtheta + st * st * s.dphi * s.dphi)) /
                  factor);

  return s;
}

void RayTracer::derivatives(const State &s, State &out) {
  // Geodesic equations for Schwarzschild metric
  // \ddot{x}^\mu = -\Gamma^\mu_{\alpha\beta} \dot{x}^\alpha \dot{x}^\beta

  double r = s.r;
  double theta = s.theta;

  double r_2 = r * r;
  double r_3 = r_2 * r;
  double st = sin(theta);
  double ct = cos(theta);
  double st_2 = st * st;

  // Non-zero Christoffel symbols (M=1 assumed in formulas, but using variable
  // M) G^t_tr = M / (r(r-2M)) G^r_tt = M(r-2M) / r^3 G^r_rr = -M / (r(r-2M))
  // G^r_thth = -(r-2M)
  // G^r_phph = -(r-2M)sin^2theta
  // G^th_rth = 1/r
  // G^th_phph = -sin(theta)cos(theta)
  // G^ph_rph = 1/r
  // G^ph_thph = cot(theta)

  double G_t_tr = M / (r * (r - 2.0 * M));
  double G_r_tt = M * (r - 2.0 * M) / r_3;
  double G_r_rr = -M / (r * (r - 2.0 * M));
  double G_r_thth = -(r - 2.0 * M);
  double G_r_phph = -(r - 2.0 * M) * st_2;
  double G_th_rth = 1.0 / r;
  double G_th_phph = -st * ct;
  double G_ph_rph = 1.0 / r;
  double G_ph_thph = ct / st;

  // \ddot{t} = -2 * G^t_tr * dt * dr
  double ddt = -2.0 * G_t_tr * s.dt * s.dr;

  // \ddot{r} = - (G^r_tt dt^2 + G^r_rr dr^2 + G^r_thth dth^2 + G^r_phph dph^2)
  double ddr = -(G_r_tt * s.dt * s.dt + G_r_rr * s.dr * s.dr +
                 G_r_thth * s.dtheta * s.dtheta + G_r_phph * s.dphi * s.dphi);

  // \ddot{theta} = - (2 * G^th_rth * dr * dth + G^th_phph * dph^2)
  double ddtheta =
      -(2.0 * G_th_rth * s.dr * s.dtheta + G_th_phph * s.dphi * s.dphi);

  // \ddot{phi} = - (2 * G^ph_rph * dr * dph + 2 * G^ph_thph * dth * dph)
  double ddphi =
      -(2.0 * G_ph_rph * s.dr * s.dphi + 2.0 * G_ph_thph * s.dtheta * s.dphi);

  out.t = s.dt;
  out.r = s.dr;
  out.theta = s.dtheta;
  out.phi = s.dphi;

  out.dt = ddt;
  out.dr = ddr;
  out.dtheta = ddtheta;
  out.dphi = ddphi;
}

void RayTracer::integrate(State &s, double h) {
  State k1, k2, k3, k4;
  State temp;

  // k1
  derivatives(s, k1);

  // k2
  temp = s;
  temp.t += k1.t * h * 0.5;
  temp.r += k1.r * h * 0.5;
  temp.theta += k1.theta * h * 0.5;
  temp.phi += k1.phi * h * 0.5;
  temp.dt += k1.dt * h * 0.5;
  temp.dr += k1.dr * h * 0.5;
  temp.dtheta += k1.dtheta * h * 0.5;
  temp.dphi += k1.dphi * h * 0.5;
  derivatives(temp, k2);

  // k3
  temp = s;
  temp.t += k2.t * h * 0.5;
  temp.r += k2.r * h * 0.5;
  temp.theta += k2.theta * h * 0.5;
  temp.phi += k2.phi * h * 0.5;
  temp.dt += k2.dt * h * 0.5;
  temp.dr += k2.dr * h * 0.5;
  temp.dtheta += k2.dtheta * h * 0.5;
  temp.dphi += k2.dphi * h * 0.5;
  derivatives(temp, k3);

  // k4
  temp = s;
  temp.t += k3.t * h;
  temp.r += k3.r * h;
  temp.theta += k3.theta * h;
  temp.phi += k3.phi * h;
  temp.dt += k3.dt * h;
  temp.dr += k3.dr * h;
  temp.dtheta += k3.dtheta * h;
  temp.dphi += k3.dphi * h;
  derivatives(temp, k4);

  // Update
  s.t += h * (k1.t + 2 * k2.t + 2 * k3.t + k4.t) / 6.0;
  s.r += h * (k1.r + 2 * k2.r + 2 * k3.r + k4.r) / 6.0;
  s.theta += h * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta) / 6.0;
  s.phi += h * (k1.phi + 2 * k2.phi + 2 * k3.phi + k4.phi) / 6.0;

  s.dt += h * (k1.dt + 2 * k2.dt + 2 * k3.dt + k4.dt) / 6.0;
  s.dr += h * (k1.dr + 2 * k2.dr + 2 * k3.dr + k4.dr) / 6.0;
  s.dtheta += h * (k1.dtheta + 2 * k2.dtheta + 2 * k3.dtheta + k4.dtheta) / 6.0;
  s.dphi += h * (k1.dphi + 2 * k2.dphi + 2 * k3.dphi + k4.dphi) / 6.0;
}
