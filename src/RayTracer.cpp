#include "RayTracer.h"
#include "ImageWriter.h"
#include "Noise.h"
#include <cmath>
#include <iostream>

Noise noise_gen;

RayTracer::RayTracer(int width, int height) : width(width), height(height) {}

void RayTracer::render(const std::string &filename) {
  std::vector<unsigned char> image(width * height * 3);

  Vec3 cam_pos(0, 1.0, -15);
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
      double x = (2 * (i + 0.5) / (double)width - 1) * aspect_ratio * scale;
      double y = (1 - 2 * (j + 0.5) / (double)height) * scale;

      Vec3 dir = unit_vector(forward + x * right + y * up);
      State s = initial_state(cam_pos, dir);

      bool hit_horizon = false;
      Vec3 accumulated_color(0, 0, 0);

      for (int step = 0; step < 2000;
           ++step) {
        integrate(s, step_size, accumulated_color);

        if (s.r < horizon * 1.01) {
          hit_horizon = true;
          break;
        }
        if (s.r > max_dist) {
          break;
        }
      }

      Vec3 color = accumulated_color;
      if (!hit_horizon) {
        color += get_color(s.theta, s.phi);
      }

      Vec3 white(1, 1, 1);
      Vec3 denom = color + white;
      color = Vec3(color.x / denom.x, color.y / denom.y, color.z / denom.z);
      color = Vec3(pow(color.x, 1 / 2.2), pow(color.y, 1 / 2.2),
                   pow(color.z, 1 / 2.2));

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
  double u = phi / (2 * M_PI);
  double v = theta / M_PI;

  double galaxy = std::exp(-10.0 * (v - 0.5) * (v - 0.5));

  bool grid = (int(u * 20) % 2 == 0) ^ (int(v * 10) % 2 == 0);

  Vec3 base = grid ? Vec3(0.1, 0.1, 0.2) : Vec3(0.05, 0.05, 0.1);

  base += Vec3(0.8, 0.7, 0.5) * galaxy;

  return base;
}

State RayTracer::initial_state(const Vec3 &origin, const Vec3 &direction) {
  State s;

  double r = origin.length();
  double theta = acos(origin.z / r);
  double phi = atan2(origin.y, origin.x);

  s.t = 0;
  s.r = r;
  s.theta = theta;
  s.phi = phi;

  double st = sin(theta);
  double ct = cos(theta);
  double sp = sin(phi);
  double cp = cos(phi);

  Vec3 e_r(st * cp, st * sp, ct);
  Vec3 e_theta(ct * cp, ct * sp, -st);
  Vec3 e_phi(-sp, cp, 0);

  s.dr = dot(direction, e_r);
  s.dtheta = dot(direction, e_theta) / r;
  s.dphi = dot(direction, e_phi) / (r * st);

  double metric_rr = 1.0 / (1.0 - 2.0 * M / r);
  double metric_ang = r * r * (s.dtheta * s.dtheta + st * st * s.dphi * s.dphi);
  double rhs = metric_rr * s.dr * s.dr + metric_ang;

  double factor = (1.0 - 2.0 * M / r);
  s.dt = sqrt(rhs / factor) / sqrt(factor);

  s.dt = sqrt((s.dr * s.dr) / pow(factor, 2) +
              (r * r * (s.dtheta * s.dtheta + st * st * s.dphi * s.dphi)) /
                  factor);

  return s;
}

void RayTracer::derivatives(const State &s, State &out) {
  double r = s.r;
  double theta = s.theta;

  double r_2 = r * r;
  double r_3 = r_2 * r;
  double st = sin(theta);
  double ct = cos(theta);
  double st_2 = st * st;

  double G_t_tr = M / (r * (r - 2.0 * M));
  double G_r_tt = M * (r - 2.0 * M) / r_3;
  double G_r_rr = -M / (r * (r - 2.0 * M));
  double G_r_thth = -(r - 2.0 * M);
  double G_r_phph = -(r - 2.0 * M) * st_2;
  double G_th_rth = 1.0 / r;
  double G_th_phph = -st * ct;
  double G_ph_rph = 1.0 / r;
  double G_ph_thph = ct / st;

  double ddt = -2.0 * G_t_tr * s.dt * s.dr;

  double ddr = -(G_r_tt * s.dt * s.dt + G_r_rr * s.dr * s.dr +
                 G_r_thth * s.dtheta * s.dtheta + G_r_phph * s.dphi * s.dphi);

  double ddtheta =
      -(2.0 * G_th_rth * s.dr * s.dtheta + G_th_phph * s.dphi * s.dphi);

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

Vec3 RayTracer::get_disk_color(double r, double phi) {
  double temp = 1.0 / (r - 1.5 * M);

  double n = noise_gen.turbulence(r * 2.0, phi * 5.0, 0.0, 4);
  double intensity = temp * (0.5 + 0.5 * n);

  Vec3 color(1.0, 0.8, 0.6);
  color *= intensity * 0.5;

  return color;
}

void RayTracer::integrate(State &s, double h, Vec3 &accumulated_color) {
  State k1, k2, k3, k4;
  State temp;

  double old_theta = s.theta;

  derivatives(s, k1);

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

  s.t += h * (k1.t + 2 * k2.t + 2 * k3.t + k4.t) / 6.0;
  s.r += h * (k1.r + 2 * k2.r + 2 * k3.r + k4.r) / 6.0;
  s.theta += h * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta) / 6.0;
  s.phi += h * (k1.phi + 2 * k2.phi + 2 * k3.phi + k4.phi) / 6.0;

  s.dt += h * (k1.dt + 2 * k2.dt + 2 * k3.dt + k4.dt) / 6.0;
  s.dr += h * (k1.dr + 2 * k2.dr + 2 * k3.dr + k4.dr) / 6.0;
  s.dtheta += h * (k1.dtheta + 2 * k2.dtheta + 2 * k3.dtheta + k4.dtheta) / 6.0;
  s.dphi += h * (k1.dphi + 2 * k2.dphi + 2 * k3.dphi + k4.dphi) / 6.0;

  double new_theta = s.theta;

  if ((old_theta - M_PI_2) * (new_theta - M_PI_2) <= 0) {
    double f = std::abs(old_theta - M_PI_2) /
               (std::abs(old_theta - M_PI_2) + std::abs(new_theta - M_PI_2));
    double r_cross = s.r * f + (s.r - k1.r * h) * (1 - f);
    r_cross = s.r;

    if (r_cross > disk_inner && r_cross < disk_outer) {
      Vec3 disk_col = get_disk_color(r_cross, s.phi);

      double Omega = sqrt(M / pow(r_cross, 3));
      double u_t = 1.0 / sqrt(1.0 - 3.0 * M / r_cross);
      double u_phi = Omega * u_t;

      double g_tt = -(1.0 - 2.0 * M / r_cross);
      double g_phph = r_cross * r_cross;

      double E_emit = -(g_tt * u_t * s.dt + g_phph * u_phi * s.dphi);
      double E_obs = 1.0;

      double redshift = E_emit / E_obs;
      double doppler = 1.0 / redshift;

      double beaming = pow(doppler, 4.0);

      beaming = std::min(beaming, 20.0);

      accumulated_color += disk_col * beaming * 0.3;
    }
  }
}
