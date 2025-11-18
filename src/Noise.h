#ifndef NOISE_H
#define NOISE_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

class Noise {
  std::vector<int> p;
  std::vector<double> r;

public:
  Noise(unsigned int seed = 2023) {
    p.resize(512);
    r.resize(256);
    std::iota(p.begin(), p.begin() + 256, 0);
    std::mt19937 gen(seed);
    std::shuffle(p.begin(), p.begin() + 256, gen);
    for (int i = 0; i < 256; ++i)
      p[256 + i] = p[i];

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < 256; ++i)
      r[i] = dist(gen);
  }

  double noise(double x, double y, double z) const {
    int X = (int)floor(x) & 255;
    int Y = (int)floor(y) & 255;
    int Z = (int)floor(z) & 255;

    x -= floor(x);
    y -= floor(y);
    z -= floor(z);

    double u = fade(x);
    double v = fade(y);
    double w = fade(z);

    int A = p[X] + Y, AA = p[A] + Z, AB = p[A + 1] + Z;
    int B = p[X + 1] + Y, BA = p[B] + Z, BB = p[B + 1] + Z;

    return lerp(
        w,
        lerp(v, lerp(u, grad(p[AA], x, y, z), grad(p[BA], x - 1, y, z)),
             lerp(u, grad(p[AB], x, y - 1, z), grad(p[BB], x - 1, y - 1, z))),
        lerp(v,
             lerp(u, grad(p[AA + 1], x, y, z - 1),
                  grad(p[BA + 1], x - 1, y, z - 1)),
             lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                  grad(p[BB + 1], x - 1, y - 1, z - 1))));
  }

  double turbulence(double x, double y, double z, int octaves = 4) const {
    double t = 0;
    double scale = 1.0;
    for (int i = 0; i < octaves; ++i) {
      t += noise(x * scale, y * scale, z * scale) / scale;
      scale *= 2.0;
    }
    return t;
  }

private:
  double fade(double t) const { return t * t * t * (t * (t * 6 - 15) + 10); }
  double lerp(double t, double a, double b) const { return a + t * (b - a); }
  double grad(int hash, double x, double y, double z) const {
    int h = hash & 15;
    double u = h < 8 ? x : y;
    double v = h < 4 ? y : h == 12 || h == 14 ? x : z;
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
  }
};

#endif
