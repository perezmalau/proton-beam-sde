#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <vector>

struct proton_path {

  proton_path(const double e0, const std::vector<double> x0,
              const std::vector<double> w0, const int n)
      : energy(n, e0), s(n, 0), x(n, x0), omega(n, w0), u(3, 0), z(3, 0),
        w(3, 0) {}

  void reset(const double e0, const std::vector<double> x0,
             const std::vector<double> w0) {
    x[0] = x0;
    omega[0] = w0;
    energy[0] = e0;
    s[0] = 0;
    return;
  }

  double log_a(const int k, const int m, const double theta) {
    double ret = log(theta + 2 * k - 1) + gsl_sf_lnpoch(theta + m, k - 1) -
                 gsl_sf_lnfact(m) - gsl_sf_lnfact(k - m);
    return ret;
  }

  double b(const int k, const int m, const double t, const double theta) {
    double ret = 1;
    if (k > 0) {
      ret = exp(log_a(k, m, theta) - k * (k + theta - 1) * t / 2);
    }
    return ret;
  }

  int c(const int m, const double t, const double theta) {
    int i = 0;
    double b_curr = b(m, m, t, theta);
    double b_next = b(m + 1, m, t, theta);
    while (b_next >= b_curr) {
      i++;
      b_curr = b_next;
      b_next = b(i + m + 1, m, t, theta);
    }
    return i;
  }

  int number_of_blocks(const double t, gsl_rng *gen) {
    int m = 0;
    double theta = 1;
    if (t < 0.07) {
      double mu = 2 / t;
      double sigma = sqrt(2 / (3 * t));
      m = round(mu + sigma * gsl_ran_gaussian_ziggurat(gen, 1));
    } else {
      std::vector<int> k(1, 0);
      bool proceed = true;
      double u = gsl_rng_uniform_pos(gen);
      double smin = 0, smax = 0, increment = 0;
      while (proceed) {
        k[m] = ceil(c(m, t, theta) / 2.0);
        for (int i = 0; i < k[m]; i++) {
          increment = b(m + 2 * i, m, t, theta) - b(m + 2 * i + 1, m, t, theta);
          smin += increment;
          smax += increment;
        }
        increment = b(m + 2 * k[m], m, t, theta);
        smin += increment - b(m + 2 * k[m] + 1, m, t, theta);
        smax += increment;
        while (smin < u && u < smax) {
          for (int i = 0; i <= m; i++) {
            k[i]++;
            increment = b(i + 2 * k[i], i, t, theta);
            smax = smin + increment;
            smin += increment - b(i + 2 * k[i] + 1, i, t, theta);
          }
        }
        if (smin > u) {
          proceed = false;
        } else {
          k.push_back(0);
          m++;
        }
      }
    }
    return m;
  }

  double wright_fisher(const double r, gsl_rng *gen) {
    double y;
    if (r > 1e-9) {
      int m = number_of_blocks(r, gen);
      y = gsl_ran_beta(gen, 1, 1 + m);
    } else {
      y = r / 2;
      y = fabs(gsl_ran_gaussian_ziggurat(gen, sqrt(r * y * (1 - y))));
    }
    return y;
  }

  double dist(const std::vector<double> &y, const std::vector<double> &z) {
    double ret = 0;
    for (unsigned int i = 0; i < y.size(); i++) {
      ret += (z[i] - y[i]) * (z[i] - y[i]);
    }
    return sqrt(ret);
  }

  double bethe_bloch(const double e) {
    double rho = 1;         // density of medium, g / cm^3
    double Z_oxy = 8;       // atomic number
    double Z_hyd = 1;       // atomic number
    double A_oxy = 16;      // atomic mass
    double A_hyd = 1;       // atomic mass
    double mecsq = 0.511;   // mass of electron * speed of light squared, MeV
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    // mean excitation energy of medium, MeV, fitted by eye
    double I = 60 * 1e-6;
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double b_oxy = 0.3072 * Z_oxy * rho *
                   (log(2 * mecsq * betasq / I * (1 - betasq)) - betasq) /
                   (A_oxy * betasq); // MeV / cm
    double b_hyd = 0.3072 * Z_hyd * rho *
                   (log(2 * mecsq * betasq / I * (1 - betasq)) - betasq) /
                   (A_hyd * betasq); // MeV / cm
    double ret = (A_oxy * b_oxy + 2 * A_hyd * b_hyd) / (A_oxy + 2 * A_hyd);
    return ret;
  }

  double multiple_scattering_sd(const double e, const double dt) {
    double Z_oxy = 8;       // atomic number
    double Z_hyd = 1;       // atomic number
    double A_oxy = 16;      // atomic mass
    double A_hyd = 1;       // atomic mass
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    // radiation length of oxygen, g / cm^2
    double x_oxy =
        A_oxy * 716.4 / (Z_oxy * (Z_oxy + 1) * log(287 / sqrt(Z_oxy)));
    // radiation length of hydrogen, g / cm^2
    double x_hyd =
        A_hyd * 716.4 / (Z_hyd * (Z_hyd + 1) * log(287 / sqrt(Z_hyd)));
    // radiation length of water via Bragg additivity rule
    double x0 = (2 * A_hyd + A_oxy) * x_oxy * x_hyd /
                (A_oxy * x_hyd + 2 * A_hyd * x_oxy);
    double pv = (2 * mpcsq + e) * e / (mpcsq + e);
    double ret = 14.1 * sqrt(dt / x0) * (1 + log10(dt / x0) / 9) / pv;
    return ret;
  }

  double energy_straggling_sd() {
    double alpha = 1 / 137.0;
    double log_hbar = -21 * log(10) + log(4.136) - log(2 * M_PI); // MeV * s
    double log_c = log(29979245800);                              // cm / s
    double log_water_density = 22 * log(10) + log(3.345); // molecules / cm^3
    double z = 10; // electrons per water molecule
    double ret = 4 * M_PI * z *
                 exp(2 * (log(alpha) + log_hbar + log_c) + log_water_density);
    return sqrt(ret);
  }

  void spherical_bm(const double time_increment, int &ix, gsl_rng *gen) {
    z[0] = sin(omega[ix - 1][0]) * cos(omega[ix - 1][1]);
    z[1] = sin(omega[ix - 1][0]) * sin(omega[ix - 1][1]);
    z[2] = cos(omega[ix - 1][0]);

    double y = wright_fisher(
        pow(multiple_scattering_sd(energy[ix - 1], time_increment), 2), gen);
    double theta = 2 * M_PI * gsl_rng_uniform(gen);
    // Set up defaults for when z is near (0, 0, 1)
    u[0] = 1;
    u[1] = 1;
    u[2] = 1;
    if (z[0] > 0) {
      u[0] = -1;
    }
    if (z[1] > 0) {
      u[1] = -1;
    }
    double denom = sqrt(z[0] * z[0] + z[1] * z[1] + (z[2] - 1) * (z[2] - 1));
    if (denom > 1e-10) {
      u[0] = -z[0] / denom;
      u[1] = -z[1] / denom;
      u[2] = (1 - z[2]) / denom;
    }
    z[0] = 2 * sqrt(y * (1 - y)) * cos(theta);
    z[1] = 2 * sqrt(y * (1 - y)) * sin(theta);
    z[2] = 1 - 2 * y;
    w[0] = (1 - 2 * u[0] * u[0]) * z[0] - 2 * u[0] * u[1] * z[1] -
           2 * u[0] * u[2] * z[2];
    w[1] = (1 - 2 * u[1] * u[1]) * z[1] - 2 * u[0] * u[1] * z[0] -
           2 * u[1] * u[2] * z[2];
    w[2] = (1 - 2 * u[2] * u[2]) * z[2] - 2 * u[0] * u[2] * z[0] +
           2 * u[1] * u[2] * z[1];
    if (ix == int(omega.size())) {
      omega.resize(2 * omega.size(), omega.back());
      energy.resize(2 * energy.size(), energy.back());
      x.resize(2 * x.size(), x.back());
      s.resize(2 * s.size(), s.back());
    }
    omega[ix][0] = acos(w[2]);
    omega[ix][1] = atan2(w[1], w[0]);
    double v0 = omega[ix - 1][0];
    double v1 = omega[ix][0];
    double w0 = omega[ix - 1][1];
    double w1 = omega[ix][1];
    double dt = time_increment;
    x[ix][0] = x[ix - 1][0] -
               dt *
                   ((v0 - v1) * (cos(v0) * cos(w0) - cos(v1) * cos(w1)) +
                    (w0 - w1) * (sin(v0) * sin(w0) - sin(v1) * sin(w1))) /
                   ((v0 - v1 + w0 - w1) * (v0 - v1 - w0 + w1));
    x[ix][1] = x[ix - 1][1] +
               dt *
                   ((w0 - w1) * (cos(w0) * sin(v0) - cos(w1) * sin(v1)) -
                    (v0 - v1) * (cos(v0) * sin(w0) - cos(v1) * sin(w1))) /
                   ((v0 - v1 + w0 - w1) * (v0 - v1 - w0 + w1));
    if (fabs(v0 - v1) > 1e-9) {
      x[ix][2] = x[ix - 1][2] + (sin(v0) - sin(v1)) * dt / (v0 - v1);
    } else {
      x[ix][2] = x[ix - 1][2] - (cos(v0) + cos(v1)) / 2;
    }
    energy[ix] = energy[ix - 1] -
                 bethe_bloch(energy[ix - 1]) * dist(x[ix - 1], x[ix]) +
                 sqrt(dist(x[ix - 1], x[ix])) * energy_straggling_sd() *
                     gsl_ran_gaussian_ziggurat(gen, 1);
    energy[ix] = fmax(energy[ix], 0);
    s[ix] = (energy[ix - 1] - energy[ix]);
    ix++;
    return;
  }

  double nonelastic_cross_section(const double e, const CS_1d &nonelastic_cs) {
    double ret = nonelastic_cs.evaluate(e);
    double log_water_density = 22 * log(10) + log(3.345); // molecules / cm^3
    double log_barns_to_cmsq = -24 * log(10);
    double A_oxy = 16;
    double A_hyd = 1;
    double c = A_oxy * exp(log_barns_to_cmsq + log_water_density) /
               (A_oxy + 2 * A_hyd);
    return c * ret; // rate per cm
  }

  double elastic_cross_section(const double e, const double lb) {
    double Z_oxy = 8;  // atomic number
    double Z_hyd = 1;  // atomic number
    double A_oxy = 16; // atomic mass
    double A_hyd = 1;  // atomic mass
    double log_ahbarc = log(197.3 / 137) - 13 * log(10);
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double pv = (2 * mpcsq + e) * e / (mpcsq + e);
    double log_density = 22 * log(10) + log(3.345); // molecules / cm^3
    double sig =
        exp(2 * (log_ahbarc + log(cos(lb / 2)) - log(pv) - log(sin(lb / 2))) +
            log_density) *
        M_PI;
    double ret = sig * (A_oxy * pow(Z_oxy, 2) + 2 * A_hyd * pow(Z_hyd, 2)) /
                 (A_oxy + 2 * A_hyd);
    return ret;
  }

  void elastic_scatter(std::vector<double> &ang, const double lb,
                       gsl_rng *gen) {
    double beta = 2 * M_PI * gsl_rng_uniform(gen);
    double u = gsl_rng_uniform(gen);
    double alpha = acos((cos(lb) - u * pow(cos(lb / 2), 2)) /
                        (1 - u * pow(cos(lb / 2), 2)));
    if (0 <= beta && beta < M_PI / 2) {
      ang[0] -= atan(sin(beta) * tan(alpha));
      ang[1] -= atan(cos(beta) * tan(alpha));
    } else if (M_PI / 2 <= beta && beta < M_PI) {
      ang[0] -= atan(sin(M_PI - beta) * tan(alpha));
      ang[1] += atan(cos(M_PI - beta) * tan(alpha));
    } else if (M_PI <= beta && beta < 3 * M_PI / 2) {
      ang[0] += atan(sin(3 * M_PI / 2 - beta) * tan(alpha));
      ang[1] += atan(cos(3 * M_PI / 2 - beta) * tan(alpha));
    } else {
      ang[0] += atan(sin(2 * M_PI - beta) * tan(alpha));
      ang[1] -= atan(cos(2 * M_PI - beta) * tan(alpha));
    }
    return;
  }

  double nonelastic_scatter(std::vector<double> &ang, const double e,
                            gsl_rng *gen, const CS_2d &angle_cdf) {
    double beta = 2 * M_PI * gsl_rng_uniform(gen);
    double alpha = angle_cdf.sample(e, gen);
    if (0 <= beta && beta < M_PI / 2) {
      ang[0] -= atan(sin(beta) * tan(alpha));
      ang[1] -= atan(cos(beta) * tan(alpha));
    } else if (M_PI / 2 <= beta && beta < M_PI) {
      ang[0] -= atan(sin(M_PI - beta) * tan(alpha));
      ang[1] += atan(cos(M_PI - beta) * tan(alpha));
    } else if (M_PI <= beta && beta < 3 * M_PI / 2) {
      ang[0] += atan(sin(3 * M_PI / 2 - beta) * tan(alpha));
      ang[1] += atan(cos(3 * M_PI / 2 - beta) * tan(alpha));
    } else {
      ang[0] += atan(sin(2 * M_PI - beta) * tan(alpha));
      ang[1] -= atan(cos(2 * M_PI - beta) * tan(alpha));
    }
    return alpha;
  }

  int simulate(const double dt, const double absorption_energy, gsl_rng *gen,
               const CS_1d &nonelastic_cs, const CS_2d &angle_cdf,
               const CS_3d &energy_cdf) {
    double nonelastic_jump_rate, elastic_jump_rate, alpha, elastic_min_scatter,
        exit_cos;
    int ix = 1;
    while (energy[ix - 1] > absorption_energy) {
      spherical_bm(dt, ix, gen);
      if (energy[ix - 1] > absorption_energy) {
        nonelastic_jump_rate =
            nonelastic_cross_section(energy[ix - 1], nonelastic_cs);
        elastic_min_scatter = 2.5 * multiple_scattering_sd(energy[ix - 1], dt);
        elastic_jump_rate = 0;
        if (elastic_min_scatter < M_PI) {
          elastic_jump_rate =
              elastic_cross_section(energy[ix - 1], elastic_min_scatter);
        }
        alpha = elastic_jump_rate + nonelastic_jump_rate;
        if (gsl_rng_uniform(gen) < 1 - exp(-alpha * dt)) {
          if (gsl_rng_uniform(gen) < elastic_jump_rate / alpha) {
            // Undo change in angle due to diffusion here because it is
            // accounted for in the Rutherford cross section in elastic_scatter.
            omega[ix - 1][0] = omega[ix - 2][0];
            omega[ix - 1][1] = omega[ix - 2][1];
            elastic_scatter(omega[ix - 1], elastic_min_scatter, gen);
          } else {
            exit_cos = nonelastic_scatter(omega[ix - 1], energy[ix - 1], gen,
                                          angle_cdf);
            energy[ix - 1] = energy_cdf.sample(energy[ix - 1], exit_cos, gen);
            s[ix - 1] = (energy[ix - 2] - energy[ix - 1]);
          }
        }
      }
    }
    return ix;
  }

  void print() {
    for (unsigned int i = 0; i < x.size(); i++) {
      std::cout << x[i][0] << " " << x[i][1] << " " << x[i][2] << " "
                << omega[i][0] << " " << omega[i][1] << " " << energy[i] << " "
                << s[i] << std::endl;
    }
    return;
  }

  std::vector<double> energy, s;
  std::vector<std::vector<double>> x, omega;
  // Dummy vectors for spherical BM
  std::vector<double> u, z, w;
};
