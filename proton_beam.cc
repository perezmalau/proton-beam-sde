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
      : energy(n, e0), s(n, 0), x(n, x0), omega(n, w0), u(3, 0), z(3, 0), w(3, 0) {
  }

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

  double dose_deposition(const double e) {
    double p = 1.77;
    double a = 2.2 * 1e-2; // in mm / MeV
    double ret = pow(e, 1 - p) / (a * p);
    return ret;
  }

  double diffusion_coeff(const double e) {
    double ret = sigma_elastic(e);
    return ret;
  }

  void spherical_bm(const double time_increment, const double absorption_energy,
                    int &ix, gsl_rng *gen) {
    z[0] = sin(omega[ix - 1][0]) * cos(omega[ix - 1][1]);
    z[1] = sin(omega[ix - 1][0]) * sin(omega[ix - 1][1]);
    z[2] = cos(omega[ix - 1][0]);

    double y = wright_fisher(time_increment * diffusion_coeff(energy[ix - 1]), gen);
    double theta = 2 * M_PI * gsl_rng_uniform(gen);
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
    omega[ix][0] = acos(w[2]);
    omega[ix][1] = atan2(w[1], w[0]);
    x[ix][0] = x[ix - 1][0] + sin(omega[ix][0]) * cos(omega[ix][1]) * time_increment;
    x[ix][1] = x[ix - 1][1] + sin(omega[ix][0]) * sin(omega[ix][1]) * time_increment;
    x[ix][2] = x[ix - 1][2] + cos(omega[ix][0]) * time_increment;
    energy[ix] = fmax(absorption_energy, energy[ix - 1] - time_increment * dose_deposition(energy[ix - 1]));
    s[ix] = (energy[ix - 1] - energy[ix]) / dist(x[ix - 1], x[ix]);
    ix++;
    return;
  }

  double sigma_inelastic_ub() {
    double e_peak = 20; 
    double ret = sigma_inelastic(e_peak);
    return ret;
  }

  double sigma_inelastic(const double e) {
    double ret = 0;
    double e_min = 6;
    double e_peak = 20;
    double e_plateau = 100;
    double plateau = 0.3;
    double peak = 0.55;
    double grad;
    if (e_min < e && e <= e_peak) {
      grad = peak / (e_peak - e_min);
      ret = grad * (e - e_min);
    } else if (e_peak < e && e <= e_plateau) {
      grad = (plateau - peak) / (e_plateau - e_peak);
      ret = grad * (e - e_peak) + peak;
    } else {
      ret = plateau;
    }
    double c = 1e-3;
    return c * ret;
  }

  double sigma_elastic(const double e) {
    double c = 1e-2;
    double ret = c / (e * e);
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

  void inelastic_scatter(std::vector<double> &ang, gsl_rng *gen) {
    double alpha = 2;
    double beta = 2;
    if (0 < ang[0] && ang[0] <= M_PI / 2) {
      beta = M_PI / ang[0];
    } else {
      alpha = M_PI / (M_PI - ang[0]);
    }
    ang[0] = M_PI * gsl_ran_beta(gen, alpha, beta);

    if (ang[1] < 0) {
      ang[1] += 2 * M_PI;
    }
    if (0 < ang[1] && ang[1] <= M_PI) {
      alpha = 2;
      beta = 2 * M_PI / ang[1];
    } else {
      alpha = 2 * M_PI / (2 * M_PI - ang[1]);
      beta = 2;
    }
    ang[1] = 2 * M_PI * gsl_ran_beta(gen, alpha, beta);
    return;
  }

  int simulate(const double dt, const double lb, const double absorption_e, gsl_rng *gen) {
    double inelastic_jump_rate, elastic_jump_rate, e;
    double jump_rate, jump_time, alpha;
    int ix = 1;
    do {
      inelastic_jump_rate = sigma_inelastic_ub();
      elastic_jump_rate = sigma_elastic(energy[ix - 1]);
      jump_rate = inelastic_jump_rate + elastic_jump_rate;
      jump_time = gsl_ran_exponential(gen, 1 / jump_rate);

      int nsteps = ceil(jump_time / dt);
      int m = 0;
      while (m < nsteps && energy[ix - 1] > absorption_e) {
        spherical_bm(jump_time / nsteps, absorption_e, ix, gen);
        m++;
      }
      e = energy[ix - 1];
      alpha = exp(log(sigma_inelastic(e) + sigma_elastic(e)) - log(jump_rate));
      if (e > absorption_e && gsl_rng_uniform(gen) < alpha) {
        alpha = exp(log(sigma_elastic(e)) - log(sigma_inelastic(e) + sigma_elastic(e)));
        if (gsl_rng_uniform(gen) < alpha) {
          elastic_scatter(omega[ix - 1], lb, gen);
        } else {
          inelastic_scatter(omega[ix - 1], gen);
          energy[ix - 1] *= gsl_rng_uniform(gen);
          s[ix - 1] = (energy[ix - 2] - energy[ix - 1]) / dist(x[ix - 1], x[ix - 2]);
        }
      }
    } while (energy[ix - 1] > absorption_e);
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
