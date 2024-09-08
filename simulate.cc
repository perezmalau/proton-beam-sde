#include "grid_2d.cc"
#include "grid_3d.cc"
#include "proton_beam.cc"
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <libconfig.h++>
#include <string>
#include <vector>

// Member function of proton_beam copied here for access by ODE track length
// solver
double dose_deposition(const double e) {
  double p = 1.77;
  double a = 2.2 * 1e-2; // in mm / MeV
  double ret = pow(e, 1 - p) / (a * p);
  return ret;
}

int solve_track_length(const double e0, const double dt,
                       const double absorption_e) {
  int ret = 0;
  double e = e0;
  while (e > absorption_e) {
    ret++;
    e -= dt * dose_deposition(e);
  }
  return ret;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cout << "Call " << argv[0] << " config-path" << std::endl;
    return 1;
  }
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));

  double nozzle_radius, e0;
  std::vector<double> x(3, 0), w(2, 0);
  w[0] = M_PI / 2;

  libconfig::Config cfg;
  cfg.readFile(argv[1]);
  cfg.lookupValue("nozzle_radius", nozzle_radius);
  double theta = 2 * M_PI * gsl_rng_uniform(gen);
  double r = nozzle_radius * gsl_rng_uniform(gen);
  x[1] = sqrt(r * nozzle_radius) * cos(theta);
  x[2] = sqrt(r * nozzle_radius) * sin(theta);

  double initial_e_min, initial_e_max;
  cfg.lookupValue("initial_e_min", initial_e_min);
  cfg.lookupValue("initial_e_max", initial_e_max);
  e0 = gsl_ran_flat(gen, initial_e_min, initial_e_max);

  double dt, scatter_angle_lb, absorption_e;
  int nrep;
  cfg.lookupValue("max_step_size", dt);
  cfg.lookupValue("min_elastic_scatter", scatter_angle_lb);
  cfg.lookupValue("absorption_energy", absorption_e);
  cfg.lookupValue("replicates", nrep);

  double grid_dx;
  cfg.lookupValue("grid_dx", grid_dx);

  int output_1d, output_2d, output_3d, output_2d_slice;
  std::string path_1d, path_2d, path_3d, path_slice;
  cfg.lookupValue("output_1d", output_1d);
  cfg.lookupValue("output_2d", output_2d);
  cfg.lookupValue("output_3d", output_3d);
  cfg.lookupValue("output_2d_slice", output_2d_slice);
  if (output_1d == 1) {
    cfg.lookupValue("out_path_1d", path_1d);
  }
  if (output_2d == 1) {
    cfg.lookupValue("out_path_2d", path_2d);
  }
  if (output_3d == 1) {
    cfg.lookupValue("out_path_3d", path_3d);
  }
  if (output_2d_slice == 1) {
    cfg.lookupValue("out_path_2d_slice", path_slice);
  }

  int n = solve_track_length(initial_e_max, dt, absorption_e);
  // Extra 10 slots to be safe
  proton_path p(e0, x, w, n + 10);
  // We need all grids to be initialised for code to compile
  // but make unneeded ones tiny
  int tmp = n + 10;
  if (output_1d == 0 && output_2d == 0) {
    tmp = 1;
  }
  grid_2d g2d(tmp, grid_dx);
  tmp = n + 10;
  if (output_3d == 0) {
    tmp = 1;
  }
  grid_3d g3d(tmp, grid_dx);
  tmp = n + 10;
  if (output_2d_slice == 0) {
    tmp = 1;
  }
  grid_2d g2d_slice(tmp, grid_dx);

  int len;
  for (int i = 0; i < nrep; i++) {
    len = p.simulate(dt, scatter_angle_lb, absorption_e, gen);
    if (output_1d == 1 || output_2d == 1) {
      g2d.add(p.x, p.s, len);
    }
    if (output_3d == 1) {
      g3d.add(p.x, p.s, len);
    }
    if (output_2d_slice == 1) {
      g2d_slice.add_to_slice(p.x, p.s, len);
    }
    e0 = gsl_ran_flat(gen, initial_e_min, initial_e_max);
    theta = 2 * M_PI * gsl_rng_uniform(gen);
    r = nozzle_radius * gsl_rng_uniform(gen);
    x[1] = sqrt(r * nozzle_radius) * cos(theta);
    x[2] = sqrt(r * nozzle_radius) * sin(theta);
    p.reset(e0, x, w);
  }
  if (output_1d == 1) {
    g2d.print_1d(path_1d);
  }
  if (output_2d == 1) {
    g2d.print(path_2d);
  }
  if (output_3d == 1) {
    g3d.print_2d(path_3d);
  }
  if (output_2d_slice == 1) {
    g2d_slice.print(path_slice);
  }
  gsl_rng_free(gen);
  return 1;
}
