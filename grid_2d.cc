#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct grid_2d {

  grid_2d(const int n, const double dx_) : dx(dx_), x_pos(n), x_neg(n) {
    std::vector<double> tmp(1, 0);
    for (int i = 0; i < n; i++) {
      x_pos[i] = tmp;
      x_neg[i] = tmp;
    }
  }

  void print_1d(const std::string filename) {
    std::ofstream outfile(filename);
    double val;
    for (unsigned int j = 0; j < x_pos.size(); j++) {
      val = 0;
      for (unsigned int k = 0; k < x_pos[j].size(); k++) {
        val += x_pos[j][k];
      }
      for (unsigned int k = 0; k < x_neg[j].size(); k++) {
        val += x_neg[j][k];
      }
      if (val > 0) {
        outfile << dx * j << " " << val << std::endl;
      }
    }
    outfile.close();
    return;
  }

  void print(const std::string filename) {
    std::ofstream outfile(filename);
    for (unsigned int j = 0; j < x_pos.size(); j++) {
      for (unsigned int k = 0; k < x_pos[j].size(); k++) {
        if (x_pos[j][k] > 0) {
          outfile << dx * j << " " << dx * k << " " << x_pos[j][k] << std::endl;
        }
      }
      for (unsigned int k = 0; k < x_neg[j].size(); k++) {
        if (x_neg[j][k] > 0) {
          outfile << dx * j << " " << -dx * (k + 1) << " " << x_neg[j][k]
                  << std::endl;
        }
      }
    }
    outfile.close();
    return;
  }

  void add_to_slice(const std::vector<std::vector<double>> &y,
                    const std::vector<double> &s, const int len) {
    int ix;
    unsigned int iy;
    for (int i = 0; i < len; i++) {
      ix = floor(y[i][0] / dx);
      iy = floor(fabs(y[i][1]) / dx);
      if (ix > 0 && y[i][2] >= 0 && y[i][2] < dx) {
        if (y[i][1] >= 0) {
          if (iy >= x_pos[ix].size()) {
            x_pos[ix].resize(iy + 1, 0);
          }
          x_pos[ix][iy] += s[i];
        } else {
          if (iy >= x_neg[ix].size()) {
            x_neg[ix].resize(iy + 1, 0);
          }
          x_neg[ix][iy] += s[i];
        }
      }
    }
    return;
  }

  void add(const std::vector<std::vector<double>> &y,
           const std::vector<double> &s, const int len) {
    int ix;
    unsigned int iy;
    for (int i = 0; i < len; i++) {
      ix = floor(y[i][0] / dx);
      iy = floor(fabs(y[i][1]) / dx);
      if (ix > 0) {
        if (y[i][1] >= 0) {
          if (iy >= x_pos[ix].size()) {
            x_pos[ix].resize(iy + 1, 0);
          }
          x_pos[ix][iy] += s[i];
        } else {
          if (iy >= x_neg[ix].size()) {
            x_neg[ix].resize(iy + 1, 0);
          }
          x_neg[ix][iy] += s[i];
        }
      }
    }
    return;
  }

  double dx;
  std::vector<std::vector<double>> x_pos, x_neg;
};
