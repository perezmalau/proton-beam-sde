#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct CS_3d {

  CS_3d(const std::string filename) : energy(), angle(), exit_energy(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    getline(file, line);
    std::stringstream iss2;
    iss2 << line;
    while (getline(iss2, token, ' ')) {
      angle.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      cdf.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss4;
      iss4 << line;
      while (getline(iss4, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      exit_energy.push_back(tmp_vec);
    }
  }

  void print() const {
    for (unsigned int i = 0; i < energy.size(); i++) {
      std::cout << energy[i] << " ";
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < angle.size(); i++) {
      std::cout << angle[i] << " ";
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < cdf.size(); i++) {
      for (unsigned int j = 0; j < cdf[i].size(); j++) {
        std::cout << cdf[i][j] << " ";
      }
      std::cout << std::endl;
      for (unsigned int j = 0; j < exit_energy[i].size(); j++) {
        std::cout << exit_energy[i][j] << " ";
      }
      std::cout << std::endl;
    }
    return;
  }

  int find_energy_index(const double e) const {
    int ret = 0;
    while (ret < int(energy.size()) && e > energy[ret]) {
      ret++;
    }
    return ret;
  }

  int find_angle_index(const double ang) const {
    int ret = 0;
    while (ret < int(angle.size()) && ang > angle[ret]) {
      ret++;
    }
    return ret;
  }

  double sample(const double e, const double ang, gsl_rng *gen) const {
    double u = gsl_rng_uniform(gen);
    int e_row = find_energy_index(e);
    int ang_row = find_angle_index(ang);
    double ret, ret2, ret3, ret4;
    // CDF column counter initialised at 1 because all CDFs begin with a zero
    // entry
    int c = 1;
    int r;
    if (e_row == 0) {
      if (ang_row == 0) {
        while (u > cdf[0][c]) {
          c++;
        }
        ret = ((cdf[0][c] - u) * fmin(e, exit_energy[0][c - 1]) +
               (u - cdf[0][c - 1]) * fmin(e, exit_energy[0][c])) /
              (cdf[0][c] - cdf[0][c - 1]);
      } else if (ang_row == int(angle.size())) {
        r = ang_row - 1;
        while (u > cdf[r][c]) {
          c++;
        }
        ret = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
               (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
              (cdf[r][c] - cdf[r][c - 1]);
      } else {
        r = ang_row;
        while (u > cdf[r][c]) {
          c++;
        }
        ret = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
               (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
              (cdf[r][c] - cdf[r][c - 1]);
        c = 1;
        while (u > cdf[r - 1][c]) {
          c++;
        }
        ret2 = ((cdf[r - 1][c] - u) * fmin(e, exit_energy[r - 1][c - 1]) +
                (u - cdf[r - 1][c - 1]) * fmin(e, exit_energy[r - 1][c])) /
               (cdf[r - 1][c] - cdf[r - 1][c - 1]);
        ret = ((angle[r] - ang) * ret2 + (ang - angle[r - 1]) * ret) /
              (angle[r] - angle[r - 1]);
      }
    } else if (e_row == int(energy.size())) {
      if (ang_row == 0) {
        r = (e_row - 1) * angle.size();
        while (u > cdf[r][c]) {
          c++;
        }
        ret = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
               (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
              (cdf[r][c] - cdf[r][c - 1]);
      } else if (ang_row == int(angle.size())) {
        r = (e_row - 1) * angle.size() + ang_row - 1;
        while (u > cdf[r][c]) {
          c++;
        }
        ret = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
               (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
              (cdf[r][c] - cdf[r][c - 1]);
      } else {
        r = (e_row - 1) * angle.size() + ang_row;
        while (u > cdf[r][c]) {
          c++;
        }
        ret = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
               (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
              (cdf[r][c] - cdf[r][c - 1]);
        c = 1;
        while (u > cdf[r - 1][c]) {
          c++;
        }
        ret2 = ((cdf[r - 1][c] - u) * fmin(e, exit_energy[r - 1][c - 1]) +
                (u - cdf[r - 1][c - 1]) * fmin(e, exit_energy[r - 1][c])) /
               (cdf[r - 1][c] - cdf[r - 1][c - 1]);
        ret =
            ((angle[ang_row] - ang) * ret2 + (ang - angle[ang_row - 1]) * ret) /
            (angle[ang_row] - angle[ang_row - 1]);
      }
    } else {
      if (ang_row == 0) {
        r = e_row * angle.size();
        while (u > cdf[r][c]) {
          c++;
        }
        ret = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
               (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
              (cdf[r][c] - cdf[r][c - 1]);
        c = 1;
        r = (e_row - 1) * angle.size();
        while (u > cdf[r][c]) {
          c++;
        }
        ret2 = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
                (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
               (cdf[r][c] - cdf[r][c - 1]);
        ret = ((energy[e_row] - e) * ret2 + (e - energy[e_row - 1]) * ret) /
              (energy[e_row] - energy[e_row - 1]);
      } else if (ang_row == int(angle.size())) {
        r = e_row * angle.size() + ang_row - 1;
        while (u > cdf[r][c]) {
          c++;
        }
        ret = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
               (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
              (cdf[r][c] - cdf[r][c - 1]);
        c = 1;
        r = (e_row - 1) * angle.size() + ang_row - 1;
        while (u > cdf[r][c]) {
          c++;
        }
        ret2 = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
                (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
               (cdf[r][c] - cdf[r][c - 1]);
        ret = ((energy[e_row] - e) * ret2 + (e - energy[e_row - 1]) * ret) /
              (energy[e_row] - energy[e_row - 1]);
      } else {
        r = (e_row - 1) * angle.size() + ang_row - 1;
        while (u > cdf[r][c]) {
          c++;
        }
        ret = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
               (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
              (cdf[r][c] - cdf[r][c - 1]);
        r = (e_row - 1) * angle.size() + ang_row;
        c = 1;
        while (u > cdf[r][c]) {
          c++;
        }
        ret2 = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
                (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
               (cdf[r][c] - cdf[r][c - 1]);
        ret =
            ((angle[ang_row] - ang) * ret + (ang - angle[ang_row - 1]) * ret2) /
            (angle[ang_row] - angle[ang_row - 1]);
        r = e_row * angle.size() + ang_row - 1;
        c = 1;
        while (u > cdf[r][c]) {
          c++;
        }
        ret3 = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
                (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
               (cdf[r][c] - cdf[r][c - 1]);
        r = e_row * angle.size() + ang_row;
        c = 1;
        while (u > cdf[r][c]) {
          c++;
        }
        ret4 = ((cdf[r][c] - u) * fmin(e, exit_energy[r][c - 1]) +
                (u - cdf[r][c - 1]) * fmin(e, exit_energy[r][c])) /
               (cdf[r][c] - cdf[r][c - 1]);
        ret3 = ((angle[ang_row] - ang) * ret3 +
                (ang - angle[ang_row - 1]) * ret4) /
               (angle[ang_row] - angle[ang_row - 1]);
        ret = ((energy[e_row] - e) * ret + (e - energy[e_row - 1]) * ret3) /
              (energy[e_row] - energy[e_row - 1]);
      }
    }
    return ret;
  }

  std::vector<double> energy, angle;
  std::vector<std::vector<double>> exit_energy, cdf;
};

struct CS_2d {

  CS_2d(const std::string filename) : energy(), angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    getline(file, line);
    std::stringstream iss2;
    iss2 << line;
    while (getline(iss2, token, ' ')) {
      angle.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      cdf.push_back(tmp_vec);
    }
  }

  void print() const {
    for (unsigned int i = 0; i < energy.size(); i++) {
      std::cout << energy[i] << " ";
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < angle.size(); i++) {
      std::cout << angle[i] << " ";
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < cdf.size(); i++) {
      for (unsigned int j = 0; j < cdf[i].size(); j++) {
        std::cout << cdf[i][j] << " ";
      }
      std::cout << std::endl;
    }
    return;
  }

  double sample(const double e, gsl_rng *gen) const {
    double u = gsl_rng_uniform(gen);
    int r, c = 1;
    double ret, ret2;
    if (e <= energy[0]) {
      while (u > cdf[0][c]) {
        c++;
      }
      ret = ((cdf[0][c] - u) * angle[c - 1] + (u - cdf[0][c - 1]) * angle[c]) /
            (cdf[0][c] - cdf[0][c - 1]);
    } else if (e >= energy.back()) {
      r = energy.size() - 1;
      while (u > cdf[r][c]) {
        c++;
      }
      ret = ((cdf[r][c] - u) * angle[c - 1] + (u - cdf[r][c - 1]) * angle[c]) /
            (cdf[r][c] - cdf[r][c - 1]);
    } else {
      r = 1;
      while (e > energy[r]) {
        r++;
      }
      while (u > cdf[r][c]) {
        c++;
      }
      ret = ((cdf[r][c] - u) * angle[c - 1] + (u - cdf[r][c - 1]) * angle[c]) /
            (cdf[r][c] - cdf[r][c - 1]);
      c = 1;
      while (u > cdf[r - 1][c]) {
        c++;
      }
      ret2 = ((cdf[r - 1][c] - u) * angle[c - 1] +
              (u - cdf[r - 1][c - 1]) * angle[c]) /
             (cdf[r - 1][c] - cdf[r - 1][c - 1]);
      ret = ((energy[r] - e) * ret2 + (e - energy[r - 1]) * ret) /
            (energy[r] - energy[r - 1]);
    }
    return ret;
  }

  std::vector<double> energy, angle;
  std::vector<std::vector<double>> cdf;
};

struct CS_1d {

  CS_1d(const std::string filename) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    std::stringstream iss;

    getline(file, line);
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }

    getline(file, line);
    std::stringstream iss2;
    iss2 << line;
    while (getline(iss2, token, ' ')) {
      rate.push_back(atof(token.c_str()));
    }
  }

  void print() const {
    for (unsigned int i = 0; i < energy.size(); i++) {
      std::cout << energy[i] << " ";
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < rate.size(); i++) {
      std::cout << rate[i] << " ";
    }
    std::cout << std::endl;
    return;
  }

  double evaluate(const double e) const {
    int r;
    double ret;
    if (e <= energy[0]) {
      ret = rate[0];
    } else if (e >= energy.back()) {
      ret = rate.back();
    } else {
      r = 1;
      while (e > energy[r]) {
        r++;
      }
      ret = ((energy[r] - e) * rate[r - 1] + (e - energy[r - 1]) * rate[r]) /
            (energy[r] - energy[r - 1]);
    }
    return ret;
  }

  std::vector<double> energy, rate;
};
