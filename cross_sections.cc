#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef CS
#define CS
// necessary coef for proton-proton is ia = oa = 1 iz = oz = 0,
// this class is incident/outgoing particle dependent generally,
// preparing setup for general case
// defining here instead of materials as
// appears in ENDF cross section struct as input to member function
struct AtomConsts {
  const int a, z;
  const double Mi, Mo, ANU, BNU, iNU, oNU, Si, So;
  double AngleScatterCoef(const int AC, const int NC, const int ZC,
                          const int AX, const int NX, const int ZX,
                          const double sep) const {
    double temp = 15.68 * (AC - AX) -
                  (28.07 * ((pow(NC - ZC, 2) / AC) - (pow(NX - ZX, 2) / AX))) -
                  (18.56 * (pow(AC, 2.0 / 3.0) - pow(AX, 2.0 / 3.0))) +
                  (33.22 * ((pow(NC - ZC, 2) / pow(AC, 4.0 / 3.0)) -
                            (pow(NX - ZX, 2) / pow(AX, 4.0 / 3.0)))) -
                  (0.717 * ((pow(ZC, 2) / pow(AC, 1.0 / 3.0)) -
                            (pow(ZX, 2) / pow(AX, 1.0 / 3.0)))) +
                  (1.211 * ((pow(ZC, 2) / AC) - (pow(ZX, 2) / AX))) - sep;
    return temp;
  }
  AtomConsts(const int a0, const int z0, const int ia, const int iz,
             const int oa, const int oz, const int mi, const int mo,
             const double sepi, const double sepo)
      : a(a0), z(z0), Mi(mi), Mo(mo),
        ANU(((a0 - z0) * 1.0087 + z0 * 1.0073) / 1.0087),
        BNU(((a0 + ia - oa - z0 - iz + oz) * 1.0087 + (z0 + iz - oz) * 1.0073) /
            1.0087),
        iNU(((ia - iz) * 1.0087 + iz * 1.0073) / 1.0087),
        oNU(((oa - oz) * 1.0087 + oz * 1.0073) / 1.0087),
        Si(AngleScatterCoef(a0 + ia, a0 + ia - z0 - iz, z0 + iz, a0, a0 - z0,
                            z0, sepi)),
        So(AngleScatterCoef(a0 + oa, a0 + oa - z0 - oz, z0 + oz, a0, a0 - z0,
                            z0, sepo)) {}

  AtomConsts(const AtomConsts &other)
      : a(other.a), z(other.z), Mi(other.Mi), Mo(other.Mo), ANU(other.ANU),
        BNU(other.BNU), iNU(other.iNU), oNU(other.oNU), Si(other.Si),
        So(other.So) {}
  AtomConsts() : a(), z(), Mi(), Mo(), ANU(), BNU(), iNU(), oNU(), Si(), So() {}
};
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

  CS_3d(const CS_3d &other)
      : energy(other.energy), angle(other.angle),
        exit_energy(other.exit_energy), cdf(other.cdf) {}

  CS_3d() : energy(), angle(), exit_energy(), cdf() {}

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

  CS_2d(const CS_2d &other)
      : energy(other.energy), angle(other.angle), cdf(other.cdf) {}

  CS_2d() : energy(), angle(), cdf() {}

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

  CS_1d(const CS_1d &other) : energy(other.energy), rate(other.rate) {}

  CS_1d() : energy(), rate() {}

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
    double ret = 0;
    if (energy.size() > 0) {
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
    }
    return ret;
  }

  std::vector<double> energy, rate;
};

struct CS_3dENDF {

  CS_3dENDF(const std::string filename)
      : energy(), exit_energy(), cdf(), rvalue() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      exit_energy.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss4;
      iss4 << line;
      while (getline(iss4, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      cdf.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss5;
      iss5 << line;
      while (getline(iss5, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      rvalue.push_back(tmp_vec);
    }
    file.close();
  }
  CS_3dENDF(const CS_3dENDF &other)
      : energy(other.energy), exit_energy(other.exit_energy), cdf(other.cdf),
        rvalue(other.rvalue) {}

  CS_3dENDF() : energy(), exit_energy(), cdf(), rvalue() {}

  int find_energy_index(const double e) const {
    int ret = 0;
    while (ret < int(energy.size()) && e > energy[ret]) {
      ret++;
    }
    return ret;
  }
  std::vector<double> sample(const double e, gsl_rng *gen,
                             AtomConsts &coef) const {
    int enindex = find_energy_index(e);
    int densindex = 0;
    double densrng = gsl_rng_uniform(gen);
    if (enindex == 0) {
      while (cdf[enindex][densindex] < densrng) {
        densindex++;
      }
    } else {
      double densrng2 = gsl_rng_uniform(gen);
      double difcheck =
          (e - energy[enindex - 1]) / (energy[enindex] - energy[enindex - 1]);
      if (densrng2 < difcheck) {
        while (cdf[enindex][densindex] < densrng) {
          densindex++;
        }
      } else {
        enindex--;
        while (cdf[enindex][densindex] < densrng) {
          densindex++;
        }
      }
    }
    double outenergycm = exit_energy[enindex][densindex];
    double outrvalue = rvalue[enindex][densindex];
    double ea = e * (coef.ANU / (coef.ANU + coef.iNU));
    double eb = outenergycm * ((coef.BNU + coef.oNU) / (coef.BNU));
    double Ea = ea + coef.Si;
    double Eb = eb + coef.So;
    double X1 = fmin(Ea, 130) * Eb / Ea;
    double X3 = fmin(Ea, 41) * Eb / Ea;
    double aval = (0.04 * X1) + (1.8 * 1e-6 * pow(X1, 3)) +
                  (6.7 * 1e-7 * coef.Mi * coef.Mo * pow(X3, 4));
    double cdfc2 = outrvalue * cosh(aval) - sinh(aval);
    double cdfc1 = 2 * sinh(aval);
    double u2 = gsl_rng_uniform(gen);
    double z = cdfc1 * u2 + cdfc2;
    double z2 = (z + pow(pow(z, 2) - pow(outrvalue, 2) + 1, 1.0 / 2.0)) /
                (outrvalue + 1);
    double outanglecm = log(z2) / aval;
    outanglecm = fmax(outanglecm, -1);
    outanglecm = fmin(outanglecm, 1);
    double outenergylab =
        outenergycm + (e * coef.iNU * coef.oNU / pow(coef.ANU + coef.iNU, 2)) +
        ((2 * pow(coef.iNU * coef.oNU * outenergycm * e, 1.0 / 2.0) *
          outanglecm) /
         (coef.ANU + coef.iNU));
    double outanglelab =
        (pow(outenergycm / outenergylab, 1.0 / 2.0) * outanglecm) +
        (pow(coef.iNU * coef.oNU * e / outenergylab, 1.0 / 2.0) /
         (coef.ANU + coef.iNU));
    if (fabs(outanglelab) > 1) {
      outanglelab = fmin(1, outanglelab);
      outanglelab = fmax(-1, outanglelab);
    }
    std::vector<double> returnarray = {outenergylab, outanglelab};
    return returnarray;
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_energy, cdf, rvalue;
};

#endif
