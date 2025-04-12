import math
import numpy as np
import random
from wright_fisher import wright_fisher
from cross_sections import cross_section_1d
from cross_sections import cross_section_2d
from cross_sections import cross_section_3d


class proton_path:

    def __init__(self, x0, omega0, E0, dt):
        self.dt = dt  # time step
        self.points = [
            [x0, omega0, E0, 0]
        ]  # initialise paths with start position, angle, energy, and deposition
        self.nonelastic_cs = cross_section_1d("../proton_nonelastic_cs.txt")
        self.elastic_cs_oxygen = cross_section_1d("../proton_elastic_oxygen_cs.txt")
        self.elastic_cs_hydrogen = cross_section_1d("../proton_elastic_hydrogen_cs.txt")
        self.nonelastic_angle_cdf = cross_section_2d("../proton_exit_angle_cdf.txt")
        self.nonelastic_energy_cdf = cross_section_3d("../proton_exit_energy_cdf.txt")
        self.elastic_angle_hydrogen_cdf = cross_section_2d(
            "../proton_elastic_angle_hydrogen_cdf.txt"
        )
        self.elastic_angle_oxygen_cdf = cross_section_2d(
            "../proton_elastic_angle_oxygen_cdf.txt"
        )
        self.scatter()

    def rutherford_cross_section(self, e0, lb):
        Z_oxy = 8
        Z_hyd = 1
        A_oxy = 16
        A_hyd = 1
        log_ahbarc = math.log(197.3 / 137) - 13 * math.log(10)
        mpcsq = 938.346  # mass of proton * speed of light squared, MeV
        pv = (2 * mpcsq + e0) * e0 / (mpcsq + e0)
        log_density = 22 * math.log(10) + math.log(3.345)
        # molecules / cm^3
        sig = (
            math.exp(
                2
                * (
                    log_ahbarc
                    + math.log(math.cos(lb / 2))
                    - math.log(pv)
                    - math.log(math.sin(lb / 2))
                )
                + log_density
            )
            * math.pi
        )
        ret = sig * (A_oxy * Z_oxy**2 + 2 * A_hyd * Z_hyd**2) / (A_oxy + 2 * A_hyd)
        return ret

    def rutherford_scatter(self, lb):
        beta = 2 * math.pi * random.random()
        u = random.random()
        alpha = math.acos(
            (math.cos(lb) - u * math.cos(lb / 2) ** 2) / (1 - u * math.cos(lb / 2) ** 2)
        )
        angle = self.points[-1][1]
        if (0 <= beta) and (beta < math.pi / 2):
            angle[0] = angle[0] - math.atan(math.sin(beta) * math.tan(alpha))
            angle[1] = angle[1] - math.atan(math.cos(beta) * math.tan(alpha))
        elif (math.pi / 2 <= beta) and (beta < math.pi):
            angle[0] = angle[0] - math.atan(math.sin(math.pi - beta) * math.tan(alpha))
            angle[1] = angle[1] + math.atan(math.cos(math.pi - beta) * math.tan(alpha))
        elif (math.pi <= beta) and (beta < 3 * math.pi / 2):
            angle[0] = angle[0] + math.atan(
                math.sin(3 * math.pi / 2 - beta) * math.tan(alpha)
            )
            angle[1] = angle[1] + math.atan(
                math.cos(3 * math.pi / 2 - beta) * math.tan(alpha)
            )
        else:
            angle[0] = angle[0] + math.atan(
                math.sin(2 * math.pi - beta) * math.tan(alpha)
            )
            angle[1] = angle[1] - math.atan(
                math.cos(2 * math.pi - beta) * math.tan(alpha)
            )
        self.points[-1][1] = angle
        return

    def nonelastic_cross_section(self, e0):
        ret = self.nonelastic_cs.evaluate(e0)
        log_water_density = 22 * math.log(10) + math.log(3.345)  # molecules / cm^3
        log_barns_to_cmsq = -24 * math.log(10)
        A_oxy = 16
        A_hyd = 1
        c = (
            A_oxy
            * math.exp(log_barns_to_cmsq + log_water_density)
            / (A_oxy + 2 * A_hyd)
        )
        return c * ret  # rate per cm

    def elastic_cross_section(self, e0):
        elastic_oxy = self.elastic_cs_oxygen.evaluate(e0)
        elastic_hyd = self.elastic_cs_hydrogen.evaluate(e0)
        log_water_density = 22 * math.log(10) + math.log(3.345)  # molecules / cm^3
        log_barns_to_cmsq = -24 * math.log(10)
        A_oxy = 16
        A_hyd = 1
        c = math.exp(log_barns_to_cmsq + log_water_density)
        ret = c * (A_oxy * elastic_oxy + 2 * A_hyd * elastic_hyd) / (A_oxy + 2 * A_hyd)
        return c, elastic_oxy, elastic_hyd

    def nonelastic_scatter(self, omega, e0):
        ang = omega
        beta = 2 * math.pi * random.random()
        alpha = self.nonelastic_angle_cdf.sample(e0)
        if (0 <= beta) and (beta < math.pi / 2):
            ang[0] = ang[0] - math.atan(math.sin(beta) * math.tan(alpha))
            ang[1] = ang[1] - math.atan(math.cos(beta) * math.tan(alpha))
        elif (math.pi / 2 <= beta) and (beta < math.pi):
            ang[0] = ang[0] - math.atan(math.sin(math.pi - beta) * math.tan(alpha))
            ang[1] = ang[1] + math.atan(math.cos(math.pi - beta) * math.tan(alpha))
        elif (math.pi <= beta) and (beta < 3 * math.pi / 2):
            ang[0] = ang[0] + math.atan(
                math.sin(3 * math.pi / 2 - beta) * math.tan(alpha)
            )
            ang[1] = ang[1] + math.atan(
                math.cos(3 * math.pi / 2 - beta) * math.tan(alpha)
            )
        else:
            ang[0] = ang[0] + math.atan(math.sin(2 * math.pi - beta) * math.tan(alpha))
            ang[1] = ang[1] - math.atan(math.cos(2 * math.pi - beta) * math.tan(alpha))
        return alpha, ang

    def elastic_scatter(self, omega, e0, elastic_oxy, elastic_hyd):
        ang = omega
        beta = 2 * math.pi * random.random()
        alpha = 0
        A_oxy = 16
        A_hyd = 1
        if random.random() < A_oxy * elastic_oxy / (
            A_oxy * elastic_oxy + 2 * A_hyd * elastic_hyd
        ):
            alpha = self.elastic_angle_oxygen_cdf.sample(e0)
        else:
            alpha = self.elastic_angle_hydrogen_cdf.sample(e0)
        if (0 <= beta) and (beta < math.pi / 2):
            ang[0] = ang[0] - math.atan(math.sin(beta) * math.tan(alpha))
            ang[1] = ang[1] - math.atan(math.cos(beta) * math.tan(alpha))
        elif (math.pi / 2 <= beta) and (beta < math.pi):
            ang[0] = ang[0] - math.atan(math.sin(math.pi - beta) * math.tan(alpha))
            ang[1] = ang[1] + math.atan(math.cos(math.pi - beta) * math.tan(alpha))
        elif (math.pi <= beta) and (beta < 3 * math.pi / 2):
            ang[0] = ang[0] + math.atan(
                math.sin(3 * math.pi / 2 - beta) * math.tan(alpha)
            )
            ang[1] = ang[1] + math.atan(
                math.cos(3 * math.pi / 2 - beta) * math.tan(alpha)
            )
        else:
            ang[0] = ang[0] + math.atan(math.sin(2 * math.pi - beta) * math.tan(alpha))
            ang[1] = ang[1] - math.atan(math.cos(2 * math.pi - beta) * math.tan(alpha))
        return alpha, ang

    def multiple_scattering_sd(self, e0):
        Z_oxy = 8
        Z_hyd = 1
        A_oxy = 16
        A_hyd = 1
        mpcsq = 938.346  # mass of proton * speed of light squared, MeV
        # radiation length of oxygen, g / cm^2
        x_oxy = A_oxy * 716.4 / (Z_oxy * (Z_oxy + 1) * math.log(287 / math.sqrt(Z_oxy)))
        # radiation length of hydrogen, g / cm^2
        x_hyd = A_hyd * 716.4 / (Z_hyd * (Z_hyd + 1) * math.log(287 / math.sqrt(Z_hyd)))
        # radiation length of water via Bragg additivity rule
        x0 = (2 * A_hyd + A_oxy) * x_oxy * x_hyd / (A_oxy * x_hyd + 2 * A_hyd * x_oxy)
        pv = (2 * mpcsq + e0) * e0 / (mpcsq + e0)
        ret = 14.1 * math.sqrt(self.dt / x0) * (1 + math.log(self.dt / x0, 10) / 9) / pv
        return ret

    def bethe_bloch(self, e0):
        rho = 1  # density of medium, g / cm^3
        Z_oxy = 8  # atomic number
        Z_hyd = 1  # atomic number
        A_oxy = 16  # atomic mass
        A_hyd = 1  # atomic mass
        mecsq = 0.511  # mass of electron * speed of light squared, MeV
        mpcsq = 938.346  # mass of proton * speed of light squared, MeV
        # mean excitation energy of medium, MeV, fitted by eye
        I = 60 * 1e-6
        betasq = (2 * mpcsq + e0) * e0 / (mpcsq + e0) ** 2
        b_oxy = (
            0.3072
            * Z_oxy
            * rho
            * (math.log(2 * mecsq * betasq / I * (1 - betasq)) - betasq)
            / (A_oxy * betasq)
        )  # MeV / cm
        b_hyd = (
            0.3072
            * Z_hyd
            * rho
            * (math.log(2 * mecsq * betasq / I * (1 - betasq)) - betasq)
            / (A_hyd * betasq)
        )  # MeV / cm
        ret = (A_oxy * b_oxy + 2 * A_hyd * b_hyd) / (A_oxy + 2 * A_hyd)
        return ret

    def energy_straggling_sd(self):
        alpha = 1 / 137  # fine structure constant
        log_hbar = (
            -21 * math.log(10) + math.log(4.136) - math.log(2 * math.pi)
        )  # MeV * s
        log_c = math.log(29979245800)  # cm / s
        log_water_density = 22 * math.log(10) + math.log(3.345)  # molecules / cm^3
        z = 10  # electrons per water molecule
        ret = (
            4
            * math.pi
            * z
            * math.exp(2 * (math.log(alpha) + log_hbar + log_c) + log_water_density)
        )
        return math.sqrt(ret)

    def scatter(self):
        absorption_e = 0.037
        while self.points[-1][2] > absorption_e:
            self.spherical_BM(
                self.points[-1][0], self.points[-1][1], self.points[-1][2], self.dt
            )
            if self.points[-1][2] > absorption_e:
                non_e_jump = self.nonelastic_cross_section(self.points[-1][2])
                e_jump, e_oxy, e_hyd = self.elastic_cross_section(self.points[-1][2])
                rutherford_lb = 2.5 * self.multiple_scattering_sd(self.points[-1][2])
                rutherford_jump = 0
                if rutherford_lb < math.pi:
                    rutherford_jump = self.rutherford_cross_section(
                        self.points[-1][2], rutherford_lb
                    )
                alpha = non_e_jump + e_jump + rutherford_jump
                if random.random() < 1 - math.exp(-self.dt * alpha):
                    u = random.random()
                    if u < rutherford_jump / alpha:
                        self.rutherford_scatter(rutherford_lb)
                    elif u < (rutherford_jump + e_jump) / alpha:
                        exit_cos, self.points[-1][1] = self.elastic_scatter(
                            self.points[-1][1], self.points[-1][2], e_oxy, e_hyd
                        )
                    else:
                        exit_cos, self.points[-1][1] = self.nonelastic_scatter(
                            self.points[-1][1], self.points[-1][2]
                        )
                        self.points[-1][2] = self.nonelastic_energy_cdf.sample(
                            self.points[-1][2], exit_cos
                        )
                        self.points[-1][3] = self.points[-2][2] - self.points[-1][2]
        return self.points

    def spherical_BM(self, x, omega, e0, dt):
        wf = wright_fisher(self.multiple_scattering_sd(e0) ** 2)
        theta = 2 * math.pi * random.random()
        Y = np.array((np.cos(theta), np.sin(theta)))

        # Set u = (e_d - z) / |e_d - z|
        z = [
            [
                np.sin(omega[0]) * np.cos(omega[1]),
                np.sin(omega[0]) * np.sin(omega[1]),
                np.cos(omega[0]),
            ]
        ]
        e_d = np.zeros(3)
        e_d[2] = 1.0
        dist = np.linalg.norm(e_d - z)
        if dist > 1e-9:
            u = (e_d - z) / dist
        else:
            u = np.ones(3)
            if z[0] > 0:
                u[0] = -1
            if z[1] > 0:
                u[1] = -1

        # Compute O(z) = I - 2uu_T
        O_z = np.identity(3) - 2 * np.outer(u, u)
        c = np.concatenate(
            (np.array([2 * np.sqrt(wf * (1 - wf)) * Y]), np.array([1 - 2 * wf])),
            axis=None,
        )
        w = np.dot(O_z, c)

        angle = [np.arccos(w[2]), np.arctan2(w[1], w[0])]
        pos = [
            x[0] + np.sin(angle[0]) * np.cos(angle[1]) * dt,
            x[1] + np.sin(angle[0]) * np.sin(angle[1]) * dt,
            x[2] + np.cos(angle[0]) * dt,
        ]
        delta = np.linalg.norm(np.array(self.points[-1][0]) - np.array(pos))
        energy_mu = self.bethe_bloch(e0) * delta
        energy_sd = math.sqrt(delta) * self.energy_straggling_sd()
        energy = e0 - random.gauss(energy_mu, energy_sd)

        self.points.append([pos, angle, energy, e0 - energy])
