import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.interpolate

from endf_parserpy import EndfParser


# Setup for proton in proton out angle/energy. Angle-energy data is stored in endf_dict[6][5]['subsection'].
# Particle type can be accessed using endf_dict[6][5]['subsection'][a], where a specifies the output particle
# (e.g. 2 is proton). These can be found in the interpreted ENDF file. Incident energies & their keys can be
# found by appending ['E'], e.g. endf_dict[6][5]['subsection'][2]['E'] for protons. Yields for incident energies
# can be found by appending ['yields']['yi'], e.g. endf_dict[6][5]['subsection'][2]['yields']['yi'] for protons.
# Data for the angle-energy distribution for a given incident energy is obtained by adding ['table'][b], where
# [b] is the associated key for the incident energy, e.g. endf_dict[6][5]['subsection'][2]['table'][b] for protons.
# For each output particle, the angle density is discretized into 19 bins, the specific value of these bins can
# be found in the interpreted ENDF file for any given output particle. To obtain the energy-angle distribution
# for a given angle, we add the [c] option for c in {1,...,19}, e.g. endf_dict[6][5]['subsection'][2]['table'][b][c]
# for a proton. The output energy is also discretized, where for a given incident energy and angle, the discretization
# is in endf_dict[6][5]['subsection'][2]['table'][b][c]['Ep'] and the associated angle-energy density in
# endf_dict[6][5]['subsection'][2]['table'][b]['f']

# The dictionary densities contains key,value pairs of the input energy + a list of lists,
# where each entry of the list is a triple containing the output angle + lists for the
# output energy and normalized densities for those angles.


def DiscreteDensityConstructor(endf_dict, a):
    densities = {}
    incident_energies = endf_dict[6][5]["subsection"][a]["E"].items()
    for key, value in incident_energies:
        angle_keys = endf_dict[6][5]["subsection"][a]["mu"][key].keys()
        total_density = sum(
            [
                sum(endf_dict[6][5]["subsection"][2]["table"][key][i]["f"])
                for i in angle_keys
            ]
        )
        normalised_density = []
        for j in angle_keys:
            unnormalised_density = endf_dict[6][5]["subsection"][a]["table"][key][j][
                "f"
            ]
            output_energy = endf_dict[6][5]["subsection"][a]["table"][key][j]["Ep"]
            output_angle = endf_dict[6][5]["subsection"][a]["mu"][key][j]
            tmp_normalised_density = [p / total_density for p in unnormalised_density]
            counter = 0
            while output_energy[-1] > value:
                del output_energy[-1]
                tmp_normalised_density[-2] += tmp_normalised_density[-1]
                del tmp_normalised_density[-1]
            normalised_density.append(
                [output_angle, output_energy, tmp_normalised_density]
            )
        densities[value] = normalised_density
    return densities


def splineConstructorAngle(k, densities):
    x = np.linspace(-1, 1, 40)
    spline_PDF = {}
    CDF = {}
    for key in densities.keys():
        discrete_pmf = [sum(densities[key][i][2]) for i in range(len(densities[key]))]
        angles = [densities[key][i][0] for i in range(len(densities[key]))]
        spline_PDF[key] = scipy.interpolate.splrep(angles, discrete_pmf, k=k)
        integrated_spline = [
            scipy.interpolate.splint(-1, n, spline_PDF[key]) for n in x
        ]
        integrated_spline = [
            integrated_spline[n] / integrated_spline[-1] for n in range(len(x))
        ]
        for s in range(1, len(integrated_spline)):
            integrated_spline[s] = max(integrated_spline[s], integrated_spline[s - 1])
        CDF[key] = integrated_spline
    return [spline_PDF, CDF]


def splineConstructorEnergy(k, densities):
    spline_PDF = {}
    CDF = {}
    for key in densities.keys():
        spline_PDF_per_angle = {}
        CDF_per_angle = {}
        # densities for output energy with input energy "key" and angle "densities[key][i][0]"
        for i in range(19):
            energies = densities[key][i][1]
            dens = densities[key][i][2]
            angle = densities[key][i][0]
            spline_PDF_per_angle[i] = scipy.interpolate.splrep(energies, dens, k=k)
            x = np.linspace(energies[0], energies[-1], 50)
            integrated_spline = [
                scipy.interpolate.splint(energies[0], n, spline_PDF_per_angle[i])
                for n in x
            ]
            for s in range(1, len(integrated_spline)):
                integrated_spline[s] = max(
                    integrated_spline[s], integrated_spline[s - 1]
                )
            integrated_spline = [
                integrated_spline[n] / integrated_spline[-1] for n in range(len(x))
            ]
            CDF_per_angle[i] = integrated_spline
        spline_PDF[key] = spline_PDF_per_angle
        CDF[key] = CDF_per_angle
    return [spline_PDF, CDF]


# reaction rate cross section for protons
parser = EndfParser()
endf_dict = parser.parsefile("JENDLEPCS.txt")
incident_energy = endf_dict[3][5]["xstable"]["E"]
cross_section = endf_dict[3][5]["xstable"]["xs"]
spline_nonelastic_cs = scipy.interpolate.splrep(incident_energy, cross_section, k=3)


x = np.linspace(incident_energy[0], incident_energy[-1], 100)
# Uncomment for diagnostic plot
# plt.plot(x, [max(z, 0) for z in scipy.interpolate.splev(x, spline_nonelastic_cs)])
# plt.plot(incident_energy, cross_section, "o")
# plt.title("Nonelastic cross section")
# plt.xlabel("Incident Energy")
# plt.ylabel("Cross section")
# plt.savefig("proton_nonelastic_cs.pdf")
# plt.clf()

np.set_printoptions(suppress=True)
f = open("proton_nonelastic_cs.txt", "w")
f.write(" ".join([str(z / 1e6) for z in x]) + "\n")
tmp = " ".join(
    [str(max(z, 0)) for z in scipy.interpolate.splev(x, spline_nonelastic_cs)]
)
f.write(tmp + "\n")
f.close()

endf_dict = parser.parsefile("JENDLEP.txt")
angles = [
    -1,
    -0.9848078,
    -0.9396926,
    -0.8660254,
    -0.7660444,
    -0.6427876,
    -0.5,
    -0.3420201,
    -0.1736482,
    0,
    0.1736482,
    0.3420201,
    0.5,
    0.6427876,
    0.7660444,
    0.8660254,
    0.9396926,
    0.9848078,
    1,
]
densities = DiscreteDensityConstructor(endf_dict, 2)

angle_density_splines = splineConstructorAngle(3, densities)
energy_density_splines = splineConstructorEnergy(3, densities)
energies = [key for key in angle_density_splines[0].keys()]
spline_angle_PDF = angle_density_splines[0]
spline_energy_PDF = energy_density_splines[0]
spline_angle_CDF = angle_density_splines[1]
spline_energy_CDF = energy_density_splines[1]

# Uncomment for diagnostic plot
# x = np.linspace(-1, 1, 40)
# for n in [10, 20, 30]:
#     plt.subplot(1, 2, 1)
#     temp = [sum(densities[energies[n]][i][2]) for i in range(19)]
#     plt.plot(x, scipy.interpolate.splev(x, spline_angle_PDF[energies[n]]))
#     plt.plot(angles, temp, "o")
#     plt.xlabel("Cosine of angle")
#     plt.ylabel("Density")
#     plt.subplot(1, 2, 2)
#     plt.plot(x, spline_angle_CDF[energies[n]])
#     plt.xlabel("Cosine of angle")
#     plt.suptitle("Marginal PDF & CDF of angle at incident energy " + str(int(energies[n] / 1e6)))
#     plt.savefig("proton_marginal_angle_" + str(int(energies[n] / 1e6)) + ".pdf")
#     plt.clf()

# Uncomment for diagnostic plot
# for n in [5, 20, 35]:
#     for i in [6, 12, 18]:
#         temp_energies = densities[energies[n]][i][1]
#         x = np.linspace(temp_energies[0], temp_energies[-1], 50)
#         temp = densities[energies[n]][i][2]
#         angle = densities[energies[n]][i][0]
#         plt.subplot(1, 2, 1)
#         plt.plot(x, scipy.interpolate.splev(x, spline_energy_PDF[energies[n]][i]))
#         plt.plot(temp_energies, temp, "o")
#         plt.xlabel("Exit energy")
#         plt.ylabel("Density")
#         plt.subplot(1, 2, 2)
#         plt.plot(x, spline_energy_CDF[energies[n]][i])
#         plt.suptitle("Conditional energy PDF & CDF at incident energy " + str(int(energies[n] / 1e6)) + " & angle " + str(angle))
#         plt.xlabel("Exit energy")
#         plt.savefig("proton_conditional_energy_" + str(int(energies[n] / 1e6)) + "_" + str(angle) + ".pdf")
#         plt.clf()

f = open("proton_exit_angle_cdf.txt", "w")
tmp = " ".join([str(z / 1e6) for z in spline_angle_CDF.keys()])
f.write(tmp + "\n")
tmp = np.array2string(np.linspace(-1, 1, 40), max_line_width=float("inf"))[1:-1]
tmp = " ".join(tmp.split())
f.write(tmp + "\n")
for key in spline_angle_CDF.keys():
    f.write(" ".join(str(z) for z in spline_angle_CDF[key]) + "\n")
f.close()

f = open("proton_exit_energy_cdf.txt", "w")
tmp = " ".join([str(z / 1e6) for z in spline_energy_CDF.keys()])
f.write(tmp + "\n")
tmp = " ".join([str(z) for z in angles])
f.write(tmp + "\n")
for i in range(len(energies)):
    key = energies[i]
    for j in range(len(angles)):
        key2 = angles[j]
        tmp = [str(z) for z in spline_energy_CDF[key][j]]
        f.write(" ".join(tmp) + "\n")
        temp_energies = densities[key][j][1]
        tmp = " ".join(
            [str(z / 1e6) for z in np.linspace(temp_energies[0], temp_energies[-1], 50)]
        )
        f.write(tmp + "\n")
f.close()
