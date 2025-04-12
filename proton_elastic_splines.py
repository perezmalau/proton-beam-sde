import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.interpolate

from endf_parserpy import EndfParser


def DiscreteDensityConstructor(endf_dict):
    densities = {}
    incident_energies = endf_dict[6][2]["subsection"][1]["E"].items()
    for key, value in incident_energies:
        AngleProb = endf_dict[6][2]["subsection"][1]["A"][key].items()
        temp1 = []  # store angles
        temp2 = []  # store P
        for key2, value2 in AngleProb:
            if key2 % 2 == 1:
                temp1.append(value2)
            else:
                temp2.append(max(value2, 0))
        temp1[-1] = 1  # final angle is slightly away from one, converting to 1 for ease
        densities[value] = [temp1, temp2]
    return densities


def splineConstructorAngle(k, densities):
    x = np.linspace(-1, 1, 40)
    spline_PDF = {}
    CDF = {}
    for key in densities.keys():
        splt = scipy.interpolate.splrep(densities[key][0], densities[key][1], k=3)
        spline_PDF[key] = splt
        integrated_spline = []
        for n in x:  # interpolates integral of spline to obtain discrete CDF
            integrated_spline.append(scipy.interpolate.splint(-1, n, splt))
        for s in range(1, len(integrated_spline)):
            integrated_spline[s] = max(integrated_spline[s], integrated_spline[s - 1])
        integrated_spline = [
            integrated_spline[n] / integrated_spline[-1] for n in range(len(x))
        ]
        CDF[key] = integrated_spline
    return [spline_PDF, CDF]


parser = EndfParser()

endf_dictCS = parser.parsefile("JENDLEELCSO16.txt")
incident_energy = endf_dictCS[3][2]["xstable"]["E"]
cross_section = endf_dictCS[3][2]["xstable"]["xs"]
cross_section = [max(x, 0) for x in cross_section]
spline_elastic_cs = scipy.interpolate.splrep(incident_energy, cross_section, k=3)

x = np.linspace(incident_energy[0], incident_energy[-1], 100)
# Uncomment for diagnostic plot
# plt.xlim(0, 3e8)
# plt.plot(x, [max(z, 0) for z in scipy.interpolate.splev(x, spline_elastic_cs)])
# plt.plot(incident_energy, cross_section, "o")
# plt.title("Nuclear elastic cross section")
# plt.xlabel("Incident Energy")
# plt.ylabel("Cross section")
# plt.savefig("proton_elastic_cs.pdf")
# plt.clf()

np.set_printoptions(suppress=True)
f = open("proton_elastic_oxygen_cs.txt", "w")
f.write(" ".join([str(z / 1e6) for z in x]) + "\n")
tmp = " ".join([str(max(z, 0)) for z in scipy.interpolate.splev(x, spline_elastic_cs)])
f.write(tmp + "\n")
f.close()

endf_dict = parser.parsefile("JENDLEELO16.txt")
densities = DiscreteDensityConstructor(endf_dict)
angle_density_splines = splineConstructorAngle(3, densities)
energies = [key for key in angle_density_splines[0].keys()]
spline_angle_PDF = angle_density_splines[0]
spline_angle_CDF = angle_density_splines[1]

x = np.linspace(-1, 1, 40)
# Uncomment for diagnostic plot
# for n in range(len(energies)):
#    plt.subplot(1, 2, 1)
#    plt.plot(x, scipy.interpolate.splev(x, spline_angle_PDF[energies[n]]))
#    plt.plot(densities[energies[n]][0], densities[energies[n]][1], "o")
#    plt.xlabel("Cosine of angle")
#    plt.ylabel("Density")
#    plt.subplot(1, 2, 2)
#    plt.plot(x, spline_angle_CDF[energies[n]])
#    plt.xlabel("Cosine of angle")
#    plt.suptitle("PDF & CDF of angle at incident energy " + str(int(energies[n] / 1e6)))
#    plt.savefig("proton_angle_" + str(int(energies[n] / 1e6)) + ".pdf")
#    plt.clf()

f = open("proton_elastic_angle_oxygen_cdf.txt", "w")
tmp = " ".join([str(z / 1e6) for z in spline_angle_CDF.keys()])
f.write(tmp + "\n")
tmp = np.array2string(np.linspace(-1, 1, 40), max_line_width=float("inf"))[1:-1]
tmp = " ".join(tmp.split())
f.write(tmp + "\n")
for key in spline_angle_CDF.keys():
    f.write(" ".join(str(z) for z in spline_angle_CDF[key]) + "\n")
f.close()

endf_dict = parser.parsefile("JENDLEELH1.txt")
total_cs = []
raw_CS = {}
PDF = {}
CDF = {}
incident_energies = endf_dict[6][2]["subsection"][1]["E"].items()
x = np.linspace(-1, 1, 40)
energies = []
for j, en in incident_energies:
    coeffs = endf_dict[6][2]["subsection"][1]["A"][j].items()
    cross_section = []
    for k in range(len(x)):
        tmp = 0
        for l, v in coeffs:
            tmp = (
                tmp
                + (4 * (l - 1) + 1)
                * v
                * scipy.special.eval_legendre(2 * (l - 1), x[k])
                / 2
            )
        cross_section.append(tmp)
    raw_CS[en] = cross_section
    spline_elastic_cs = scipy.interpolate.splrep(x, cross_section, k=3)
    PDF[en] = spline_elastic_cs
    integrated_spline = []
    for n in x:
        integrated_spline.append(scipy.interpolate.splint(-1, n, spline_elastic_cs))
    total_cs.append(integrated_spline[-1])
    energies.append(en)
    for s in range(1, len(x)):
        integrated_spline[s] = max(integrated_spline[s], integrated_spline[s - 1])
    integrated_spline = [
        integrated_spline[n] / integrated_spline[-1] for n in range(len(x))
    ]
    CDF[en] = integrated_spline
spline_elastic_cs = scipy.interpolate.splrep(energies, total_cs, k=3)

y = np.linspace(energies[0], energies[-1], 100)
# Uncomment for diagnostic plot
# plt.xlim(0, 3e8)
# plt.plot(y, [max(z, 0) for z in scipy.interpolate.splev(y, spline_elastic_cs)])
# plt.plot(energies, total_cs, "o")
# plt.title("Nuclear elastic cross section (H)")
# plt.xlabel("Incident Energy")
# plt.ylabel("Cross section")
# plt.savefig("proton_elastic_hydrogen_cs.pdf")
# plt.clf()

np.set_printoptions(suppress=True)
f = open("proton_elastic_hydrogen_cs.txt", "w")
f.write(" ".join([str(z / 1e6) for z in y]) + "\n")
tmp = " ".join([str(max(z, 0)) for z in scipy.interpolate.splev(y, spline_elastic_cs)])
f.write(tmp + "\n")
f.close()

x = np.linspace(-1, 1, 40)
# Uncomment for diagnostic plot
# for n in range(len(energies)):
#    plt.subplot(1, 2, 1)
#    plt.plot(x, scipy.interpolate.splev(x, PDF[energies[n]]))
#    plt.plot(x, raw_CS[energies[n]], "o")
#    plt.xlabel("Cosine of angle")
#    plt.ylabel("Density")
#    plt.subplot(1, 2, 2)
#    plt.plot(x, CDF[energies[n]])
#    plt.xlabel("Cosine of angle")
#    plt.suptitle("PDF & CDF of angle at incident energy " + str(int(energies[n] / 1e6)))
#    plt.savefig("proton_angle_hydrogen_" + str(int(energies[n] / 1e6)) + ".pdf")
#    plt.clf()

f = open("proton_elastic_angle_hydrogen_cdf.txt", "w")
tmp = " ".join([str(z / 1e6) for z in CDF.keys()])
f.write(tmp + "\n")
tmp = np.array2string(np.linspace(-1, 1, 40), max_line_width=float("inf"))[1:-1]
tmp = " ".join(tmp.split())
f.write(tmp + "\n")
for key in CDF.keys():
    f.write(" ".join(str(z) for z in CDF[key]) + "\n")
f.close()
