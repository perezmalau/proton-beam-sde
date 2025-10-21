from SDE_vs_G4 import (retrieve_sde_output, retrieve_g4_output, plot_lateral_profiles,
                       plot_slice, pymedphys_gamma, compare_bragg_peaks, include_new_material,
                       range_comparison, plot_bragg_peak)
import numpy as np

R = 0.04
BS = 0.04
N_primaries = "1E6"
# suffix = "_bone+h2o" # Heterogeneous, bone and water
suffix = "" # Homogeneous
# suffix = "hetero" # Air, bone and water --this doesn't have a 150 MeV version!!
energy1 = 100 # MeV
energy2 = 150 # MeV

voxel_volume = 0.001 # cm3
density_matrix = np.ones((200, 200, 200)) # density of water 1 g/cm3
cuts_100 = [30, 50, 75] # These are the depths at which we plot lateral profiles
cuts_150 = [50, 100, 150]
# Adjust parameters according to case:
if suffix == "_bone+h2o":
    density_matrix = include_new_material(density_matrix, tmin=0, tmax=20, rho=1.45) # add bone
    cuts_100 = [10, 40, 65]
    cuts_150 = [10, 50, 145]
elif suffix == "hetero":
    density_matrix = include_new_material(density_matrix, tmin=0, tmax=50, rho=0.001225) # add air
    density_matrix = include_new_material(density_matrix, tmin=50, tmax=60, rho=1.45)  # add bone
    cuts_100 = [25, 55, 120]

mass_matrix = density_matrix * voxel_volume

g4_dose1 = (retrieve_g4_output(f"Geant4_data/G4_TotalEneDepMap_{N_primaries}_{energy1}MeV{suffix}.csv")
            / (float(N_primaries)*mass_matrix))
g4_dose2 = (retrieve_g4_output(f"Geant4_data/G4_TotalEneDepMap_{N_primaries}_{energy2}MeV{suffix}.csv")
            / (float(N_primaries)*mass_matrix))

sde_dose1 = (retrieve_sde_output(f"SDE_B={BS}_R={R}_{N_primaries}_{energy1}MeV{suffix}.txt")
             / (float(N_primaries)*mass_matrix))
sde_dose2 = (retrieve_sde_output(f"SDE_B={BS}_R={R}_{N_primaries}_{energy2}MeV{suffix}.txt")
             / (float(N_primaries)*mass_matrix))

# Plots for Energy1:
f1 = plot_slice(g4_dose1, sde_dose1, 'z', minval=-5, xmin=0, xmax=9, ymin=-5, ymax=5, zmin=-5, zmax=5)
f2 = plot_lateral_profiles(g4_dose1, sde_dose1, xmin=0, xmax=10, ymin=-5, ymax=5, zmin=-5, zmax=5, lowlim=-1,
                               uplim=1, cuts=cuts_100)
print("Gamma values for 100 MeV: ")
_ = pymedphys_gamma(g4_dose1, sde_dose1, dta=2, dd=3, th_percent=10, is_local=False, xmin=0, xmax=9,
                       ymin=-2.5, ymax=2.5, zmin=-5, zmax=5)
_ = pymedphys_gamma(g4_dose1, sde_dose1, dta=2, dd=3, th_percent=10, is_local=True, xmin=0, xmax=9,
                       ymin=-2.5, ymax=2.5, zmin=-5, zmax=5)
f3 = pymedphys_gamma(g4_dose1, sde_dose1, dta=2, dd=3, th_percent=1, is_local=True, xmin=0, xmax=9,
                        ymin=-2.5, ymax=2.5, zmin=-5, zmax=5)
_ = pymedphys_gamma(g4_dose1, sde_dose1, dta=1, dd=2, th_percent=1, is_local=True, xmin=0, xmax=9,
                       ymin=-2.5, ymax=2.5, zmin=-5, zmax=5)

f1.savefig(f"results/2DSlice_z_{energy1}MeV{suffix}.png")
f2.savefig(f"results/lateral_profiles_{energy1}MeV{suffix}.png")
f3.savefig(f"results/gammaTest_{energy1}MeV{suffix}.png")

# Plots for Energy2:
f4 = plot_slice(g4_dose2, sde_dose2, 'z', minval=-5, xmin=0, xmax=18, ymin=-10, ymax=10, zmin=-10, zmax=10)
f5 = plot_lateral_profiles(g4_dose2, sde_dose2, xmin=0, xmax=18, ymin=-5, ymax=5, zmin=-5, zmax=5, lowlim=-1, uplim=1,
                        cuts=cuts_150)
print("Gamma values for 150 MeV: ")
_ = pymedphys_gamma(g4_dose2, sde_dose2, dta=2, dd=3, th_percent=10, is_local=False, xmin=0, xmax=18,
ymin=-5, ymax=5, zmin=-5, zmax=5)
_ = pymedphys_gamma(g4_dose2, sde_dose2, dta=2, dd=3, th_percent=10, is_local=True, xmin=0, xmax=18,
ymin=-5, ymax=5, zmin=-5, zmax=5)
f6 = pymedphys_gamma(g4_dose2, sde_dose2, dta=2, dd=3, th_percent=1, is_local=True, xmin=0, xmax=18,
                       ymin=-5, ymax=5, zmin=-5, zmax=5)
_ = pymedphys_gamma(g4_dose2, sde_dose2, dta=1, dd=2, th_percent=1, is_local=True, xmin=0, xmax=18,
ymin=-5, ymax=5, zmin=-5, zmax=5)

f4.savefig(f"results/2DSlice_z_{energy2}MeV{suffix}.png")
f5.savefig(f"results/lateral_profiles_{energy2}MeV{suffix}.png")
f6.savefig(f"results/gammaTest_{energy2}MeV{suffix}.png")

# # Comparative plots in 1D:
f7, f8 = compare_bragg_peaks(sde_dose1, g4_dose1, sde_dose2, g4_dose2, energy1, energy2)
f7.savefig(f"results/1DProj_{energy1}_vs_{energy2}MeV{suffix}.png")
f8.savefig(f"results/1DSlice_{energy1}_vs_{energy2}MeV{suffix}.png")

# Proton range numbers:
r90_1 = range_comparison(g4_dose1, sde_dose1, 0.9)
r50_1 = range_comparison(g4_dose1, sde_dose1, 0.5)
print(f"For {energy1} MeV, proton range agreement is within {r90_1*10:.2f} mm (R90) and {r50_1*10:.2f} mm (R50)")
r90_2 = range_comparison(g4_dose2, sde_dose2, 0.9)
r50_2 = range_comparison(g4_dose2, sde_dose2, 0.5)
print(f"For {energy2} MeV, proton range agreement is within {r90_2*10:.2f} mm (R90) and {r50_2*10:.2f} mm (R50)")

# Comparison with other physics lists, only if phantom is homogeneous (i.e. file suffix is ""):
if suffix == "":
    g4_dose_emy1 = (retrieve_g4_output(f"Geant4_data/G4_TotalEneDepMap_{N_primaries}_{energy1}MeV_emy.csv")
                / (float(N_primaries)*mass_matrix))
    g4_dose_bert1 = (retrieve_g4_output(f"Geant4_data/G4_TotalEneDepMap_{N_primaries}_{energy1}MeV_bert.csv")
                / (float(N_primaries)*mass_matrix))
    g4_dose_emy2 = (retrieve_g4_output(f"Geant4_data/G4_TotalEneDepMap_{N_primaries}_{energy2}MeV_emy.csv")
                / (float(N_primaries)*mass_matrix))
    g4_dose_bert2 = (retrieve_g4_output(f"Geant4_data/G4_TotalEneDepMap_{N_primaries}_{energy2}MeV_bert.csv")
                / (float(N_primaries)*mass_matrix))
    f9 = plot_bragg_peak(g4_dose1, g4_dose_emy1, g4_dose_bert1, sde_dose1, how='proj',
                         names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
                        xmax=8, ref_name='QGSP_BIC_EMZ', maxdif=10)
    f10 = plot_bragg_peak(g4_dose1, g4_dose_emy1, g4_dose_bert1, sde_dose1, how='slice',
                          names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
                         xmax=8, ref_name='QGSP_BIC_EMZ')

    f11 = plot_bragg_peak(g4_dose2, g4_dose_emy2, g4_dose_bert2, sde_dose2, how='proj',
                          names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
                        xmax=17, ref_name='QGSP_BIC_EMZ', maxdif=25)
    f12 = plot_bragg_peak(g4_dose2, g4_dose_emy2, g4_dose_bert2, sde_dose2, how='slice',
                          names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
                         xmax=17, ref_name='QGSP_BIC_EMZ', maxdif=30)

    f9.savefig(f"results/1DProj_multiComparison_{energy1}MeV.png")
    f10.savefig(f"results/1DSlice_multiComparison_{energy1}MeV.png")
    f11.savefig(f"results/1DProj_multiComparison_{energy2}MeV.png")
    f12.savefig(f"results/1DSlice_multiComparison_{energy2}MeV.png")