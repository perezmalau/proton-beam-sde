# Geant4 Energy Deposition in a Cubic Phantom

This repository contains a **Geant4 simulation** used to obtain 3D energy deposition distributions from **monoenergetic proton beams** in a water phantom.
The simulations serve as the **benchmark reference for the SDE model** used in the accompanying study.

This simulation is part of the MaThRad collaboration work on stochastic differential equation (SDE) models for proton therapy radiation transport.

---

## Overview

- The setup consists of a **20 × 20 × 20 cm³ water phantom** irradiated by a **monoenergetic proton beam**.
- **Energy deposition** is scored using a **custom sensitive detector** that records hit positions and deposited energy.
- Positions are converted to voxel indices to generate a 3D dose map.
- The geometry can be modified to include **slabs of air** and/or **bone** of specific thicknesses using UI commands in the macro file.
- Optional tracking of **secondary particles** and **hadElastic processes** is implemented but **disabled by default**.

---

## Optional Features

Additional information about **secondary particles** and **hadronic elastic interactions** can be stored in a ROOT file.  
These features are **commented out by default** in:

- `TrackingAction.cc`
- `RunAction.cc`
- `SteppingAction.cc`

Uncomment the relevant sections to enable these outputs. To change the root file name, add `/output/rootFileName test.root` in your macro file.

---

## Physics Lists

To change the physics list, edit `ProtonTherapy.cc` by uncommenting the desired list. The default is QGSP_BIC_EMZ:

```cpp
// Physics lists — uncomment the one you need
auto* physicsFactory = new G4PhysListFactory();

// Reference
auto physicsList = physicsFactory->GetReferencePhysList("QGSP_BIC_EMZ");

// Other physics lists (comment line above and uncomment one of these)
// auto physicsList = physicsFactory->GetReferencePhysList("QGSP_BIC_EMY");
// auto physicsList = physicsFactory->GetReferencePhysList("QGSP_BERT");
```
## Simulation Macros

Five main macros reproduce the configurations used in the study, plus a test run:

| Macro file | Description |
|-------------|--------------|
| `run0.mac` | Test run |
| `run_air.mac` | 5 cm air gap + 2 cm bone + water, 100 MeV protons |
| `run_bone_100MeV.mac` | Bone slab only, 100 MeV protons |
| `run_bone_150MeV.mac` | Bone slab only, 150 MeV protons |
| `run_water_100MeV.mac` | Pure water phantom, 100 MeV protons |
| `run_water_150MeV.mac` | Pure water phantom, 150 MeV protons |

---

## Compilation

Ensure your Geant4 environment is sourced correctly before compiling:

```bash
source /path/to/geant4-install/bin/geant4.sh
```

Then create a `build/` directory and compile from that folder:
```bash
mkdir build
cd build
cmake ../
make -j$(nproc)
```
## Running the Simulation

After compilation, run the executable by specifying the macro file and number of threads:

```bash
./ProtonTherapy -m <macrofile.mac> -t <nthreads>
```

Each run will generate output `.csv` files containing the **3D energy deposition data** in the phantom for post-processing using Python scripts. 
The output files are written to the `Output/` directory. To change it, modify the path of the output file names in the macro files.

## Notes

- This code was developed and tested using **Geant4 version 11.3.0**. Detailed instructions on how to install Geant4 and its prerequisites can be found in https://geant4.web.cern.ch/download/11.3.0.html.
- ROOT output may grow large, disable secondary tracking if not required.
- To record simulation time, run by adding `time` right before the main running command. Example: `time ./ProtonTherapy -m run_bone_150MeV.mac -t 8`.