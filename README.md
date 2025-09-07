# Proton-beam-SDE
This repository houses a simulator for an SDE describing a the path of a proton in proton beam therapy.
It is under active development and **not fit for standalone use**.

# Installation

In Unix-like environments, simply call `make` in the project root. You will need:

- the g++ compiler (tested on version 13.2.0),
- the [GNU Scientific Library](https://www.gnu.org/software/gsl/) (tested on version 2.7.1),
- the [libconfig](https://hyperrealm.github.io/libconfig/) package (tested on version 1.5).

# Usage

Simulation and output parameters are specified in the config file `sim.cfg`. In particular, you can
modify:
- the number of protons to simulate,
- the radius of the nozzle from which protons are emitted,
- the Gaussian distribution of the initial position of the protons, conditioned to lie in the nozzle,
- the Gaussian distribution of the initial energy of the protons,
- the change points between materials,
- the index of the material in each gap between change points, based on the ordering in Materials/materials.txt,
- the simulator time-step,
- the energy at which a proton is absorbed and its path ends,
- the side length of output voxelation,
- the path to which a sparse 3d grid of dose deposition is stored.

After compilation, run the simulation by calling `./simulate sim.cfg` at the project root.
You can also replace the `sim.cfg` argument with the path to any other appropriate config file.
**Warning**: The default number of protons in `sim.cfg` is set to one million, and voxel sizes and time
steps are also set to a moderate resolution. As a result, the simulation can take 1-2 minutes.

# Post-processing

Default output paths are specified into the `Output` folder, which contains the `plot.py` script for output
visualisation. To use it you will need:

- [Python 3](https://www.python.org/) (tested on version 3.12.3), as well as the following Python modules:
- numpy (tested on version 1.26.4),
- matplotlib (tested on version 3.10.0),
- itk (tested on version 5.4.4).

If you have used the default number of protons, as well as output modes and paths in `sim.cfg`, you can simply
run `plot.py` to produce heatmaps of the integrated 2d dose, 2d slices through the beam in all three coordinates,
the integrated 1d dose, as well as the 1d dose through the central voxel of the beam.

# Python

The `Python` folder contains an equivalent Python implementation, which still requires precomputed cross
sections stored in the proect root. The Python implementation is much slower (100k replicates in about an hour
as opposed to a million in 2 minutes), but useful for prototyping. See `Python/plots_3d.py` for example usage.

# Nuclear cross sections

The `Splines` folder contains splines fitted to nuclear cross section data for hydroge, oxygen, and nitrogen.
Eventually these will get replaced with a Python script to generate the same cross sections.

# Materials

The interface for specifying materials is via the text files in the `Materials` folder.

- The `atoms.txt` file lists all the atoms needed for all desired materials. Each row specifies an atom name as
well as its atomic number and mass.
- The `materials.txt` lists names of all desired materials. Each material name must be accompanied by a corresponding
`name.txt` file.
- Each `name.txt` file lists three pieces of information. The first two rows contain the material density in `g / cm^3`
and mean excitation energy in `eV`. Mean excitation energies is not well-determined for all material, and here they
also subsume some correction terms to the Bethe-Bloch formula. Hence, they should be regarded as a tuning parameter
rather than a physical constant. Rows beyond the third specify the chemical composition of the material. Each of these
rows consists of an atom ID and a corresponding stoichiometric amount. The numerical IDs refer to the corresponding rows
of `atoms.txt`, counting from zero. The amounts provide the mean amount of that atom per unit volume.

For example, if the rows of `atoms.txt` are:

hydrogen 1 1\
nitrogen 7 14\
oxygen 8 16

then the stoichiometric rows of water are specified as

0 2.0\
2 1.0

while a good approximation for air is

1 1.57\
2 0.43

because air is 78% N<sub>2</sub> and 21% O<sub>2</sub>, and $0.78 \times 2 \approx 1.57$, allowing for slight rounding to
account for the missing 1% of matter.
