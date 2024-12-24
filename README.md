# Proton-beam-SDE
This repository houses a simulator for an SDE describing a the path of a proton in proton beam therapy.

# Installation

In Unix-like environments, simply call `make` in the project root. You will need:

- the g++ compiler (tested on version 13.2.0),
- the [GNU Scientific Library](https://www.gnu.org/software/gsl/) (tested on version 2.7.1),
- the [libconfig](https://hyperrealm.github.io/libconfig/) package (tested on version 1.5).

# Precomputing cross sections

Run `proton_nonelastic_splines.py` in the project root. The root folder must contain nonelastic cross
section data for collision of protons in water in the ENDF-6 format, available from [ENDF](https://www.nndc.bnl.gov/endf/).
Once the resulting cross section splines have been computed, this step need not be executed again.

# Usage

Simulation and output parameters are specified in the config file `sim.cfg`. In particular, you can
modify:
- the number of protons to simulate,
- the radius of the nozzle from which protons are emitted,
- the distribution of the initial energy of the protons,
- the simulator time-step,
- the energy at which a proton is absorbed and its path ends,
- the side length of output voxelation (**Warning**: very small voxels result in a large memory cost,
  especially if 3d output is desired),
- desired output modes:
  * `output_3d = 1` outputs the full 3d voxelation of dose deposition,
  * `output_2d = 1` and `output_1d = 1` output the dose with 1 or 2 dimensions integrated out,
  * `output_2d_slice = 1` returns the 2d dose deposition along a one-voxel slice at `z = 0`.
  * Setting any output mode to 0 means that output is not produced. A 3d voxel grid is only created if
    `output_3d = 1`, which slows down the simulation and increases memory cost significantly.
  * The `out_path` variables specify files in which the outputs should be stored.

After compilation, run the simulation by calling `./simulate sim.cfg` at the project root.
You can also replace the `sim.cfg` argument with the path to any other appropriate config file.
**Warning**: The default number of protons in `sim.cfg` is set to one million, and voxel sizes and time
steps are also set to a moderate resolution. As a result, the simulation can take around two minutes.

# Post-processing

Visualisation tools for `output_2d` `output_1d`, and `output_2d_slice` are provided in `plots.R`. To use
them you will need:

- [R](https://www.r-project.org/) (tested on version 4.3.3),
- the [viridis](https://cran.r-project.org/web/packages/viridis/index.html) package (tested on version 0.6.5),
- while not strictly necessary, [RStudio](https://posit.co/) or a similar IDE will also make life easier
  (tested on RStudio version 2023.06.1).

If you have used the default number of protons, as well as output modes and paths in `sim.cfg`, you can simply
run `plots.R` to produce heatmaps of the integrated 2d dose, the 2d dose slice, and the integrated 1d dose.
Make sure the R working directory is set to the project root (or the location of your simulation output files,
if different).

# Python

The `Python` folder contains an equivalent Python implementation, which still requires precomputed cross
sections stored in the proect root. The Python implementation is much slower (100k replicates in about an hour
as opposed to a million in 2 minutes), but useful for prototyping. See `Python/plots_3d.py` for example usage.
