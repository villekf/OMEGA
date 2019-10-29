Welcome to the OMEGA wiki! This wiki is under construction and will get more documentation as the time goes on.

# Getting Started

## Installation

To install OMEGA, extract the contents of the zip-file (if you downloaded the release) and put the OMEGA-folder and all its subfolders to MATLAB/Octave path. Finally run install_mex. When using git clone, the process is the same. Installation on Octave is identical to that of MATLAB.

## Examples

GATE users should use the [`gate_main.m`](https://github.com/villekf/OMEGA/blob/master/gate_main.m) file to reconstruct GATE data. For non-GATE users, the file you should start with is [`main_nongate.m`](https://github.com/villekf/OMEGA/blob/master/main_nongate.m). For computing the forward and/or backward projections use `forward_backward_projections_example.m`. For custom (gradient-based) priors, use custom_prior_test_main.m. A more simplified main-file for GATE data (simple OSEM reconstruction) is available in `gate_main_simple.m`. Inveon users should use `Inveon_PET_main.m`.

A GATE example with GATE macros is available in exampleGATE-folder. Simply run the GATE macros as a GATE simulation (GATE material database needs to be in the same folder) and then run the gate_main_example-file to reconstruct the data. By default, ASCII data is used in the reconstruction. This example is based on both benchPET and the cylindrical PET example found from https://opengate.readthedocs.io/en/latest/defining_a_system_scanner_ct_pet_spect_optical.html#cylindricalpet

When using GATE data, all the output files of the specified format will be read in the specified folder. E.g. if you select ASCII data, all .dat-files with Coincidences in the file name will be loaded from the specified folder, with LMF all .ccs files and with ROOT all .root files.

Example MAT-files for non-GATE situation can be found from example-folder. These files are based on the above GATE-example.

Example Inveon data is available from: . The `Inveon_PET_main.m` file can be used automatically for this data.

# Help

For a short tutorial in image reconstruction in OMEGA see [Tutorial](https://github.com/villekf/OMEGA/wiki/Tutorial)

For help on using the individual main-fies or the various functions see [Function help](https://github.com/villekf/OMEGA/wiki/Function-help)

If you want to extract GATE scatter, randoms and/or trues data to MATLAB see [Extracting GATE scatter, randoms and trues data](https://github.com/villekf/OMEGA/wiki/Extracting-GATE-scatter,-randoms-and-trues-data)

For recommendations and things to watch out, see [Useful information](https://github.com/villekf/OMEGA/wiki/Useful-information)

## Contact

If you prefer using e-mail for contact, use the following address:

![Contact](https://github.com/villekf/OMEGA/blob/master/docs/contact.png)
