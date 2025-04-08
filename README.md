# OMEGA
Open-source multi-dimensional tomographic reconstruction software for MATLAB, GNU Octave and Python.

## Purpose

The purpose of OMEGA is twofold. First it is designed to allow easy, fast and efficient reconstruction of any ray-tracing based tomographic imaging. Secondly, it is intended for easy algorithmic development as it allows easy matrix-free implementation of the forward (`A * x`)  and backward (`A^T * y`) projections. While OMEGA allows the use of any ray-tracing based tomographic data, it is optimized for positron emission tomography (PET), computed tomography (CT) and single emission computed tomography (SPECT).

## Introduction

OMEGA is a software for [MATLAB](https://www.mathworks.com/), [GNU Octave](https://www.gnu.org/software/octave/) and [Python](https://www.python.org/) to reconstruct tomographic data. See [Features](https://omega-doc.readthedocs.io/en/latest/features.html) for more information on available features. See Known Issues and Limitations below for software limitations. If you wish to add your own code (e.g. reconstruction algorithm) see [Contributing code to OMEGA](https://github.com/villekf/OMEGA/wiki/Contributing-code-to-OMEGA).

Documentation for the current version is available at https://omega-doc.readthedocs.io/en/latest/index.html

The algorithms implemented so far include:
### Projector models
- Improved Siddon's ray tracer algorithm for the system matrix creation (code for regular Siddon available, but not used) [1,2]
- Orthogonal distance-based ray tracer [3]
- Volume of intersection ray tracer (THOR) [28].
- Interpolation-based ray tracer [31]
- Branchless distance-driven ray tracer [32,33]
- Rotation-based projector [34]
### Reconstruction algorithms
- Maximum Likelihood Expectation Maximization (MLEM) [4,5]
- Ordered Subsets Expectation Maximization (OSEM) [6]
- Complete-data Ordered Subsets Expectation Maximization (COSEM) [7]
- Enhanced COSEM (ECOSEM) [8]
- Accelerated COSEM (ACOSEM) [9]
- Row-Action Maximum Likelihood Algorithm (RAMLA) [10]
- Relaxed OSEM (ROSEM)
- Rescaled Block-Iterative EM (RBI) [11]
- Dynamic RAMLA (DRAMA) [12]
- Modified RAMLA (MRAMLA), aka modified BSREM-2 [13] 
- Block Sequential Regularized Expectation Maximization (BSREM) [14]
- One-step-late algorithm (OSL) [15]
- Preconditioned Krasnoselskii-Mann algorithm (PKMA) [29]
- Primal-dual hybrid gradient [35]
- Condat-Vu [36,37]
- Primal-dual Davis-Yin [38]
- FISTA [39]
- FDK [40]
- LSQR [46]
- CGLS [47]
- (OS-)SART [48,49]
- ASD-POCS [50]
- Barzilai-Borwein method [51]
### Prior models
- Quadratic prior (Gibbs prior with quadratic potential function)
- Huber prior [45]
- Median Root Prior (MRP) [16]
- L-filter (MRP-L) prior [17]
- Finite Impulse Response Median Hybrid (MRP-FMH) prior [17,18]
- Weighted mean prior [19,20]
- Total variation (TV) [21, 22, 23]
- Weighted TV [42]
- Total generalized variation (TGV) [24]
- Anisotropic diffusion (AD) Median Root Prior
- Asymmetric parallel levels sets prior (APLS) [22]
- Non-local means prior (NLM), including non-local TV [25,26,27]
- Relative difference prior [30]
- Generalized Gaussian Markov random field prior [41]
- Modified hyperbolic prior [23,43]
- Modified Lange prior [29,44]


## Installation

For additional install help, see [installation help](https://omega-doc.readthedocs.io/en/latest/installation.html).

Pre-built libraries are supplied in the [releases](https://github.com/villekf/OMEGA/releases), however, you can also manually compile everything. 

For manual compilation you're going to need a C++ compiler in order to compile the MEX-files/libraries and use this software. Visual Studio and GCC have been tested to work and are recommended depending on your platform (Visual Studio in Windows, GCC in Linux, clang should work in MacOS). Specifically, Visual Studio 2019 and 2022 have been tested to work in Windows 10 and as well as g++ 9.3 and g++ 10.5 on Ubuntu 22.04. MinGW++ also works though it is unable to compile ArrayFire OpenCL reconstructions (implementation 2) in Windows by default. Octave supports only MinGW++ in Windows and as such implementation 2 in Windows is only supported if you manually compile ArrayFire from source with MinGW (for instructions, see [here](https://github.com/villekf/OMEGA/wiki/Building-ArrayFire-with-Mingw-on-Windows)).

MinGW++ for MATLAB can be downloaded from [here](https://se.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler).

Visual Studio can be downloaded from [here](https://visualstudio.microsoft.com/). For Visual studio you'll only need "Desktop development with C++".

On Ubuntu you can install g++ with `sudo apt install build-essential`.

To install the OMEGA software, either simply extract the release/master package, download the MATLAB toolbox file from [releases](https://github.com/villekf/OMEGA/releases) (`OMEGA.-.Open-source.MATLAB.emission.tomography.software.mltbx`) or obtain the source code through git:  
`git clone https://github.com/villekf/OMEGA`
and then add the OMEGA folder and subfolders to MATLAB/Octave path (this is done automatically if you install with the mltbx-file) or /path/to/OMEGA/source/Python to PYTHONPATH. 
Finally, run `install_mex` in the source folder to build the necessary MEX-files or compile.py with Python. ROOT, OpenCL and CUDA support will be installed, if the corresponding files are found. 
Possible compilation errors can be seen with `install_mex(1)`. OpenCL include and library paths, ArrayFire path and ROOT path can also be set manually with 
`install_mex(0, OpenCL_include_path, OpenCL_lib_path, AF_PATH, ROOT_PATH)`. `OpenCL_include_path` should be the folder where `cl.h` is located, `OpenCL_lib_path` the folder where `OpenCL.lib/libOpenCL.so` (Windows/Linux) is located, 
`AF_PATH` the path to ArrayFire installation location and `ROOT_PATH` to ROOT installation location. For Python the paths can be input with "python compile.py -R /path/to/ROOT -A /path/to/arrayfire -O /path/to/OpenCL".

Certain features on Octave (such as normalization calculation) require packages io and statistics. You can install them from the Octave user interface with the following commands (io has to be installed first):

`pkg install -forge io`

`pkg install -forge statistics`

and then you need to load the statistics package:

`pkg load statistics`

Python only requires NumPy, though to load mat-files you need `pymatreader` and for multi-resolution reconstruction `scikit-image`.

In order to enable OpenCL support (implementations 2 and 3), you're going to need an OpenCL SDK/library and (for implementation 2) ArrayFire (see below). 
in Linux you can alternatively just install the OpenCL headers and library. Below examples are for Ubuntu, but the packages should exist for other distros as well.

Headers (required only when manually building):
`sudo apt-get install opencl-headers`

and then the library:  
`sudo apt-get install ocl-icd-opencl-dev`

Alternative libraries in case the above one fails:
`sudo apt-get install nvidia-opencl-dev`
or
`sudo apt-get install intel-opencl-icd`

In case the above doesn't work or you use Windows then you need to obtain an OpenCL SDK. The SDK can be any (or all) of the following: CUDA Toolkit, Intel OpenCL SDK, OCL-SDK, AMD APP SDK. 
On all cases, the OpenCL library and header files (only when manually building) need to be on your system's PATH. By default, the install_mex-file assumes that you have installed CUDA toolkit (Linux and Windows), 
AMD APP SDK v3.0 (Linux and Windows), OCL-SDK (Windows), AMD GPU Pro drivers (Linux) or Intel SDK (Linux and Windows). If you get an error message like "CL/cl.h: No such file or directory", the headers could not be found. 
You can manually add custom OpenCL paths with `install_mex(0, '/path/to/cl.h', '/path/to/libOpenCL.so')`. On Ubuntu you can use command `find / -iname cl.h 2>/dev/null` to find the required cl.h file and 
`find / -iname libOpenCL.so 2>/dev/null` to find the required library file. See `install_mex.m` for further details. `compile.py` functions similarly.

CUDA functionality requires CUDA toolkit.

**All library paths needs to be on system path when running the mex-files or otherwise the required libraries will not be found.**

Links:  
https://software.intel.com/en-us/intel-opencl  
https://developer.nvidia.com/cuda-toolkit  
https://github.com/GPUOpen-LibrariesAndSDKs/OCL-SDK/releases  

Once you have the header and library files, you need drivers/OpenCL runtimes for your device(s). If you have GPUs/APUs then simply having the vendor drivers should be enough. 
For Intel CPUs without an integrated GPU you need CPU runtimes (see the link below). 

For AMD CPUs it seems that the AMD drivers released around the summer 2018 and after no longer support CPUs so you need an older driver in order to get CPU support or use an alternative runtime. 
One possibility is to use PoCL http://portablecl.org/ and another is to try the Intel runtimes (link below).

Intel runtimes can be found here:
https://software.intel.com/en-us/articles/opencl-drivers


This software also uses ArrayFire library for the GPU/OpenCL implementation. You can find AF binaries from here:  
https://arrayfire.com/binaries
and the source code from here:  
https://github.com/arrayfire/arrayfire

Installing/building ArrayFire to the default location (`C:\Program Files\ArrayFire` in Windows, `/opt/arrayfire/` in Linux/MacOS) should cause `install_mex` and `compile.py` to automatically locate everything. 
However, in both cases you need to add the library paths to the system PATH. In Windows you will be prompted for this during the installation, for Linux you need to add `/opt/arrayfire/lib` or 
`/opt/arrayfire/lib64` (depending which exists) to the library path (e.g. `sudo ldconfig /opt/arrayfire/lib64/` or `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/arrayfire/lib64` if you don't have sudo permissions). 
Alternatively, in Linux, you can also build/install ArrayFire directly into the `/usr/local/` folder.

Using CUDA code instead of OpenCL requires the CUDA toolkit. On both cases the CUDA folder should be on the system path. `install_mex` and `compile.py` always attempts to build the CUDA code as well so no 
additional input is required from the user if all the header and library data is found. By default `install_mex` and `compile.py` looks for CUDA in `/usr/local/cuda/` in Linux. 
In Windows, CUDA location is determined from the environmental variables (PATH).

For additional install help, see [installation help](https://omega-doc.readthedocs.io/en/latest/installation.html).

Portions of version 2 were tested with the following GPUs: Nvidia Tesla A100, AMD Instinct MI100, Nvidia Tesla P100, Nvidia Tesla M40, Nvidia GeForce RTX 4060, Nvidia GeForce RTX 4090, AMD Radeon 7900 XT, Nvidia Titan RTX, Nvidia Quadro A6000 Ada, and Intel Arc A380.
All the GPUs were tested on Linux except AMD Radeon 7900 XT which was tested on Windows.

## Getting Started

For detailed installation instructions, see https://omega-doc.readthedocs.io/en/latest/installation.html

Precompiled libraries are included in [releases](https://github.com/villekf/OMEGA/releases). However, in case those do not work you can also manually compile all the necessary files. When using MATLAB or GNU Octave, run `install_mex` first. For Python, you need run `compile.py` located in /path/to/OMEGA/source/Python.

For basic usage, see https://omega-doc.readthedocs.io/en/latest/usage.html

Examples for MATLAB/GNU Octave are in main-files folder, while for Python in the aforementioned Python-folder.

## Features

See [Features](https://omega-doc.readthedocs.io/en/latest/features.html) for more information on available features

### Additional features

These features can be used as independent functions without any input needed from any other OMEGA files

- Save images (matrices) in MATLAB/Octave in NIfTI, MetaImage, Interfile, Analyze 7.5, DICOM and raw binary formats ([saveImage.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/saveImage.m))
- Import NIfTI, MetaImage, Interfile, Analyze 7.5, DICOM and raw binary formats into MATLAB/Octave ([importData.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/importData.m))
- Save images (matrices) in MATLAB/Octave in Interfile ([saveInterfile.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/saveInterfile.m)) or MetaImage formats ([saveMetaimage.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/saveMetaimage.m))
- Convert CT-attenuation coefficients into 511 keV attenuation coefficients ([attenuationCT_to_511.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/attenuationCT_to_511.m))
- (Experimental) Convert CT-attenuation coefficients directly from CT DICOM images into 511 keV attenuation coefficients ([create_atten_matrix_CT.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/m-files/create_atten_matrix_CT.m))
- Convert COO (Coordinate list) sparse matrix row indices into CSR (Compressed sparse row) indices ([coo_to_csr.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/coo_to_csr.m))
- Convert voxelized phantoms/sources into GATE compatible files ([Voxelized_phantom_handle.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/Voxelized_phantom_handle.m), [Voxelized_source_handle.m](https://github.com/villekf/OMEGA/blob/master/source/m-files/Voxelized_source_handle.m))



## System Requirements

MATLAB R2009a or later is mandatory. Following versions are guaranteed to work: 2022a and 2023b.

For Octave, any version above 5.0 should be fine. io, statistics and image packages are required for some features.

For Python 3.8 and above should work, though most likely earlier versions will work too.

C++11 compiler is required when manually compiling.

For Windows Visual Studio 2022 or 2019 is recommended with "Desktop development with C++", no other options are required. https://visualstudio.microsoft.com/

For Linux it is recommended to use g++ which usually comes bundled with the system. The version can matter only with MATLAB and it is recommended to use the one supported by your MATLAB version: https://www.mathworks.com/support/requirements/supported-compilers-linux.html

On MacOS Xcode should be used https://apps.apple.com/us/app/xcode/id497799835?mt=12.

OpenCL library is required for OpenCL functionality.

ArrayFire is required for implementation 2 (required for Python!).

For OpenCL, an OpenCL 1.2 compatible device is required. For CUDA, compute capability of 2.0 or higher is required.

The following third-party MATLAB codes are NOT required, but can be useful in certain specialized cases as they can be optionally used:  
https://www.mathworks.com/matlabcentral/fileexchange/27076-shuffle (Shuffle, used by random subset sampling)
https://www.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2 (vol3d v2, used for 3D visualization)  
https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image (Tools for NIfTI and ANALYZE image, to load/save Analyze files and also NIfTI files in absence of image processing toolbox).



## Known Issues and Limitations

### Python & MATLAB & Octave

Moving bed is not supported at the moment (needs to be step-and-shoot and the different bed positions need to be handled as separate cases). Though it should be possible to manually achieve a moving bed examination.

Only cylindrical symmetric scanners are supported inherently for PET, for other types of scanners the user has to input the detector coordinates or use index-based reconstruction.

For CT, only cone beam flat panel scanners are supported. For other types of scanners, the user has to input the detector coordinates or modify the data such that it is approximately flat panel.

### MATLAB & Octave

LMF output currently has to contain the time stamp (cannot be removed in GATE) and detector indices. The source location needs to be included if it was selected in the main-file, same goes for the scatter data. If you have any other options selected in the LMF output in GATE, then you will not get any sensible detector data. Source locations and/or scatter data can be deselected. LMF data, with different format than in GATE, are not supported.

LMF source information is a lot more unreliable than the ASCII or ROOT version. LMF support has been deprecated in version 2.0.

ROOT or ASCII data is not yet supported with GATE CT data.

ECAT PET geometry is supported only with ASCII data. ROOT data might also work (untested).

If you are experiencing crashes when using implementation 2, it might be caused by the graphics features of ArrayFire (AF). In this case I recommend renaming/removing the libforge.so files from the ArrayFire library folder (e.g. `/opt/arrayfire/lib64/`). Alternatively you can install the no-gl AF:  
http://arrayfire.s3.amazonaws.com/index.html (3.6.2 is the latest). Finally, you can also simply build AF from source, preferably without building Forge. This seems to apply only to Linux and affects both MATLAB and Octave. Python is unaffected!

Implementation 3 doesn't support TOF data. In general, implementation 3 is not recommended anymore and will probably be deprecated in a future release.

### MATLAB

If you get GLIBCXX_3.4.XX/CXXABI_1.3.XX not found error (or similar with a different version number) or an error about "undefined reference to dlopen/dlclose/dlsomethingelse" when building or running files, this should be fixed with one of the methods presented here:  
https://se.mathworks.com/matlabcentral/answers/329796-issue-with-libstdc-so-6

Or see the solutions in [installation help](https://omega-doc.readthedocs.io/en/latest/installation.html#linux).

If you are using ROOT data with ROOT 6.16.00 or newer you might receive the following error message: "undefined symbol: _ZN3tbb10interface78internal20isolate_within_arenaERNS1_13delegate_baseEl". This is caused by the `libtbb.so.2` used by MATLAB (located in `/matlabroot/bin/glnxa64`). Same solutions apply as with the above case (e.g. renaming the file). See [installation help](https://omega-doc.readthedocs.io/en/latest/installation.html#linux) for details.

ROOT data import is unstable in MATLAB R2018b and earlier versions due to a library incompatibility between the Java virtual machine in MATLAB and ROOT. in Linux you will experience MATLAB crashes when importing ROOT data. There is a workaround for this by using MATLAB in the no Java mode (e.g `matlab -nojvm`), though you won't have any GUI or graphic features. MATLAB R2019a and up are unaffected. It is recommended to use `nojvm` for data load only (set `options.only_sinos = true` to load only the data). The new desktop might not have this issue, but this is currently untested.

### Octave

When using Windows, implementation 2 (ArrayFire matrix free OpenCL) can only be enabled by manually building ArrayFire with Mingw. Instructions are provided [here](https://github.com/villekf/OMEGA/wiki/Building-ArrayFire-with-Mingw-on-Windows). Note that CUDA won't work even with manual building.

Implementations 3 and 5 fail to build with Octave 9.2 in Windows. Octave 8.3 should work fine.

Almost all MATLAB-based code runs significantly slower compared to MATLAB (this is due to the slowness of loops in Octave). Reconstructions are unaffected.

MAT-files that are over 2 GB are not supported by Octave and such large data sets cannot be saved in Octave at the moment.

### Python

Only implementation 2 is supported.

Status messages, such as the current iteration number, might be displayed only after the computation is already done.

### Intel

Intel GPUs do not support forward and/or backward projection masks. 


## Upcoming Features


Here is a list of features that should appear in future releases:

- Additional SPECT features
- Additional CT features
- PET scatter correction based on SSS
- Improved dual-layer PET support


## Reporting Bugs and Feature Requests

For any bug reports I recommend posting an issue on GitHub. For proper analysis I need the main-file that you have used and if you have used GATE data then also the macros. Preferably also all possible .mat files created, especially if the problem occurs in the reconstruction phase.

For feature requests, post an issue on GitHub. I do not guarantee that a specific feature will be added in the future.


## Citations

If you wish to use this software in your work, cite this paper: V-V Wettenhovi et al 2021 Phys. Med. Biol. 66 065010. The peer reviewed (open access) paper on OMEGA can be found from https://doi.org/10.1088/1361-6560/abe65f.

If you use some specific algorithm or prior, please cite one of references here or some other original paper!


## Acknowledgments

Original versions of COSEM, ACOSEM, ECOSEM, RAMLA, MRAMLA, MRP, L-filter, FMH, weighted mean, quadratic prior, sinogram coordinate and sinogram creation MATLAB codes were written by Samuli Summala. Normalization coefficient and variance reduction codes were written by Anssi Manninen. Initial work on TOF was done by Jonna Kangasniemi. Initial work on SPECT was done by Matti Kortelainen and Akuroma George. First version of Volume3Dviewer was done by Nargiza Djurabekova. The ray tracer projectors for SPECT were implemented by [Niilo Saarlemo](https://github.com/saarlemo). All other codes were written by Ville-Veikko Wettenhovi. Some pieces of code were copied from various websites (Stack Overflow, MATLAB Answers), the original sources of these codes can be found in the source files.

This work was supported by a grant from [Jane and Aatos Erkko foundation](https://jaes.fi/en/), [Instrumentarium Science Foundation](http://instrufoundation.fi/en.php), 
[Jenny and Antti Wihuri Foundation](https://wihurinrahasto.fi/?lang=en) and [The Finnish Research Impact Foundation](https://www.vaikuttavuussaatio.fi/en/). 
This work has been supported by [University of Eastern Finland](https://www.uef.fi/en) and Academy of Finland. 
This work was supported by the Research Council of Finland ([Flagship of Advanced Mathematics for Sensing Imaging and Modelling](https://fameflagship.fi/) grant 358944).


## References

1. Siddon, R. L. (1985), Fast calculation of the exact radiological path for a three dimensional CT array. Med. Phys., 12: 252-255. https://doi.org/10.1118/1.595715

2. Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path through a Pixel or Voxel Space. Journal of computing and information technology, 6 (1), 89-94.

3. Aguiar, P. , Rafecas, M. , Ortuño, J. E., Kontaxakis, G. , Santos, A. , Pavía, J. and Ros, D. (2010), Geometrical and Monte Carlo projectors in 3D PET reconstruction. Med. Phys., 37: 5691-5702. https://doi.org/10.1118/1.3501884

4. Dempster, A., Laird, N., & Rubin, D. (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. Journal of the Royal Statistical Society. Series B (Methodological), 39(1), 1-38. http://www.jstor.org/stable/2984875

5. L. A. Shepp and Y. Vardi, "Maximum Likelihood Reconstruction for Emission Tomography," IEEE Transactions on Medical Imaging, vol. 1, no. 2, pp. 113-122, Oct. 1982. https://doi.org/0.1109/TMI.1982.4307558

6. H. M. Hudson and R. S. Larkin, "Accelerated image reconstruction using ordered subsets of projection data," IEEE Transactions on Medical Imaging, vol. 13, no. 4, pp. 601-609, Dec. 1994. https://doi.org/10.1109/42.363108

7. Ing-Tsung Hsiao, Ing-Tsung Hsiao, Anand Rangarajan, Anand Rangarajan, Gene R. Gindi, Gene R. Gindi. "Provably convergent OSEM-like reconstruction algorithm for emission tomography", Proc. SPIE 4684, Medical Imaging 2002: Image Processing, (9 May 2002); https://doi.org/10.1117/12.467144

8. Ing-Tsung Hsiao and Anand Rangarajan and Parmeshwar Khurd and Gene Gindi. "An accelerated convergent ordered subsets algorithm for emission tomography", Physics in Medicine & Biology, vol. 49, no. 11, pp. 2145-2156, 2004, https://doi.org/10.1088/0031-9155/49/11/002

9. Ing-Tsung Hsiao and Hsuan-Ming Huang, "An accelerated ordered subsets reconstruction algorithm using an accelerating power factor for emission tomography", Physics in Medicine & Biology, vol. 55, no. 3, pp. 599-614, 2010, https://doi.org/10.1088/0031-9155/55/3/003

10. J. Browne and A. B. de Pierro, "A row-action alternative to the EM algorithm for maximizing likelihood in emission tomography," IEEE Transactions on Medical Imaging, vol. 15, no. 5, pp. 687-699, Oct. 1996. https://doi.org/10.1109/42.538946

11. C. L. Byrne, "Block-iterative methods for image reconstruction from projections," in IEEE Transactions on Image Processing, vol. 5, no. 5, pp. 792-794, May 1996. https://doi.org/10.1109/83.499919

12. Eiichi Tanaka and Hiroyuki Kudo, "Subset-dependent relaxation in block-iterative algorithms for image reconstruction in emission tomography," 2003 Phys. Med. Biol. 48 1405, https://doi.org/10.1088/0031-9155/48/10/312

13. Sangtae Ahn and J. A. Fessler, "Globally convergent image reconstruction for emission tomography using relaxed ordered subsets algorithms," in IEEE Transactions on Medical Imaging, vol. 22, no. 5, pp. 613-626, May 2003. https://doi.org/10.1109/TMI.2003.812251

14. A. R. De Pierro and M. E. B. Yamagishi, "Fast EM-like methods for maximum 'a posteriori' estimates in emission tomography," IEEE Trans. Med. Imag., vol. 20, pp. 280–288, Apr. 2001, https://doi.org/10.1109/42.921477

15. P. J. Green, "Bayesian reconstructions from emission tomography data using a modified EM algorithm," IEEE Transactions on Medical Imaging, vol. 9, no. 1, pp. 84-93, March 1990. https://doi.org/10.1109/42.52985

16. Alenius, Sakari and Ruotsalainen, Ulla, "Bayesian image reconstruction for emission tomography based on median root prior", European Journal of Nuclear Medicine, 1997, vo. 24, no. 3, pp. 258-265, https://doi.org/10.1007/BF01728761

17. S. Alenius and U. Ruotsalainen, "Improving the visual quality of median root prior images in PET and SPECT   reconstruction," 2000 IEEE Nuclear Science Symposium. Conference Record (Cat. No.00CH37149), Lyon, France, 2000, pp. 15/216-15/223 vol.2. https://doi.org/10.1109/NSSMIC.2000.950105
     
18. J. Astola and P. Kuosmanen, "Fundamentals of nonlinear digital filtering," CRC Press, Boca Raton, 1997.

19. K. Lange, "Convergence of EM image reconstruction algorithms with Gibbs smoothing," in IEEE Transactions on Medical Imaging, vol. 9, no. 4, pp. 439-446, Dec. 1990. https://doi.org/10.1109/42.61759

20. S. Alenius and U. Ruotsalainen, "Generalization of median root prior reconstruction," in IEEE Transactions on Medical Imaging, vol. 21, no. 11, pp. 1413-1420, Nov. 2002. https://doi.org/10.1109/TMI.2002.806415

21. Wettenhovi, VV., Kolehmainen, V., Huttunen, J. et al., "State Estimation with Structural Priors in fMRI," J Math Imaging Vis (2018) 60: 174. https://doi.org/10.1007/s10851-017-0749-x

22. M. J. Ehrhardt et al., "PET Reconstruction With an Anatomical MRI Prior Using Parallel Level Sets," in IEEE Transactions on Medical Imaging, vol. 35, no. 9, pp. 2189-2199, Sept. 2016. https://doi.org/10.1109/TMI.2016

23. Lijun Lu et al, "Anatomy-guided brain PET imaging incorporating a joint prior model," 2015 Phys. Med. Biol. 60 2145, https://doi.org/10.1088/0031-9155/60/6/2145

24. Kristian Bredies, Karl Kunisch, and Thomas Pock, "Total Generalized Variation," SIAM Journal on Imaging Sciences 2010 3:3, 492-526, https://doi.org/10.1137/090769521

25. Antoni Buades, Bartomeu Coll, Jean-Michel Morel, "A review of image denoising algorithms, with a new one," SIAM Journal on Multiscale Modeling and Simulation: A SIAM Interdisciplinary Journal, 2005, 4 (2), pp.490-530, https://doi.org/10.1137/040616024

26. Xiaoqing Cao et al, "A regularized relaxed ordered subset list-mode reconstruction algorithm and its preliminary application to undersampling PET imaging," 2015 Phys. Med. Biol. 60 49, https://doi.org/10.1088/0031-9155/60/1/49

27. Zhang, Hao et al. “Applications of nonlocal means algorithm in low-dose X-ray CT image processing and reconstruction: A review.” Medical physics vol. 44,3 (2017): 1168-1185. https://doi.org/10.1002/mp.12097

28. A Lougovski et al. "A volume of intersection approach for on-the-fly system matrix calculation in 3D PET image reconstruction," 2014 Phys. Med. Biol. 59 561, https://doi.org/10.1088/0031-9155/59/3/561

29. Y. Lin, C. R. Schmidtlein, Q. Li, S. Li and Y. Xu, "A Krasnoselskii-Mann Algorithm With an Improved EM Preconditioner for PET Image Reconstruction," in IEEE Transactions on Medical Imaging, vol. 38, no. 9, pp. 2114-2126, Sept. 2019, https://doi.org/10.1109/TMI.2019.2898271

30. J. Nuyts, D. Beque, P. Dupont and L. Mortelmans, "A Concave Prior Penalizing Relative Differences for Maximum-a-Posteriori Reconstruction in Emission Tomography,", in IEEE Transactions on Nuclear Science, vol. 49, no. 1, pp. 56 - 60, Aug. 2002, https://doi.org/10.1109/TNS.2002.998681

31. Xun Jia, Bin Dong, Yifei Lou, and Steve B Jiang, "GPU-based iterative cone-beam CT reconstruction using tight frame regularization", 2011 Physics in Medicine and Biology, Vol. 56, No. 13, https://doi.org/10.1088/0031-9155/56/13/004

32. Basu, Samit and De Man, Bruno, "Branchless distance driven projection and backprojection", SPIE Proceedings 2006, https://doi.org/10.1117/12.659893

33. Liu, Rui and Fu, Lin and De Man, Bruno and Yu, Hengyong, "GPU-Based Branchless Distance-Driven Projection and Backprojection", 2017 IEEE Transactions on Computational Imaging, Vol. 3, No. 4, https://doi.org/10.1109/TCI.2017.2675705

34. Di Bella, E.V.R. and Barclay, A.B. and Eisner, R.L. and Schafer, R.W., "A comparison of rotation-based methods for iterative reconstruction algorithms", 1996 IEEE Transactions on Nuclear Science, Vol. 43, No. 6, https://doi.org/10.1109/23.552756

35. Chambolle, Antonin and Pock, Thomas, "A First-Order Primal-Dual Algorithm for Convex Problems with Applications to Imaging", 2011 Journal of Mathematical Imaging and Vision, Vol. 40, No. 1, https://doi.org/10.1007/s10851-010-0251-1

36. Condat, Laurent, "A primal-dual splitting method for convex optimization involving Lipschitzian, proximable and linear composite terms", 2013 Journal of Optimization Theory and Applications, Vol. 158, No. 2, https://doi.org/10.1007/s10957-012-0245-9

37. Bằng Công Vũ, "A splitting algorithm for dual monotone inclusions involving cocoercive operators", 2011 Advances in Computational Mathematics, Vol. 38, No. 3, https://doi.org/10.1007/s10444-011-9254-8

38. Adil Salim and Laurent Condat and Konstantin Mishchenko and Peter Richtárik, "Dualize, Split, Randomize: Toward Fast Nonsmooth Optimization Algorithms", 2022 Journal of Optimization Theory and Applications, Vol. 195, No. 1, https://doi.org/10.1007/s10957-022-02061-8

39. Amir Beck and Marc Teboulle, "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems", 2009 SIAM Journal on Imaging Sciences, Vol. 2, No. 1, https://doi.org/10.1137/080716542

40. L. A. Feldkamp and L. C. Davis and J. W. Kress, "Practical cone-beam algorithm", 1984 J. Opt. Soc. Am. A, Vol. 1, No. 6, https://doi.org/10.1364/JOSAA.1.000612

41. Jean-Baptiste Thibault and Ken D. Sauer and Charles A. Bouman and Jiang Hsieh, "A three-dimensional statistical approach to improved image quality for multislice helical CT", 2007 Medical Physics, Vol. 34, No. 11, https://doi.org/10.1118/1.2789499

42. Liu Y, Ma J, Fan Y, Liang Z. "Adaptive-weighted total variation minimization for sparse data toward low-dose x-ray computed tomography image reconstruction," Phys Med Biol. 2012 Dec 7;57(23):7923-56. https://doi.org/10.1088/0031-9155/57/23/7923

43. Charbonnier P, Blanc-Feraud L, Aubert G, Barlaud M. "Deterministic edge-preserving regularization in computed imaging," IEEE Trans Image Process. 1997;6(2):298-311. https://doi.org/10.1109/83.551699

44. K. Lange, "Convergence of EM image reconstruction algorithms with Gibbs smoothing," in IEEE Transactions on Medical Imaging, vol. 9, no. 4, pp. 439-446, Dec. 1990, https://doi.org/10.1109/42.61759

45. Huber, Peter J. and Ronchetti, Elvezio M. "Robust Statistics," 2009, John Wiley & Sons, https://doi.org/10.1002/9780470434697

46. Christopher C. Paige and Michael A. Saunders. 1982. "LSQR: An Algorithm for Sparse Linear Equations and Sparse Least Squares," ACM Trans. Math. Softw. 8, 1 (March 1982), 43–71. https://doi.org/10.1145/355984.355989

47. Hestenes, Magnus R. and Eduard Stiefel. "Methods of conjugate gradients for solving linear systems," Journal of research of the National Bureau of Standards 49 (1952): 409-435. https://doi.org/10.6028/jres.049.044

48. A.H. Andersen, A.C. Kak, "Simultaneous Algebraic Reconstruction Technique (SART): A superior implementation of the ART algorithm," Ultrasonic Imaging, Volume 6, Issue 1, 1984, Pages 81-94, https://doi.org/10.1016/0161-7346(84)90008-7

49. G. Wang, and M. Jiang, "Ordered-Subset simultaneous algebraic reconstruction techniques (OS-SART)," Journal of X-Ray Science and Technology, 2004, 12, pp. 169-177.

50. Emil Y. Sidky and Xiaochuan Pan, "Image reconstruction in circular cone-beam computed tomography by constrained, total-variation minimization," Phys. Med. Biol., 2008, 53, 4777, http://dx.doi.org/10.1088/0031-9155/53/17/021

51. J. Barzilai and J. M. Borwein, “Two-Point Step Size Gradient Methods,” IMA Journal of Numerical Analysis, vol. 8, no. 1, pp. 141–148, 1988, https://doi.org/10.1093/imanum/8.1.141.