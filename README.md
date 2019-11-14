# OMEGA
Open-source MATLAB Emission Tomography Software

## Introduction

OMEGA is a software for MATLAB and Octave to reconstruct data obtained with a positron emission tomography device. This software also allows to easily reconstruct ASCII, LMF or ROOT data obtained from GATE simulations. See Features section below for more information on available features and Known Issues and Limitations for software limitations. If you wish to add your own code (e.g. reconstruction algorithm) see [Contributing code to OMEGA](https://github.com/villekf/OMEGA/wiki/Contributing-code-to-OMEGA).

The algorithms implemented so far are:
- Improved Siddon's algorithm for the system matrix creation (code for regular Siddon available, but not used) [1,2]
- Orthogonal distance-based ray tracer [3]
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
- Quadratic prior (Gibbs prior with quadratic potential function)
- Median Root Prior (MRP) [16]
- L-filter (MRP-L) prior [17]
- Finite Impulse Response Median Hybrid (MRP-FMH) prior [17,18]
- Weighted mean prior [19,20]
- Total variation (TV) [21, 22, 23]
- Total generalized variation (TGV) [24]
- Anisotropic diffusion (AD) Median Root Prior
- Asymmetric parallel levels sets prior (APLS) [22]
- Non-local means prior (NLM), including non-local TV [25,26,27]



## Getting Started

First you need to put the OMEGA-folder and all its subfolders to MATLAB/Octave path and then run `install_mex`. 

GATE users should use the `gate_main.m` file to reconstruct GATE data. For non-GATE users, the file you should start with is `main_nongate.m`. For computing the forward and/or backward projections use `forward_backward_projections_example.m`. For custom (gradient-based) priors, use `custom_prior_test_main.m`. A more simplified main-file for GATE data (simple OSEM reconstruction) is available in `gate_main_simple.m`.

A GATE example with GATE macros is available in exampleGATE-folder. Simply run the GATE macros as a GATE simulation (the GATE material database needs to be in the same folder as the macros) and then run the `gate_main_example.m` to load and reconstruct the data. By default, ASCII data is used for compatibility.

Example MAT-files for non-GATE situation can be found from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3522199.svg)](https://doi.org/10.5281/zenodo.3522199). These files are based on the above GATE-example. The original simulated GATE data can be found from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3526859.svg)](https://doi.org/10.5281/zenodo.3526859).

For more information see the [wiki](https://github.com/villekf/OMEGA/wiki).

Sinograms created with v0.9 need to be transposed for them to work.

## Features

The following features are currently present:

- Supports both MATLAB and Octave
- Reconstruct any PET sinogram/raw list-mode data
- Reconstruction with MLEM, OSEM, COSEM, ECOSEM, ACOSEM, RAMLA, MRAMLA, RBI, ROSEM, BSREM, MBSREM, DRAMA, MRP, Quadratic prior, L-filter, FMH, weighted mean, TV, TGV, AD, APLS, NLM algorithms in MATLAB/Octave (OpenCL support, except for NLM)
- Import GATE LMF, ASCII or ROOT data into MATLAB/Octave and either reconstruct them in their raw list-mode format or in the user specified sinogram format (see Known issues and limitations for LMF and ROOT limitations)
- Extract GATE scatter, randoms and/or trues data and optionally reconstruct it
- Compare the reconstructed image with the actual "true" GATE source image (i.e. error analysis)
- Matrix-free reconstruction possible, with a pure CPU version (OpenMP parallelization), pure OpenCL version (multidevice support, e.g. multiple GPUs or heterogenous computing) or OpenCL version utilizing ArrayFire libraries
- Include attenuation correction, normalization, scatter correction and/or randoms correction into the reconstruction (either user-made or OMEGA made data)
- Compute normalization coefficients
- Perform variance reduction on randoms data and/or smoothing on randoms/scatter data
- Perform GATE Monte Carlo scatter correction
- Perform corrections either to the measurement data (excluding attenuation) or during the reconstruction phase
- Optionally allows to obtain only the system/observation matrix used in PET reconstruction
- All the data (e.g. sinograms, system matrix) can be used with your own algorithms
- Supports machines with pseudo detectors
- Parallel and matrix free forward and back projection functions
- Ready-made function for custom gradient-based priors
- Two different projectors available, one being the improved Siddon's algorithm with up to 5 rays and the other the orthogonal distance-based ray tracer (2D or 3D mode)
- Support for dynamic (time-varying) imaging (time-series of images)
- Several different subset selection methods, including random sampling, golden angle sampling, every nth measurement, etc.
- Support for Siemens Inveon PET list-mode data



## Installation

You're going to need C++ compiler in order to compile the MEX-files and use this software. Visual Studio and GCC have been tested to work so I recommend those depending on your platform. Specifically, Visual Studio 2015, 2017 and 2019 have been tested to work on Windows 7/10 and as well as G++ 5.5, 6.4 and 7.3 on Ubuntu 16.04. MinGW++ also works though the ArrayFire OpenCL reconstructions (implementation 2) is not supported. Octave supports only MinGW++ and as such implementation 2 is not supported on Windows.

To install the software, either simply extract the package or obtain the source code through git:  
`git clone https://github.com/villekf/OMEGA`  
and then add the OMEGA folder and subfolders to MATLAB/Octave path. Finally, run `install_mex` in the source folder to build the necessary MEX-files. Both ROOT and OpenCL support will be installed, if the corresponding files are found. ROOT is, however, only supported on Unix-platforms. Possible compilation errors can be seen with `install_mex(1)`. OpenCL include and library paths, ArrayFire path and ROOT path can also be set manually with `install_mex(0, OpenCL_include_path, OpenCL_lib_path, AF_PATH, ROOT_PATH)`.

In order to enable OpenCL support, you're going to need an OpenCL SDK and (for implementation 2) ArrayFire (see below). 

On Linux you can alternatively just install the OpenCL headers (NOTE: UNTESTED):  
`sudo apt-get install opencl-headers`

and then the library:  
`sudo apt-get install ocl-icd-libopencl1`

The SDK can be any (or all) of the following CUDA Toolkit, Intel OpenCL SDK, OCL-SDK, AMD APP SDK. On all cases, the OpenCL library and header files need to be on your system's PATH. By default, the install_mex-file assumes that you have installed CUDA toolkit (Linux and Windows), AMD APP SDK v3.0 (Linux and Windows), OCL-SDK (Windows), AMD GPU Pro drivers (Linux) or Intel SDK (Linux and Windows). If you get an error message like "CL/cl.h: No such file or directory", the headers could not be found. You can also add these manually to `install_mex` by adding `-I/path/to/CL` and `-L/path/to/OpenCLlib` before the .cpp file (simply replace the CUDA paths with the correct ones). On Ubuntu you can use command `find / -iname cl.h 2>/dev/null` to find the required cl.h file and `find / -iname libOpenCL.so 2>/dev/null` to find the required library file. See `install_mex.m` for further details.

Links:  
https://software.intel.com/en-us/intel-opencl  
https://developer.nvidia.com/cuda-toolkit  
https://github.com/GPUOpen-LibrariesAndSDKs/OCL-SDK/releases  

Once you have the header and library files, you need drivers/OpenCL runtimes for your device(s). If you have GPUs/APUs then simply having the vendor drivers should be enough. For Intel CPUs without an integrated GPU you need CPU runtimes. For AMD CPUs it seems that the AMD drivers released around the summer 2018 and after no longer support CPUs so you need an older driver in order to get CPU support. Intel runtimes can be found here:  
https://software.intel.com/en-us/articles/opencl-drivers


This software also uses ArrayFire library for the GPU/OpenCL implementation. You can find AF binaries from here:  
https://arrayfire.com/download/
and the source code from here:  
https://github.com/arrayfire/arrayfire

Installing/building ArrayFire to the default location (`C:\Program Files\ArrayFire` on Windows, `/opt/arrayfire/` on Linux) should cause `install_mex` to automatically locate everything.



## System Requirements

MATLAB R2009a or later is mandatory. Following versions are guaranteed to work: 2017a, 2017b, 2018b and 2019a.

For Octave, 5.1. works. 4.4. should also work but is untested.

C++11 compiler is required.

OpenCL SDK/headers are required for OpenCL functionality.

ArrayFire is required for implementation 2.

The following third-party MATLAB codes are NOT required, but can be useful as they can be optionally used:  
https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle (Shuffle, used by random subset sampling)
https://se.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2 (vol3d v2, used for 3D visualization)  
https://github.com/stefanengblom/stenglib (FSPARSE, used when creating sparse matrices)



## Known Issues and Limitations

### MATLAB & Octave

Submodules are not supported.

Raw list-mode data with non-GATE data is still experimental.

Multi-device reconstruction only supports OSEM and MLEM.

Implementation 4 (OpenMP CPU) supports only one prior at a time and only a select number of algorithms.

LMF output currently has to contain the time stamp (cannot be removed in GATE) and detector indices. The source location needs to be include if it was selected in the main-file, same goes for the scatter data. If you have any other options selected in the LMF output in GATE, then you will not get any sensible detector data. Source locations and/or scatter data can be deselected.

LMF source information is a lot more unreliable than the ASCII or ROOT version.

Non-square images are untested.

Only machines with a total number of detectors of up to 65536 are supported. I.e. if you have a machine with more detectors than 65536 then nothing will work. This can be easily fixed though, if necessary, since it is simply caused by using 16-bit unsigned integers. Put up an issue on the GitHub page or send me an e-mail if you need a version with support for higher number of detectors.

Due to the same reason as above, maximum number of counts per pixel is 65535 (applies only to GATE data).

Moving bed is not supported at the moment (needs to be step-and-shoot and the different bed positions need to be handled as seprate cases).

Only cylindrical symmetric devices are supported.

Attenuation correction can be applied only with attenuation images (e.g. CT images scaled to 511 keV).

ECAT geometry is supported only with ASCII data. ROOT data might also work (untested).

If you get GLIBCXX_3.4.XX/CXXABI_1.3.XX not found error or an error about "undefined reference to dlopen/dlclose/dlsomethingelse" when building or running files, this should be fixed with one of the methods presented here:  
https://se.mathworks.com/matlabcentral/answers/329796-issue-with-libstdc-so-6

If you are using ROOT data with ROOT 6.18.00 or newer you might receive the following error message: "undefined symbol: _ZN3tbb10interface78internal20isolate_within_arenaERNS1_13delegate_baseEl". This is caused by the `libtbb.so.2` used by MATLAB (located in `/matlabroot/bin/glnxa64`). Same solutions apply as with the above case (e.g. renaming the file).

If you are experiencing crashes at the end of your computations when using implementation 2, it might be caused by the graphics features of ArrayFire. In this case I recommend installing (or building) the no-gl AF:  
http://arrayfire.s3.amazonaws.com/index.html (use the latest version available). This should currently only apply to older versions of ArrayFire.

### MATLAB

ROOT data import is unstable in MATLAB R2018b and older due to a library incompatibility between the Java virtual machine in MATLAB and ROOT. On Linux you will experience MATLAB crashes when importing ROOT data. There is a workaround for this by using MATLAB in the no Java mode (e.g `matlab -nojvm`), though you won't have any GUI features. MATLAB R2019a and up are unaffected.

ROOT is not supported on Windows, though it should, theoretically, work if you use ROOT with 32-bit MATLAB, but this is untested.

### Octave

Implementation 2 (ArrayFire matrix free OpenCL) is not supported on Windows due to a compiler incompatability between MinGW and ArrayFire.

Status messages usually only appear after the function has finished.

Almost all MATLAB-based code runs significantly slower compared to MATLAB (this is due to the slowness of loops in Octave). Reconstructions are unaffected.


## Upcoming Features


Here is a list of features that should appear in future releases:

- Support for SPECT data
- TOF support


## Reporting Bugs and Feature Requests

For any bug reports I recommend posting an issue on GitHub. For proper analysis I need the main-file that you have used and if you have used GATE data then also the macros. Preferably also all possible .mat files created, especially if the problem occurs in the reconstruction phase.

For feature requests, post an issue on GitHub. I do not guarantee that a specific feature will be added in the future.


## Citations

If you wish to use this software in your work, at the moment cite the GitHub page. This will most likely change soon though so check back here later.


## Acknowledgments

Original versions of COSEM, ACOSEM, ECOSEM, RAMLA, MRAMLA, MRP, L-filter, FMH, weighted mean, quadratic prior, sinogram coordinate and sinogram creation codes were written by Samuli Summala. Normalization coefficient and variance reduction codes were written by Anssi Manninen. All other codes were written by Ville-Veikko Wettenhovi. Some pieces of code were copied from various websites (Stack Overflow, MATLAB Answers), the original sources of these codes can be found in the source files.

This work was supported by a grant from Jane and Aatos Erkko foundation.


## References

1. Siddon, R. L. (1985), Fast calculation of the exact radiological path for a three dimensional CT array. Med. Phys., 12: 252-255. doi:10.1118/1.595715

2. Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path through a Pixel or Voxel Space. Journal of computing and information technology, 6 (1), 89-94.

3. Aguiar, P. , Rafecas, M. , Ortuño, J. E., Kontaxakis, G. , Santos, A. , Pavía, J. and Ros, D. (2010), Geometrical and Monte Carlo projectors in 3D PET reconstruction. Med. Phys., 37: 5691-5702. doi:10.1118/1.3501884

4. Dempster, A., Laird, N., & Rubin, D. (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. Journal of the Royal Statistical Society. Series B (Methodological), 39(1), 1-38. http://www.jstor.org/stable/2984875

5. L. A. Shepp and Y. Vardi, "Maximum Likelihood Reconstruction for Emission Tomography," IEEE Transactions on Medical Imaging, vol. 1, no. 2, pp. 113-122, Oct. 1982. doi: 10.1109/TMI.1982.4307558

6. H. M. Hudson and R. S. Larkin, "Accelerated image reconstruction using ordered subsets of projection data," IEEE Transactions on Medical Imaging, vol. 13, no. 4, pp. 601-609, Dec. 1994. doi: 10.1109/42.363108

7. Ing-Tsung Hsiao, Ing-Tsung Hsiao, Anand Rangarajan, Anand Rangarajan, Gene R. Gindi, Gene R. Gindi. "Provably convergent OSEM-like reconstruction algorithm for emission tomography", Proc. SPIE 4684, Medical Imaging 2002: Image Processing, (9 May 2002); doi: 10.1117/12.467144

8. Ing-Tsung Hsiao and Anand Rangarajan and Parmeshwar Khurd and Gene Gindi. "An accelerated convergent ordered subsets algorithm for emission tomography", Physics in Medicine & Biology, vol. 49, no. 11, pp. 2145-2156, 2004.

9. Ing-Tsung Hsiao and Hsuan-Ming Huang, "An accelerated ordered subsets reconstruction algorithm using an accelerating power factor for emission tomography", Physics in Medicine & Biology, vol. 55, no. 3, pp. 599-614, 2010.

10. J. Browne and A. B. de Pierro, "A row-action alternative to the EM algorithm for maximizing likelihood in emission tomography," IEEE Transactions on Medical Imaging, vol. 15, no. 5, pp. 687-699, Oct. 1996. doi:10.1109/42.538946

11. C. L. Byrne, "Block-iterative methods for image reconstruction from projections," in IEEE Transactions on Image Processing, vol. 5, no. 5, pp. 792-794, May 1996. doi:10.1109/83.499919

12. Eiichi Tanaka and Hiroyuki Kudo, "Subset-dependent relaxation in block-iterative algorithms for image reconstruction in emission tomography," 2003 Phys. Med. Biol. 48 1405

13. Sangtae Ahn and J. A. Fessler, "Globally convergent image reconstruction for emission tomography using relaxed ordered subsets algorithms," in IEEE Transactions on Medical Imaging, vol. 22, no. 5, pp. 613-626, May 2003. doi:10.1109/TMI.2003.812251

14. A. R. De Pierro and M. E. B. Yamagishi, "Fast EM-like methods for maximum 'a posteriori' estimates in emission tomography," IEEE Trans. Med. Imag., vol. 20, pp. 280–288, Apr. 2001. 

15. P. J. Green, "Bayesian reconstructions from emission tomography data using a modified EM algorithm," IEEE Transactions on Medical Imaging, vol. 9, no. 1, pp. 84-93, March 1990. doi: 10.1109/42.52985

16. Alenius, Sakari and Ruotsalainen, Ulla, "Bayesian image reconstruction for emission tomography based on median root prior", European Journal of Nuclear Medicine, 1997, vo. 24, no. 3, pp. 258-265.

17. S. Alenius and U. Ruotsalainen, "Improving the visual quality of median root prior images in PET and SPECT   reconstruction," 2000 IEEE Nuclear Science Symposium. Conference Record (Cat. No.00CH37149), Lyon, France, 2000, pp. 15/216-15/223 vol.2. doi:10.1109/NSSMIC.2000.950105
     
18. J. Astola and P. Kuosmanen, "Fundamentals of nonlinear digital filtering," CRC Press, Boca Raton, 1997.

19. K. Lange, "Convergence of EM image reconstruction algorithms with Gibbs smoothing," in IEEE Transactions on Medical Imaging, vol. 9, no. 4, pp. 439-446, Dec. 1990. doi:10.1109/42.61759

20. S. Alenius and U. Ruotsalainen, "Generalization of median root prior reconstruction," in IEEE Transactions on Medical Imaging, vol. 21, no. 11, pp. 1413-1420, Nov. 2002. doi:10.1109/TMI.2002.806415

21. Wettenhovi, VV., Kolehmainen, V., Huttunen, J. et al., "State Estimation with Structural Priors in fMRI," J Math Imaging Vis (2018) 60: 174. https://doi.org/10.1007/s10851-017-0749-x

22. M. J. Ehrhardt et al., "PET Reconstruction With an Anatomical MRI Prior Using Parallel Level Sets," in IEEE Transactions on Medical Imaging, vol. 35, no. 9, pp. 2189-2199, Sept. 2016. doi:10.1109/TMI.2016.

23. Lijun Lu et al, "Anatomy-guided brain PET imaging incorporating a joint prior model," 2015 Phys. Med. Biol. 60 2145.

24. Kristian Bredies, Karl Kunisch, and Thomas Pock, "Total Generalized Variation," SIAM Journal on Imaging Sciences 2010 3:3, 492-526

25. Antoni Buades, Bartomeu Coll, Jean-Michel Morel, "A review of image denoising algorithms, with a new one," SIAM Journal on Multiscale Modeling and Simulation: A SIAM Interdisciplinary Journal, 2005, 4 (2), pp.490-530

26. Xiaoqing Cao et al, "A regularized relaxed ordered subset list-mode reconstruction algorithm and its preliminary application to undersampling PET imaging," 2015 Phys. Med. Biol. 60 49

27. Zhang, Hao et al. “Applications of nonlocal means algorithm in low-dose X-ray CT image processing and reconstruction: A review.” Medical physics vol. 44,3 (2017): 1168-1185. doi:10.1002/mp.12097
