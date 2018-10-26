# OMEGA
Open-source MATLAB Emission Tomography Software

## Introduction

OMEGA is a software for MATLAB to reconstruct data obtained with a positron emission tomography device. This software also allows to easily reconstruct ASCII, LMF or Root  data obtained from GATE simulations. See Features section below for more information on available features and Known Issues and Limitations for software limitations.

The algorithms implemented so far are:
- Improved Siddon's algorithm for the system matrix creation (code for regular Siddon available, but not used) [1,2]
- Maximum Likelihood Expectation Maximization (MLEM) [3,4]
- Ordered Subsets Expectation Maximization (OSEM) [5]
- Complete-data Ordered Subsets Expectation Maximization (COSEM) [6]
- Enhanced COSEM (ECOSEM) [7]
- Accelerated COSEM (ACOSEM) [8]
- Row-Action Maximum Likelihood Algorithm (RAMLA) [9]
- Modified RAMLA (MRAMLA), aka modified BSREM [10]
- Block Sequential Regularized Expectation Maximization (BSREM) [11]
- One-step-late algorithm (OSL) [12]
- Quadratic prior (Gibbs prior with quadratic potential function)
- Median Root Prior (MRP) [13]
- L-filter (MRP-L) prior [14]
- Finite Impulse Response Median Hybrid (MRP-FMH) prior [14,15]
- Weighted mean prior [16,17]



## Getting Started

GATE users should use the gate_main.m file to reconstruct GATE data. For non-GATE users, the file you should start with is main_nongate.m.

A GATE example with GATE macros is available in exampleGATE-folder. Simply run the GATE macros as a GATE simulation and then run the gate_main_example-file to reconstruct the data. By default ASCII data is used.

Example MAT-files for non-GATE situation can be found from example-folder. These files are based on the above GATE-example.

For more information see the [wiki](https://github.com/villekf/OMEGA/wiki) or [how_to_use.txt](docs/how_to_use.txt).



## Features

The following features are currently present:

- Reconstruct any PET sinogram/raw list-mode data
- Reconstruction with MLEM, OSEM, COSEM, ECOSEM, ACOSEM, RAMLA, MRAMLA, MRP, Quadratic prior, L-filter, FMH and weighted mean algorithms in MATLAB (no OpenCL/GPU support)
- Reconstruction with OSEM with GPU/OpenCL
- Import GATE LMF, ASCII or Root (experimental) data into MATLAB and either reconstruct them in their raw list-mode format or in the user specified sinogram format (see Known issues and limitations for LMF and Root limitations)
- Compare the reconstructed image with the actual "true" GATE source image (i.e. error analysis)
- Sequential (matrix-free) OSEM or MLEM reconstruction possible, both a pure CPU version and (experimental) OpenCL version (less memory intensive than the matrix versions)
- Include attenuation correction into the reconstruction
- Optionally allows to obtain only the system/observation matrix used in PET reconstruction
- All the data (e.g. sinograms, system matrix) can be used with your own algorithms
- Supports machines with pseudo detectors



## Installation

You're going to need C++ compiler in order to compile the MEX-files and use this software. Visual Studio and GCC have been tested to work so I recommend those depending on your platform. Specifically Visual Studio 2015 and 2017 have been tested to work as well as GCC 5.4.0 on Windows 7/10 and Ubuntu 16.04.

To install the software, either simply extract the package or obtain the source code through git:  
git clone https://github.com/villekf/OMEGA  
and then add the software folder and subfolders to MATLAB path. Finally, run install_mex in the source folder to build the necessary MEX-files without OpenCL or Root support. Use install_mex(1) to install also OpenCL support, but no Root support, install_mex(0,1) for Root support, but no OpenCL support, and lastly install_mex(1,1) for both.

In order to enable OpenCL support you're going to need an OpenCL SDK and ArrayFire (see below). The SDK can be any (or all) of the following CUDA Toolkit, Intel OpenCL SDK, OCL-SDK, AMD APP SDK. On all cases, the OpenCL library and header files need to be on your system's PATH. By default, the install_mex-file assumes that you have installed CUDA toolkit (Linux and Windows) or Intel SDK (Windows). If you get an error message like "CL/cl.h: No such file or directory", the headers could not be found. You can also add these manually to install_mex-file by adding -I/path/to/CL and -L/path/to/OpenCLlib before the .cpp file (simply replace the CUDA paths with the correct ones). On Ubuntu you can use command "find / -iname cl.h 2>/dev/null" to find the required cl.h file and "find / -iname libOpenCL.so 2>/dev/null" to find the required library file. See install_mex.m-file for further details.

Links:  
https://software.intel.com/en-us/intel-opencl  
https://developer.nvidia.com/cuda-toolkit  
https://github.com/GPUOpen-LibrariesAndSDKs/OCL-SDK/releases  

On Linux you can alternatively just install the OpenCL headers:  
sudo apt-get install opencl-headers  
and then the library:  
sudo apt-get install ocl-icd-libopencl1


Once you have the header and library files, you need drivers/OpenCL runtimes for your device(s). If you have GPUs then simply having the vendor drivers should be enough. For Intel CPUs without an integrated GPU you need CPU runtimes. For AMD CPUs it seems that the AMD drivers released around this summer (2018) no longer support CPUs so you need an older driver in order to get CPU support. Intel runtimes can be found here:  
https://software.intel.com/en-us/articles/opencl-drivers


This software also uses ArrayFire library for the GPU/OpenCL implementation. You can find AF binaries from here:  
https://arrayfire.com/download/  
and the source code from here:  
https://github.com/arrayfire/arrayfire



## System Requirements

MATLAB R2015a or later is mandatory since this software uses the repelem function that was introduced in 2015a. Older versions can be used if repelem is replaced with some other function. If you need to use MATLAB versions prior to 2015a, please either post an issue on the Github page or send me an e-mail. Following versions are guaranteed to work: 2017a, 2017b and 2018b.

Median root prior reconstruction requires Image Processing Toolbox.

C++ compiler is required.

OpenCL SDK/headers are required for OpenCL functionality.

ArrayFire is required for GPU/OpenCL reconstruction.

The following third-party MATLAB codes are NOT required, but can be useful as they can be optionally used:  
https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle (Shuffle)  
https://se.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2 (vol3d v2)  
https://github.com/stefanengblom/stenglib (FSPARSE)  



## Known Issues and Limitations

Submodules are not yet supported.

Root data is still experimental. On Windows, compiling doesn't work at the moment due to a bug in VS 2017. On Linux you might experience a crash either when the Root data is being imported or after it has been imported. Any actions on the GUI (e.g. inputting text in the command window, plotting, viewing variables) seems to cause the crash. Due to this, the Root import part is less verbose than the other methods (there is only notification when the job has been finished). If you are experiencing crashes, I also suggest not to take any actions in MATLAB until the job is finished (this includes the editor).

Raw list-mode data with non-GATE data is still experimental.

OpenCL matrix reconstruction on Nvidia GPUs is currently very slow. This is caused by a bug in ArrayFire that should be fixed in the next version. Matrix-free version is not affected by the slowdown.

GPU/OpenCL reconstruction only supports OSEM. Matrix-free reconstructions only support OSEM and MLEM.

LMF output currently has to contain the time stamp (cannot be removed in GATE) and detector indices as well as the source location if it was selected in the main-file. If you have any other options selected in the LMF output in GATE, then you will not get any sensible detector data. Source locations can be deselected.

LMF source information is a lot more unreliable than the ASCII or Root version.

Dynamic imaging/reconstruction is only partially supported. This means that it MIGHT work as there are code elements in place for it, but this has not yet been tested. ASCII output is needed in order to use the experimental dynamic support.

Only square image sizes are supported (e.g. 128x128, 64x64, 256x256, etc.). Z-direction is not affected.

Only machines with a total number of detectors of up to 65536 are supported. I.e. if you have a machine with more detectors than 65536 then nothing will work. This can be easily fixed though, if necessary, since it is simply caused by the use of 16-bit unsigned integers. Put up an issue on the Github page or send me an e-mail if you need a version with support for higher number of detectors.

Due to the same reason as above, maximum number of counts per pixel is 65535 (applies only to GATE data).

Moving gantry is not supported at the moment.

Only cylindrical symmetric devices are supported.

Attenuation correction can be applied only with attenuation images (e.g. CT images scaled to 511 keV).

Crystals have to be grouped in square blocks (the crystals need to be in e.g. 6x6, 8x8, 13x13, etc. combinations).

Random coincidences are not yet supported.

ECAT geometry is supported only with ASCII data.

OpenCL files might fail to build on Linux systems with an error message about GLIBCXX_3.4.XX not found or with an error about "undefined reference to dlopen/dlclose/dlsomethingelse". This should be fixed with one of the methods presented here:  
https://se.mathworks.com/matlabcentral/answers/329796-issue-with-libstdc-so-6

If you are experiencing crashes at the end of your computations when using OpenCL, it might be caused by the graphics features of ArrayFire. In this case I recommend installing (or building) the no-gl AF:  
http://arrayfire.s3.amazonaws.com/index.html (use the latest version available)



## Upcoming Features


Here is a list of features that should appear in future releases:

- Dynamic imaging support
- Fourier rebinning algorithm
- Filtered backprojection
- All/most algorithms from method 1 to all/most other methods
- Support for SPECT data
- TV prior


## Reporting Bugs and Feature Requests

For any bug reports I recommend posting an issue on Github. For proper analysis I need the main-file that you have used and if you have used GATE data then also the macros. Preferably also all possible .mat files created, especially if the problem occurs in the reconstruction phase.

For feature requests, post an issue on Github. I do not guarantee that a specific feature will be added in the future.


## Citations

If you wish to use this software in your work, at the moment cite the Github page. This will most likely change in the near future though so check back here later.


## Acknowledgments

Original versions of COSEM, ACOSEM, ECOSEM, RAMLA, MRAMLA, MRP, L-filter, FMH, weighted mean and quadratic prior codes were written by Samuli Summala. All other codes were written by Ville-Veikko Wettenhovi. Some pieces of code were copied from various websites, the original sources of these codes can be found in the source files.

This work was supported by a grant from Jane and Aatos Erkko foundation.


## References

1. Siddon, R. L. (1985), Fast calculation of the exact radiological path for a three dimensional CT array. Med. Phys., 12: 252-255. doi:10.1118/1.595715

2. Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path through a Pixel or Voxel Space. Journal of computing and information technology, 6 (1), 89-94.

3. Dempster, A., Laird, N., & Rubin, D. (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. Journal of the Royal Statistical Society. Series B (Methodological), 39(1), 1-38. http://www.jstor.org/stable/2984875

4. L. A. Shepp and Y. Vardi, "Maximum Likelihood Reconstruction for Emission Tomography," IEEE Transactions on Medical Imaging, vol. 1, no. 2, pp. 113-122, Oct. 1982. doi: 10.1109/TMI.1982.4307558

5. H. M. Hudson and R. S. Larkin, "Accelerated image reconstruction using ordered subsets of projection data," IEEE Transactions on Medical Imaging, vol. 13, no. 4, pp. 601-609, Dec. 1994. doi: 10.1109/42.363108

6. Ing-Tsung Hsiao, Ing-Tsung Hsiao, Anand Rangarajan, Anand Rangarajan, Gene R. Gindi, Gene R. Gindi. "Provably convergent OSEM-like reconstruction algorithm for emission tomography", Proc. SPIE 4684, Medical Imaging 2002: Image Processing, (9 May 2002); doi: 10.1117/12.467144

7. Ing-Tsung Hsiao and Anand Rangarajan and Parmeshwar Khurd and Gene Gindi. "An accelerated convergent ordered subsets algorithm for emission tomography", Physics in Medicine & Biology, vol. 49, no. 11, pp. 2145-2156, 2004.

8. Ing-Tsung Hsiao and Hsuan-Ming Huang, "An accelerated ordered subsets reconstruction algorithm using an accelerating power factor for emission tomography", Physics in Medicine & Biology, vol. 55, no. 3, pp. 599-614, 2010.

9. J. Browne and A. B. de Pierro, "A row-action alternative to the EM algorithm for maximizing likelihood in emission tomography," IEEE Transactions on Medical Imaging, vol. 15, no. 5, pp. 687-699, Oct. 1996. doi: 10.1109/42.538946

10. Sangtae Ahn and J. A. Fessler, "Globally convergent image reconstruction for emission tomography using relaxed ordered subsets algorithms," in IEEE Transactions on Medical Imaging, vol. 22, no. 5, pp. 613-626, May 2003. doi: 10.1109/TMI.2003.812251

11. A. R. De Pierro and M. E. B. Yamagishi, "Fast EM-like methods for maximum 'a posteriori' estimates in emission tomography," IEEE Trans. Med. Imag., vol. 20, pp. 280–288, Apr. 2001. 

12. P. J. Green, "Bayesian reconstructions from emission tomography data using a modified EM algorithm," IEEE Transactions on Medical Imaging, vol. 9, no. 1, pp. 84-93, March 1990. doi: 10.1109/42.52985

13. Alenius, Sakari and Ruotsalainen, Ulla, "Bayesian image reconstruction for emission tomography based on median root prior", European Journal of Nuclear Medicine, 1997, vo. 24, no. 3, pp. 258-265.

14. S. Alenius and U. Ruotsalainen, "Improving the visual quality of median root prior images in PET and SPECT   reconstruction," 2000 IEEE Nuclear Science Symposium. Conference Record (Cat. No.00CH37149), Lyon, France, 2000, pp. 15/216-15/223 vol.2. doi: 10.1109/NSSMIC.2000.950105
     
15. J. Astola and P. Kuosmanen, "Fundamentals of nonlinear digital filtering," CRC Press, Boca Raton, 1997.

16. K. Lange, "Convergence of EM image reconstruction algorithms with Gibbs smoothing," in IEEE Transactions on Medical Imaging, vol. 9, no. 4, pp. 439-446, Dec. 1990. doi: 10.1109/42.61759

17. S. Alenius and U. Ruotsalainen, "Generalization of median root prior reconstruction," in IEEE Transactions on Medical Imaging, vol. 21, no. 11, pp. 1413-1420, Nov. 2002. doi: 10.1109/TMI.2002.806415
