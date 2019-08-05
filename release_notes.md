# Release notes

## 1.0

- Fixed a bug with LMF input and source coordinates 

- Fixed a bug with check for LMF and Root input 

- Fixed a bug that caused sinogram formation to fail if you had not loaded the data in the same session 

- Fixed a bug involving pseudo detectors 

  - The variable options.pseudot was also changed to NOT have the locations of the pseudo rings, but rather the number of pseudo rings 

- Redesigned (and fixed) the detector coordinate code 

- Fixed a bug involving perpendicular detectors in improved Siddon (all versions) 

- Fixed some bugs in the matrix free OpenCL kernel 

- Fixed some bugs when compiling the cpp-files with certain compilers 

- ACOSEM didn't previously work if one of the other COSEMs weren't selected as well 

- Wording of reconstruction methods changed to implementations 

- Attenuation correction wasn't working correctly when used with reconstruction implementations 3 or 4 (currently implementations 4 and 2, respectively) 

- Adjusting the sinogram rotation previously had no effect, currently it instead performs the rotation in the detector space allowing the user to rotate the reconstructed image before it is actually reconstructed, though the precompute phase needs to be performed first 

- 2D dimensions of the sinograms was flipped 

- Previously if no GATE output files were found, the data load would silently fail 

- Related to the above change, the data load now also checks the current working folder for GATE data if the folder specified in fpath-variable had no data (or doesn't exist) 

- Added ability to view other views in visualize_pet.m (e.g. sagittal) 

- Added support for dynamic imaging 

  - Supports all input formats and reconstruction implementations  

  - Added some dynamic visualization to visualize_pet.m 

- Added support for Siemens Inveon PET list-mode and sinogram data 

  - Also added support for attenuation correction when using the machine created atn-files 

- Added an error checking code 

- Added total variation prior (TV), with optional weighting based on an anatomical reference image (e.g. CT or MR image), three different weighting schemes available 

- Added anisotropic diffusion (AD) smoothing prior (more specifically, an MRP-prior where the median filtered image is replaced with the AD smoothed image) 

- Added Asymmetric Parallel Level Sets (APLS) prior 

- Added BSREM versions of all the algorithms 

- Added all algorithms previously only available in MATLAB to OpenCL matrix free implementation 

- Switched the numbers of implementations 2 and 4, implementation 2 is now the matrix free OpenCL implementation 

- Sequential reconstruction is now implementation 4 

- Implementation 1 and 4 now have parallel computation support 

  - Implementation 1 uses C++11 threads, implementation 4 OpenMP 

  - C++11 threads support also added to precomputation phase 

- Added proper support for pseudo detectors by allowing the user to optionally interpolate (fill) the gaps in the sinogram 

  - Two different gap filling implementations available, one utilizes a built-in function (1D interpolation) and the other a function from MATLAB file exchange (not included with OMEGA distribution) (2D interpolation) 

- Added ability to obtain only the true GATE coincidences and to reconstruct only the true coincidences 

- Added ability to obtain "true" GATE random events and scattered events 

  - Scattered events can include Compton scattering in the phantom and/or in the detector as well as Rayleigh scattering in the phantom and/or in the detector 

- Added ability to use the delayed coincidence window data as randoms correction 

- Added (optional) randoms variance reduction (smoothing) 

- Added normalization corrections 

  - Can be formed by using OMEGA functions or input by user (e.g. Inveon .nrm-file) 

- All corrections (except attenuation) can be performed either in the sinogram phase or in the reconstruction phase 

  - Attenuation is still only in the reconstruction phase 

- Modified the included example to have shorter coincidence window, include delayed coincidence events, scatter information and event IDs. Also changed the simulation to backtoback type (faster simulation) 

- Added more comments, especially to the .cpp and .cl files 

- Default FMH weights now follow the same pattern as in the original article 

  - Furthermore, the default weights were previously incorrectly calculated if the neighborhood size (Ndx, Ndy or Ndz) was larger than 1 

- Default L-filter weights now resemble more closely Laplace distribution as originally intended, however, only grid sizes 3 (Ndx = 1, Ndy = 0, Ndz = 0), 9 (Ndx = 1, Ndy = 1, Ndz = 0) and 25 (Ndx = 2, Ndy = 2, Ndz = 0) use the same values as in previous literature 

  - L-filter also has an optional new weighting scheme where the weights are (loosely) based on Laplace distribution, but mainly vary according to distance 

- Raw list-mode data and sinogram data now use the same reconstruction and Siddon functions 

- Slight adjustments of default weights for weighted mean prior 

- Added harmonic and geometric mean to the weighted mean prior 

- MRAMLA (and MBSREM) follows the original article now 

- COSEM had an unnecessary multiplication with the scalar value present in ACOSEM 

- ACOSEM had the first power factor in the wrong place (f*H)^(1/h) instead of f^(1/h)*H 

- RAMLA and BSREM now follow the original article 

- Related to the above change, relaxed OSEM (ROSEM) has been added, which was actually previously the RAMLA and BSREM algorithms 

- Several new subset selection implementations added 

  - For sinogram data only: Take every nth row or column from the sinogram, or every nth measurement on the row or column dimension 

  - For all data: Random selection (previously only for raw data), every nth measurement (previously only for sinogram data) or select subsets based on the angle of the LOR to the positive X-axis 

- RAMLA, BSREM and regularized ROSEM MAP implementations now use the same estimate for the first iteration 

- Added Rescaled Block Iterative EM and its MAP version 

- Added Dynamic RAMLA (DRAMA) 

- MRP works now even without image processing toolbox, though it is still the preferred implementation 

- You can now choose to ignore the normalization step in MRP, L-filter, FMH, AD and weighted mean 

- Matrix-free OpenCL implementation now pre-compiles the binary for the corresponding device on the first run (implementation 2) 

  - Separate binaries are created for each device used 

- Optimizations to the OpenCL kernels 

- Added an experimental version of the non-local means prior (implementation 1 only) 

  - Supports also patches based on an anatomical reference image and an MRP-like implementation (I.e. the median filtered image is replaced with an NLM filtered image) 

- Added a separate file to test your own priors with all the MAP algorithms and other priors 

- Added support for multi-GPU/device MLEM and OSEM (implementation 3) 

- Added separate functions to calculate forward projection and/or backprojection 

- Added new projector, the orthogonal distance ray-tracing algorithm (experimental) 

  - Supports both a 2D distance (strip/area) or 3D distance (tube/volume) 

- The number of precomputed MAT-files has been reduced and they are now stored in mat-files-folder 

- Added separate precompute function for the OpenCL implementations (uses single precision numbers instead of doubles) 

- The raw list-mode data file is now always saved initially. Sinograms can be (if required) formed from this data 

- Scatter correction can be applied if the user has a ready-made scatter correction sinogram/raw data matrix 

- Implementations 2,3 and 4 now support reconstruction without a precomputation phase 

  - Precomputation is still recommended for raw list-mode data 

- Install_mex now tries to build both the ROOT mex-files as well as the OpenCL/ArrayFire mex-files. If the building process is not successful a warning is shown, but the process is not halted. Install_mex(1) now shows the error logs if an error is encountered and end the building process on error. 

- Previously known reconstruction method 5 has been deprecated; the codes still exist (and the mex-files are built if OpenCL is found) and it can be used by selecting implementation 5, but it has not been tested and will mostly likely fail 

- 64-bit atomic operations are used if the device used supports them (OpenCL only) 

- Added specific functions to query the number of platforms and/or devices on your system (OpenCL) 

  - ArrayFire_OpenCL_device_info() for implementation 2 (devices and their numbers) 

  - OpenCL_device_info() for Implementation 3 (platforms and their numbers and the devices on each platform) 

- Implementations 2, 3 and 4 are now faster on subsequent iterations due to precomputation done on the first iteration. With OpenCL, this is only the case if there is sufficient memory available on the device 

- Previously only MATLAB 2015a and up were supported, now older versions should work as well 
