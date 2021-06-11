# Release notes

## OMEGA v1.2.1

### Bug fixes and enhancements

- Compilation fixes for Octave on Windows

- Implementation 2 (OpenCL) can now be used on Octave on Windows

## OMEGA v1.2.0

### New features

- Added support for transmission tomography data
  - Examples of computed tomography (CT) are included
  - Supports same implementations and algorithms as PET data
  - Implementation 1 always uses precomputation (this is done automatically)
  - Precomputation is turned off in all other implementations
  - Multi-ray Siddon and orthogonal distance-based ray tracer are not available
  
- CUDA backend now supports listmode reconstruction as well as the more efficient median root prior

- Added PKMA as a built-in algorithm for implementations 1, 2 and 4

- Implemented a total overhaul of the built-in algorithms with implementations 1, 2 and 4
  - Adding new built-in algorithms or priors is much easier now
  - This resulted in naming change for all non-MAP algorithms, regularization parameters and some relaxation parameters
  - Backwards compatibility is, however, maintained on all main-files except (possibly) custom prior
  - This has little to no effect on the normal user

- Added relative difference prior for implementations 1, 2 and 4
  - The current implementation can be memory intensive with large neighborhoods or image sizes

- Added support for 32-bit integer atomics (OpenCL only)
  - Should give about 20-30% faster computations at the expense of accuracy
  - Can cause integer overflow when the number of counts is high (in the range of thousands per voxel)
  - Numerical accuracy WILL be negatively affected
  - Should be used only when speed is of utmost importance
  
- GATE attenuation map from the MuMap actor can now be used as the attenuation image
  - Simply use the .mhd-file in `options.attenuation_datafile`
  - The resolution should be the same as the output image
  - Values are automatically scaled to 1/mm if they have some other units
  
### Bug fixes and enhancements

- Fixed the use of corrections when using the forward/backward projection class

- Improved Siddon should be about 10-20% faster with Nvidia GPUs (affects both OpenCL and CUDA)

- Orthogonal and volume-based ray tracers should be a little faster

- Fixed raw data sampling increase when using sampling larger than 2

- Fixed custom normalization data load when the data was not in the current working directory

- `ImageMed` now accepts color limits (similar to imagesc)

- Fixed MLEM with precomputed system matrix (implementation 1)
  - PSF was not correctly applied before
  - Not saving the intermediate iterations did not work before
  
- Fixed possible warnings when using quadratic prior

- Fixed a bug with quadratic prior when using implementations 1 or 4

- Fixed PSF with custom gradient-based priors when using implementation 2

- Fixed automatic image resize when using reference images with TV or APLS

- Fixed FMH weights when using 2D data

- Fixed TV prior when using 2D data

- Implementations 1 and 4 should now be faster
  - This should be especially visible when using Octave

- Fixed errors when using precomputed data with the MATLAB toolbox version

- Compilation of ROOT support in Linux and MacOS environments can now specify the ROOT installation directory
  - Previously this was possible only on Windows despite the warning messages suggesting the contrary
  
- ROOT supports now scanners with more than 65535 crystals when using R2019a or later (were previously supported only with R2018b and earlier)

- R2019a and newer can now use the (unstable) ROOT data load
  - When loading large ROOT files (about 2 GB), the data load might hang when using R2019a and newer
  - Using "legacy" ROOT data load (set with `options.legacyROOT = true`) fixes this, but causes the same instability as with earlier MATLAB versions
  - Use the "legacy" ROOT data load with `-nojvm` MATLAB option to prevent crashes
  - Octave is unaffected

- Octave fixes
  - Fixed multi-ray Siddon
  - Implementations 1 and 4 should be much faster now
  - Forward/backward projection class was not always functioning before
  
- Implementation 1 is faster when using R2021a or newer (not actually OMEGA related)

- Many MEX-files now use, when available, the new interleaved complex API
  - Has no effect on the user
  - Uses type-safe data access
  
- Reduced the number of TOF integration points with trapezoidal rule from 6 to 5

- Fixed MetaImage load when the (raw) filename has uppercase letters

- Fixed InterFile load when the (raw) file is in the working directory and is not i33 or img file
  
- CUDA backend should work more reliably now
  - CUDA may not work on Octave

## OMEGA v1.1.1

### New features

- Added `imageMed.m` function to easily visualize various views of an input 3D matrix

- Improved the memory efficiency of MRP with implementation 2 and also slightly improved the speed

- Custom detector/list-mode reconstruction is now supported by implementations 2 and 3
  - Sensitivity image can now be computed separately for all valid LORs

- Detector coordinates can now be extracted from Inveon list-mode data files to perform list-mode reconstruction

- `install_mex` will now try to use the supported g++ compiler if available instead of default (Linux only)

- Sinogram/raw data precorrection is now applied automatically before reconstructions with the selected corrections

### Bug fixes and enhancements

- Additional parameters are now saved in the `pz` cell matrix
  - These include the FWHM of the PSF, TOF information and whether gap filling was used
  
- Fixed possible crashes when using implementation 2 with MBSREM/MRAMLA

- Improved the accuracy of MBSREM/MRAMLA with implementation 1 and 2

- Fixed MBSREM estimates with implementation 2

- Fixed scatter correction with Inveon when using the scn-file(s)

- MBSREM, ROSEM and MRAMLA can now correctly use manual relaxation parameters
  - There is also now a check to make sure that the number of input relaxation parameters is at least the number of iterations
  
- Fixed a bug in `scaleImages`

- Fixed errors and crashes when loading only raw data

- `gate_main_simple` should now work without errors

- `main_PET` had erroneously set OSL-OSEM to true despite using implementation 4 and OSEM already

- LMF and Inveon support are now optional in the sense that if they cannot be built only a warning is displayed

- LMF and Inveon files should now compile on Windows with older versions of MinGW

- Fix loading of detector coordinates with simulated GATE data when using ROOT files

- Fix detector coordinates when transaxial multiplier is greater than 1

- PSF with implementation 3 is now identical with the other implementations

- Trapezoidal rule with TOF now uses 6 integration points
  - This can be easily increased by modifying a single value in the source code (see wiki on TOF for details)
  
- Changed the words "machine" to "scanner"

- PSF related functions now accept non-struct input values
  - options-struct is no longer necessary input variable
  
- Fixed RAMLA when using implementation 2

- Fixed subset_type 4 and 5
  - Both subset types were giving incorrect sinogram indices

- Beta values for MBSREM and BSREM were incorrectly flipped with AD and APLS priors (implementation 2)
  - MBSREM beta-values were used for BSREM and vice versa

- Fixed ADMRP with MLEM (implementations 2 and 4)

- Fixed OSL-COSEM with MRP (implementation 2)
  - OSL-COSEM with MRP was throwing an error previously

- Fixed OSL-MLEM when using TV prior and any of the subsequent priors (implementation 2)
  - Further estimates were not correctly updated

- Fixed OSL-MLEM with APLS (implementation 2)

- Fixed APLS when TV is not applied at the same time
  - APLS reconstructions were giving incorrect results

- Fixed RBI-OSL with APLS when using implementation 1

- Fixed BSREM for all supported implementations
  - BSREM was giving incorrect results for implementations 1 and 2
  - BSREM methods did not work when using implementation 4
  
- Fixed MRP with implementations 1 and 4 when not using FMH or L-filter at the same time

- Fixed Huber prior for ROSEM-MAP and BSREM when using implementation 4

- Fixed NLM when using BSREM with implementation 4

- Various fixes to the forward/backward projection class/functions

- Various fixes to the custom prior reconstruction
  - Algorithms other than OSEM-OSL should now work
  
- Various fixes to `mainPET.m` single reconstruction section
  
- GATE ASCII data load now supports 2D data

- Slightly improved the speed of sinogram coordinate creation
  - Slightly improves the speed of reconstructions
  
- Fixed increased sampling (interpolation) for raw data

- Fixed TV type 2 when using anatomical prior (implementation 2)

- Fixed TV type 3 when using anatomical prior (implementation 1 & 4)

- Fixed MLEM when using implementation 1
  - Fixed also errors when creating the system matrix
  
- Fixed COSEM-OSL when using implementation 1

- Fixed MLEM/OSL-MLEM with implementation 2 when using CUDA

- MATLAB documentation now has a search database
  - The MATLAB documentation search can now also search OMEGA documentation

## OMEGA v1.1.0

### New features

- Fixed non-local means regularization (NLM is no longer an experimental feature) 
  - NLM is now also available in implementation 2
  - Computations are significantly faster than before and also less memory intensive

- Added non-local total variation

- Added Huber prior

- Added support for MetaImage import/export 

- Added support for golden angle based subset sampling 

- Added zero padding support to `padding.m`

- Added experimental arc correction for sinogram data only
  - Only for non-precomputed data

- Added ability to increase the sampling rate of sinogram and raw data (i.e. interpolate additional rows to the sinograms) to prevent aliasing artifacts

- Added HTLM documentation

- Inveon data is now available from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3528056.svg)](https://doi.org/10.5281/zenodo.3528056)

- Redesigned the orthogonal distance-based ray tracer
  - 2D version is now 2.5D, i.e. it takes fully into account two dimensions (transaxial), but only partially axial dimension
  - 3D version is slower, but more accurate than the previous one (no more line artifacts)
  - Removed support for implementation 1 when precompute_lor = false
  
- Added a volume-of-intersection ray tracer (THOR)
  
- Added preliminary CUDA support for implementation 2 (run-time compilation, does not require NVCC)
  - Currently improved Siddon is faster than OpenCL, but orthogonal/volume-based is slower

- Added PSF reconstruction for all implementations and projectors
  - Optional deblurring phase possible

- Improved Siddon now supports N transaxial rays and M axial rays, previously only total ray counts of 1, 4 and 5 were allowed

- Implementation 2 no longer saves the binaries, but rather performs the compilations always on runtime (force_build has been removed)
  - Compilation times are significantly faster

- Implementation 4 now supports all algorithms except MBSREM and MRAMLA
  - All algorithms are also slightly faster

- COSEM and its variants are faster when using implementation 1 and use much less memory

- Quadratic prior and weighted mean use less memory, especially with larger neighborhoods
  
- Inveon data now supports Inveon CT UMAP-images as attenuation data

- Added support for Siemens Biograph mCT and Vision 64-bit list-mode data
  - 32-bit list-mode support for mCT as well
  - Closed source release ONLY, available for MATLAB R2015a and up on Windows and Linux
  
- Sinograms are automatically created during data load if raw data has not been explicitly selected
  - Raw data can be completely disabled by setting `options.store_raw_data = false`
  - Significantly faster sinogram creation, especially with dynamic data
  
- Added initial support for TOF data
  - Sinogram data only
  - Implementations 2, 3 and 4 only
  - Supports GATE data as well (no need to add temporal blurring)
  
- Allowed sinogram corrections to be applied without re-creating the sinogram (as long as the uncorrected sinogram already exists)

- Added a global correction factor that can be used to correct all LORs with the same correction factor (e.g. decay correction)

- Interaction coordinates (i.e. coordinate where the single has been absorbed in the crystal) can now be optionally saved with GATE data

- Custom detector coordinates are now easier to include

- Running main-files with error checking now produces much more information on the selected features
  - The user will be shown the selected data type (sinogram or raw), selected algorithm and prior, projector, corrections, device, etc.
  
- Sinogram reconstruction should now work as originally intended
  - Unintentional blurring was caused by incorrect transaxial coordinates with the oblique sinograms

- Added ability to add scatter correction to the system matrix as a diagonal matrix

- Added ability to visualize each algorithm in `visualize_pet.m` with their own color limits (set `color_from_algo = 0`)

- Added depth of interaction (DOI) support
  - This is the depth (mm) where the photon is assumed to have fully absorbed
  - Default value is 0, but can be set manually with `options.DOI`, e.g. `options.DOI = 2` assumes absorption at 2 mm depth
  - Precomputed data with different DOI will (most likely) not work
  
- Added functions to automatically create GATE compatible voxelized phantoms and/or sources
  - Supports either Interfile or MetaImage output
  - Input can be binary images, DICOM images or BMP/PNG/TIFF (grayscale) images
  
- Added MATLAB/Octave class to compute forward and backward projections
  - Includes operator overloading, e.g. computing `y = A * x` and `x = A' * y` operations with the `A` system matrix class object
  - Can compute the operations matrix-free with implementations 3 or 4, or with matrices with implementation 1
  - Can be used to extract the system matrix (or a subset of it) with implementation 1
  - See `forward_backward_projections_example.m` for examples
  
- Added support for different number of crystals in axial direction as well as transaxial multiplier
  - Both values are optional
  - Number of crystals in axial direction is required only if number of crystals in axial direction is different from the number of crystals in transaxial direction
  - Transaxial multiplier is required if blocks/modules/submodules are repeated in transaxial direction
  - Added an example that demonstrates the use of both parameters
  
- Added support for GATE submodules
  - The repeated component that includes the crystals can be either RSector/Block, Module or Submodule, i.e. both modules and submodules can be omitted
  
- It is now possible to save only scattered coincidences that have undergone N scatter events in the same medium and with the same type (e.g. Compton in phantom)
  - For example, by setting `options.scatter_components = [2 0 0 0];`, only scattered events that have undergone at least two Compton scatter events in the phantom will be stored in the scatter data
  
- Compilation on MacOS should now work with both OpenMP (as long as it is installed) and OpenCL

### Bug fixes and enhancements

- Renamed `main_nongate.m` to `main_PET.m`

- Fixed some bugs with various file import and export files 

- Fixed bugs with the prepass computations on OpenCL 

- Added mention of implementation 4 support to several algorithms and priors (e.g. supports OSL-OSEM with all priors, but only one prior at a time) 

- Implementation 2 now uses symmetric padding

- Implementation 2 is now as fast as implementation 3 when using one algorithm

- MLEM now works with implementation 2 and 3

- MLEM iterations are now updated with implementation 4

- Using smaller number of sinograms than the maximum now works with implementation 2 and 3

- Visualization with N_iter number of iterations now works without errors

- Visualization with vol3d now works without errors

- Fixed sinogram mashing

- Implemention 3 now correctly discards devices with less than 2GB of memory

- Fixed pseudo detector sinogram creation

- Visualization now supports `gate_main_simple.m` as well

- Gap filling should now work properly
  - fillmissing now uses 1D interpolation in two directions
  - Uses normalization data if available and selected

- Implementation 2 didn't work previously without randoms correction

- Orthogonal distance-based ray tracer now works correctly

- Compilation should produce less warnings when using Visual studio
  - Older versions (2013) of Visual studio should now work better (untested)

- Fixed implementation 2 when using custom prior

- Fixed bugs in normalization
  - Raw data with cylindrical normalization data should now work
  
- Implementation 1 should be slightly more numerically stable

- Enhanced error checking

- Random crashes when normalization is not applied with implementations 1 and 4 should be fixed

- Orthogonal and volume-based ray tracers can be computed faster by selecting the new `apply_accleration` variable
  - Gives about 30% faster speeds, but uses more device memory, especially on GPUs
  
- Added ability to specify the maximum ring difference with raw data

- Sinogram creation can now load data created from different file type than specified (e.g. if ROOT data is selected, but no ROOT raw list-mode data is found, but ASCII is available then the ASCII data is used)

- Enhanced the `main_PET` single reconstruction section

- Improved documentation

- Improved dynamic visualization

- Improved Interfile import and export

- Added a function to convert COO sparse row indices to CSR sparse row indices

- Converted OpenCL code to the C++ API
  - Does not require any actions from the user
  - Implementation 3 should now clear all memory after computations
  - Implementation 2 might still leave residual memory to the device that is cleared when MATLAB/Octave is closed
  
- Fixed corrections in the forward and backward projection case

- Gap filling can now be performed before the image reconstruction step without needing to separately redo the corrections

- Normalization correction should work better with mashed data

- Improved sinogram mashing

- Compilation should work more reliably in Unix platforms with less warnings

- Scatter data should now be obtained correctly

- You can now save only the very last iteration or all iterations
  - Previously all iterations were always saved
  
- Fixed asymmetric parallel level set

- Fixed the load of delayed coincidences with mCT 32-bit list-mode data

- Fixed various issues with Octave

## OMEGA v1.0.0

### New features

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
- Added proper support for pseudo detectors by allowing the user to optionally interpolate (fill) the gaps in the sinogram 

  - Two different gap filling implementations available, one utilizes a built-in function (1D interpolation) and the other a function from MATLAB file exchange (not included with OMEGA distribution) (2D interpolation) 

- Added ability to obtain only the true GATE coincidences and to reconstruct only the true coincidences 

- Added ability to obtain "true" GATE random events and scattered events 

  - Scattered events can include Compton scattering in the phantom and/or in the detector as well as Rayleigh scattering in the phantom and/or in the detector 

- Added ability to use the delayed coincidence window data as randoms correction 

- Added (optional) randoms variance reduction 

- Added (optional) randoms and/or scatter smoothing, with moving mean smoothing 

- Added normalization corrections 

  - Can be formed by using OMEGA functions or input by user (e.g. Inveon .nrm-file) 

  - Supports axial geometric, axial block, block profile, crystal efficiency and transaxial geometric corrections 

- All corrections (except attenuation) can be performed either in the sinogram phase or in the reconstruction phase 

  - Attenuation is still only in the reconstruction phase 

- Added harmonic and geometric mean to the weighted mean prior 

- Added Rescaled Block Iterative EM and its MAP version 

- Added Dynamic RAMLA (DRAMA) 

- RAMLA and BSREM now follow the original article 

- Related to the above change, relaxed OSEM (ROSEM) has been added, which was actually previously the RAMLA and BSREM algorithms 

- Several new subset selection implementations added 

  - For sinogram data only: Take every nth row or column from the sinogram, or every nth measurement on the row or column dimension 

  - For all data: Random selection (previously only for raw data), every nth measurement (previously only for sinogram data) or select subsets based on the angle of the LOR to the positive X-axis 

- Added an experimental version of the non-local means prior (implementation 1 only) 

  - Supports also patches based on an anatomical reference image and an MRP-like implementation (I.e. the median filtered image is replaced with an NLM filtered image) 

- Added a separate file to test your own priors with all the MAP algorithms and other priors 

- Added support for multi-GPU/device MLEM and OSEM (implementation 3) 

- Added separate functions to calculate forward projection and/or backprojection 

- Added new projector, the orthogonal distance ray-tracing algorithm (experimental) 

  - Supports both a 2D distance (strip/area) or 3D distance (tube/volume) 

- Scatter correction can be applied if the user has a ready-made scatter correction sinogram/raw data matrix. Scatter data created from GATE simulation in OMEGA can also be used. 

- Implementations 2,3 and 4 now support reconstruction without a precomputation phase 

  - Precomputation is still recommended for raw list-mode data 

- Added specific functions to query the number of platforms and/or devices on your system (OpenCL) 

  - ArrayFire_OpenCL_device_info() for implementation 2 (devices and their numbers) 

  - OpenCL_device_info() for Implementation 3 (platforms and their numbers and the devices on each platform) 

- Added support for multi-ray improved Siddon, with up to five rays 

- Added support for Octave 

- Added conversion function to convert OMEGA sinogram/reconstructed image data into NIfTI, Analyze, DICOM, Interfile or raw format 

- Measurement data can now be imported from NIfTI, Analyze, DICOM or raw format 

- Michlegorams with the specified span and ring difference can be visualized 

- Raw data can be visualized easily now 

- Added functions to convert CT/HU attenuation values to 511 keV attenuation 

- Added function to convert 122 keV attenuation to 511 keV 

### Bug fixes and enhancements

- Fixed a bug with LMF input and source coordinates 

- Fixed a bug with check for LMF and ROOT input 

- Fixed a bug that caused sinogram formation to fail if you had not loaded the data in the same session 

- Fixed a bug with sinogram segments that caused too low maximum ring difference values to fail 

- Fixed a bug that caused sinogram z-coordinate creation to fail with certain span and segment combinations 

  - Also fixed some bugs in the z-coordinate creation 

- Fixed a bug involving pseudo detectors 

  - The variable options.pseudot was also changed to NOT have the locations of the pseudo rings, but rather the number of pseudo rings 

- Redesigned (and fixed) the detector coordinate code 

- Fixed a bug involving perpendicular detectors in improved Siddon (all versions) 

- Fixed some bugs in the matrix free OpenCL kernel 

- Fixed some bugs when compiling the cpp-files with certain compilers 

- Added missing `padding.m`

- Corrupted ASCII files now work, skipping the corrupted line (with a warning message) 

- ASCII data load is faster on MATLAB R2019a and up 

- ASCII data no longer needs the user to specify each column number, but rather simply to copy-paste the coincidence mask. 

- ACOSEM didn't previously work if one of the other COSEMs weren't selected as well 

- Wording of reconstruction methods changed to implementations 

- Attenuation correction wasn't working correctly when used with reconstruction implementations 3 or 4 (currently implementations 4 and 2, respectively) 

- Adjusting the sinogram rotation previously had no effect, currently it instead performs the rotation in the detector space allowing the user to rotate the reconstructed image before it is actually reconstructed, though the precompute phase needs to be performed first 

- 2D dimensions of the sinograms was flipped 

- Previously if no GATE output files were found, the data load would silently fail 

- Related to the above change, the data load now also checks the current working folder for GATE data if the folder specified in fpath-variable had no data (or doesn't exist) 

- Switched the numbers of implementations 2 and 4, implementation 2 is now the matrix free OpenCL implementation 

- Sequential reconstruction is now implementation 4 

- Implementation 1 and 4 now have parallel computation support 

  - Implementation 1 uses C++11 threads, implementation 4 OpenMP 

  - C++11 threads support also added to precomputation phase 

- Modified the included example to have shorter coincidence window, include delayed coincidence events, scatter information and event IDs. Also changed the simulation to backtoback type (faster simulation). Additional cylinder was also added. 

- Added more comments, especially to the .cpp and .cl files 

- Default FMH weights now follow the same pattern as in the original article 

  - Furthermore, the default weights were previously incorrectly calculated if the neighborhood size (Ndx, Ndy or Ndz) was larger than 1 

- Default L-filter weights now resemble more closely Laplace distribution as originally intended, however, only grid sizes 3 (Ndx = 1, Ndy = 0, Ndz = 0), 9 (Ndx = 1, Ndy = 1, Ndz = 0) and 25 (Ndx = 2, Ndy = 2, Ndz = 0) use the same values as in previous literature 

  - L-filter also has an optional new weighting scheme where the weights are (loosely) based on Laplace distribution, but mainly vary according to distance 

- Raw list-mode data and sinogram data now use the same reconstruction and Siddon functions 

- Slight adjustments of default weights for weighted mean prior 

- MRAMLA (and MBSREM) follows the original article now 

- COSEM had an unnecessary multiplication with the scalar value present in ACOSEM 

- ACOSEM had the first power factor in the wrong place (f*H)^(1/h) instead of f^(1/h)*H 

- RAMLA, BSREM and regularized ROSEM MAP implementations now use the same estimate for the first iteration 

- MRP works now even without image processing toolbox, though image processing toolbox is still used if available 

- You can now choose to ignore the normalization step in MRP, L-filter, FMH, AD and weighted mean 

- Matrix-free OpenCL implementation now pre-compiles the binary for the corresponding device on the first run (implementation 2) 

  - Separate binaries are created for each device used and for each projector 

- Optimizations to the OpenCL kernels 

- The number of precomputed MAT-files has been reduced and they are now stored in mat-files-folder 

- Added separate precompute function for the OpenCL implementations (uses single precision numbers instead of doubles) 

- The raw list-mode data file is now always saved initially. Sinograms can be (if required) formed from this data 

- Install_mex now tries to build both the ROOT mex-files as well as the OpenCL/ArrayFire mex-files. If the building process is not successful a warning is shown, but the process is not halted. Install_mex(1) now shows the error logs if an error is encountered and end the building process on error. 

- Previously known reconstruction method 5 has been deprecated; the codes still exist (and the mex-files are built if OpenCL is found) and it can be used by selecting implementation 5, but it has not been tested and will mostly likely fail 

- 64-bit atomic operations can be used if the device used supports them (OpenCL only) 

- Implementations 2, 3 and 4 are now faster on subsequent iterations due to precomputation done on the first iteration. With OpenCL, this is only the case if there is sufficient memory available on the device 

- Previously only MATLAB 2015a and up were supported, now older versions should work as well 

- ROOT data import works without crashes on R2019a and up 
