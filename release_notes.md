# Release notes

## 1.1

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

- Added ability to increase the sampling rate of sinogram and raw list-mode data (i.e. interpolate additional rows to the sinograms) to prevent aliasing artifacts

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
  
- Sinograms are automatically created during raw data load if raw data has not been explicitly selected
  - Raw data is still saved regardless of choices
  - Slightly speeds up the sinogram creation and uses less memory
  
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
  - Input can be binary 32-bit floats, DICOM images or BMP/PNG/TIFF (grayscale) images

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
  - Raw-list mode data with cylindrical normalization data should now work
  
- Implementation 1 should be slightly more numerically stable

- Enhanced error checking

- Random crashes when normalization is not applied with implementations 1 and 4 should be fixed

- Orthogonal and volume-based ray tracers can be computed faster by selecting the new `apply_accleration` variable
  - Gives about 30% faster speeds, but uses more device memory, especially on GPUs
  
- Added ability to specify the maximum ring difference with raw list-mode data

- Sinogram creation can now load data created from different file type than specified (e.g. if ROOT data is selected, but no ROOT raw list-mode data is found, but ASCII is available then the ASCII data is used)

- Enhanced the `main_PET` single reconstruction section

- Improved documentation

- Improved dynamic visualization

- Improved Interfile import and export

- Added a function to convert COO sparse row indices to CSR sparse row indices

## 1.0

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
