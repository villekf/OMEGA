# Release notes

## OMEGA v2.1.0

### Breaking (API) changes

- SPECT data-handling behavior changed with ray-tracing projectors (projector type 1)
  - This is now closer to CT-functionality, meaning that subsets are projection-based
  - Subset types 0, 8-11 are now supported instead of 0, 1, and 3
  - You'll need to switch the subset type to a supported one after updating
  - Should be more memory efficient than the previous version

- Changed size of `options.rayShiftsSource` and `options.rayShiftsDetector` 
  - Should be of the size `(2*options.nRays, options.nRowsD, options.nColsD, options.nProjections)`
  - Allows for individual ray shifts for each measurement

- Modified SPECT variable names to agree with CT variable names
  - Detector pixel size: `options.crXY` -> `options.dPitchX` and `options.dPitchY`

- Most new features will be restricted to implementation 2 now and in the future (for Python this is the only one available)
  - All other implementations will only be partially supported from now on
  - This means that testing will be limited, but potential issues will be fixed
  
### New features

- OMEGA can now be installed in Python through pip using `pip install omegatomo`

- Added support for orthogonal distance ray-tracing projector with SPECT
  - Voxels are weighted by a Gaussian distribution sampled at the orthogonal distance, the variance of which is determined by the parallel distance between detector and voxel center

- Projector type 1 now has built-in support for SPECT parallel-hole and pinhole collimators
  - Variables for focal length: `options.colFxy` and `options.colFz`
  - Both take the value Inf for parallel-hole collimators
  - Both take the value 0 for pinhole collimators
  - Tested with SIMIND v8.0 built-in examples

- Added support for 3D masks
  - Forward and backward projection masks can now be either 2D or 3D
  - In 2D case, the same mask is used at each slice/projection/sinogram
  - In 3D case, you can specify a unique mask for each
  - Needs to be the same size as the image/measurement data
  - Probably won't work with multi-resolution reconstruction, but might work with projection extrapolation
  - Idea is the same as before, pixels/voxels with value of 1 are included, pixels/voxels with 0 are omitted
  - Implementation 2 only!
  
- Added subset support for FISTA
  - Implementation 2 only!

- Added another FISTA acceleration scheme
  - Select with `options.FISTAType`, use either 0 or 1
  - Implementation 2 only!
  
- FISTA acceleration can now be used with any algorithm
  - Enable with `options.FISTA_acceleration`
  - Can lead to divergence!
  - Implementation 2 only!
  
- Added support for stochastic subset selection
  - Subsets are selected randomly during the reconstruction
  - implementation 2 only!
  
- Added support for SAGA algorithm, both emission and transmission likelihood
  - Implementation 2 only!
  
- Added Barzilai-Borwain algorithm
  - Implementation 2 only!
  - MATLAB/Octave only at the moment
  
- Added support for attenuation correction when using projector type 6 (rotation-based SPECT projector)

- Projector type 6 now supports different FOV sizes
  - Sinograms are resampled to match voxel size

- Custom reconstructions are now supported with projector type 1 and 2 with SPECT data in Python
  - Functionality is the same as with PET, CT, or the SPECT rotation-based projector
  - Easily create your own algorithms by using the built-in operators for forward and/or backward projections
  - Full GPU support
  
- Limited Mac Metal support for SPECT
  - Supports only the forward and backward projection operators
  - MATLAB only!
  - Can be run with integrated Mac GPUs
  
- Added support for TOF with list-mode data
  - A separate uint8 vector needs to be input (`options.TOFIndices`) that contains the indices to the TOF time windows specified by `options.TOFCenter`
  - Implementation 2 only!
  
- Added support for standalone GPU-based regularizers in Python
  - Can be used with any data, as long as the input image is a vector and in PyOpenCL, ArrayFire, CuPy or PyTorch format
  - Supports RDP, non-local regularizers and gradient-based TV
  
- Added support for parallel beam reconstruction (`options.useParallelBeam`)
  - Very slow
  - Implementation 2 only!
  
- Subset type 0 now supports projection images as well
  - First subset uses the N first projections, and so on
  
- volume3Dviewer can now display colorbar
  - volume3Dviewer can also be used in Python if your UI supports it
  
- volume3Dviewer can now be used in multimodal viewing, i.e. merge CT and PET/SPECT images

- "Large dimensionality" support now works correctly with non-local regularizers, RDP, GGMRF, TV, and hyperbolic prior
  - Enable with `options.largeDim = true`
  - Only a subset of the image is reconstructed at a time, allowing any reconstructions with practically any GPU

- Large dimensionality now also supports EM and sensitivity image based preconditioners

- Major overhaul of the verbosity
  - Verbose level 1 stays the same
  - Verbose level 2 now has less messages but has more timing information, i.e. the time spent per sub-iteration and the approximate time left
  - Verbose level 3 contains some messages that were previously shown with level 2. Also, verbosity level 3 gives even more timing information, down to kernel and step level
  - I.e. with verbosity level of 3, you can now see, for example, the time spent in forward and backward projections, as well as in other steps

- Multi-resolution reconstruction now keeps the voxel size fixed in the multi-resolution volumes and adjusts the FOV accordingly, previously it was the other way around

- Measurement data can now be used in unsigned 16-bit or 8-bit integer format without needing to convert it first to single precision
  - Implementation 2 only!
  - Decreases memory usage when uint16 or uint8 data is used
  - Can be especially useful with TOF data
  
- SART supports regularization now

- Hybrid projectors are now available in broader capacity
  - Previously combinations such as 43 were unavailable, but now they work
  - BDD still only supports combinations with 1 or 4
  - Rotation-based projector doesn't support hybrid methods
  
- Improved dynamic reconstruction support
  - Should now work with any data
  
- The support for gaps between rings is now better
  - Now a specific length of the gap can be input instead of the number of gaps
  - The gap lengths need to be input for each gap
  
- Added randoms/scatter smoothing and randoms/scatter variance reduction in Python
  - Also added arc correction, but it is not particularly recommended
  
- Added an example Python script for GATE 10 PET simulation reconstruction
  
### Bug fixes and enhancements

- OpenCL devices that don't support images should now work better

- BSREM has been fixed

- Lambda and alpha values are now correctly computed if left zero when only one iteration is used

- Attenuation correction now attempts to scale the data if it's of different resolution than the reconstructed image

- "Reconstruction complete in..." now correctly includes precalculation time such as sensitivity image

- Functions using random variables now support fixed seeds with `options.seed`

- TV type 1 now works correctly with anatomic weighting

- Reference images are resized to the reconstructed image if different

- APLS has been fixed

- Projector type 4 now supports backprojection mask with PET data as well

- Fixed subset types 3, 6 and 7 with projector type 4 when using PET data

- Hybrid projector 45 wasn't using the correct interpolation length before, this has been fixed

- Fixed weighting when source and/or detector are inside the FOV

- The volume of intersection ray tracer (projector type 3) was previously incorrectly weighted

- Fixed attenuation correction for projector type 1 when using SPECT data

- Fixed projector types 2 and 3

- Gradient-based TV regularization should now work with CUDA

- When using projector type 6 (SPECT), the image can now be correctly rotated when using `options.offangle`

- MBSREM now works

- Weighted TV and non-local regularizers correctly work with CUDA

- Momentum-based preconditioner now works even if the momentum values are not provided

- Quadrature-based preconditioner now works

- Corrected default weights of various regularizers in Python

- AD-MRP now correctly uses the regularization parameter

- Projector types 2 and 3 work in CUDA with listmode data when using sensitivity image computation

- PDHG and PKMA now work with CUDA

- Subset type 4 now supports column sizes that are not divisible by the number of subsets

- GGMRF weights now work when Ndx/Ndy and Ndz are different

- MetaImage load should no longer throw warnings or errors

- Static reconstruction is now correctly displayed if partitions is empty

- If multiple algorithms are selected, MATLAB/Octave now correctly displays the selected algorithms

- Filtering-based preconditioner now works with TOF data

- TOF data now supports more algorithms, such as ACOSEM, and features

- Fixed non-local regularization when using CPU with implementation 2

- Fixed mask images when using sensitivity image with list-mode data

- If the multi-resolution volumes were previously extracted (by setting CELL true), but no multi-resolution volumes were present an error was thrown. Now the volumes are extracted only if they exist

- ASD-POCS should now work in all cases

- Fixed OSL-COSEM and OSL-ACOSEM when using implementation 2

- Fixed PSF blurring with ECOSEM when using implementation 2

- Fixed arc correction and increasing the sampling rate of sinograms
  - Arc correction is not particularly recommended though
  
- All algorithms, except BB, now work with Python

- Improved error checking

- Updated all the examples

## OMEGA v2.0.0

### Breaking changes

- Support for raw data format dropped
  - If you need the raw data format, use version 1.2.1 instead
  - Similar functionality can be achieved by using "listmode" data, i.e. inputting coordinates for each measurement or using the new index-based reconstruction
  - Might be re-added in future
  - Code still largely exists, but is not used (and is untested)
  
- Precomputation removed
  - Implementation 1 uses precomputation exclusively now and it is always on
  - Other implementations now longer support precomputation phase
  - Similar effect can be achieved by inputting a mask for the forward projection
  - This also includes precomputing the observation/system matrix
  - The class object can still be used to construct an actual system matrix
  
- Support for custom prior computations removed
  - Custom prior(s) can still be implemented by using the forward/backward projector class
  - This is most conveniently achieved nowadays by using the Python version
  
- All implementations now support only one algorithm and (optional) prior combination at a time

- Several variables now have different names
  - Backwards compatibility should be preserved
  - Many variables now only have one name, for example the regularization parameter is now only beta and is used for all priors
  - Relaxation parameters are now only defined in parameter lambda (or lambdaN in Python) which is used by all algorithms using relaxation
  
- Orthogonal and volume-of-intersection based projectors are no longer supported by implementation 1

- Multi-GPU/heterogeneous computing support has been dropped for implementation 3
  - You can, however, select the platform and the device now rather than having the device selected automatically

- Changed the name of the class objects to `projectorClass`
  - Defaults to PET data
  - CT can be used by putting `options.CT = true`
  - SPECT similarly with `options.SPECT = true`
  - Using the options-struct is recommended
  
- (CB)CT reconstructions are much more efficient now than before as long as GPU computing is used

- TOF support has been dropped from implementation 3

### New features

- Added ray-based projector for SPECT reconstructions
  - Enabled by setting projector type to 1
  - Supported in MATLAB, Octave and Python
  - Supports implementations 2 and 4
  - Currently doesn't support custom reconstruction algorithms
  - Thanks to [@saarlemo](https://github.com/saarlemo)

- Added FISTA-based acceleration for every algorithm
  - There is a risk of completely failed reconstruction
  
- Overhauled relative difference prior for implementation 2
  - Like before there are two different RDP methods
  - The default is basically the original RDP, where only the left/right/top/bottom/front/back voxels are taken into account with no weighting
  - Second method, enabled with `options.RDPIncludeCorners`, on the other hand is dependent on the neighborhood size that you specify, as well as uses the same weights as quadratic prior
  - The second method can have a neigborhood of any size
  - The second method also includes an optional "reference image" weighting
  - When using other implementations, the functionality is similar to the second method but the functionality is limited (no reference image weighting) and performance can be poor
  
- Added SART and ASD-POCS
  - Latter can use any of the non-proximal priors though functionality is not guaranteed
  - Unlike the original ASD-POCS, `options.beta` affects the regularization strength as well
  
- Added an adaptive NLM weighting
  - Any non-local regularizer can use adaptive weighting when `options.NLAdaptive` is set as true
  - See the documentation for more details: https://omega-doc.readthedocs.io/en/latest/algorithms.html#nlm
  
### Bug fixes and enhancements

- Projector type 6 should now work with filtering-based preconditioner

- TOF should now work with image-based preconditioners and power method

- Fixed memory leak with (at least) OSEM and ROSEM

- Image-based preconditioners should now work with measurement-filtering

- Fixed compilation warnings on Python

- The size of the extended FOV can now be adjusted more easily with `options.eFOVLength`

- The size of the FOV can now be changed when using implementation 1

- Verbosity in certain cases, such as when using non-local priors, has been improved (MATLAB/Octave only)

- Error checking has been slightly improved

- TOF should now work with projector type 4

- Fixed compilation warnings with latest g++

- PSF blurring now works with multi-resolution

- PDHG and its variants should now work properly with projector type 6

- PSF support has been added for implementation 5

- TOF should now work with implementation 4

- PSF should now work with implementation 3

## OMEGA v2.0.0 Release Candidate

### New features
  
- Added (limited) support for Python
  - Only a subset of features are implemented
  - All core components have been implemented though
  - Supports only implementation 2 and the class object
  - Custom algorithms can be computed fully in the GPU by using ArrayFire arrays (OpenCL), PyTorch tensors (CUDA) or CuPy arrays (CUDA)
  - The OMEGA forward and/or backward projection operators thus accept ArrayFire arrays, PyTorch tensors or CuPy arrays

- Added new implementation, implementation 5
  - Essentially the same as previously computing the forward and/or backward projections with the class object
  - Computes only the forward and/or backward projections in OpenCL, rest in MATLAB/Octave
  - Does not require the installation of ArrayFire, only OpenCL
  - Supports the same features as implementations 1 and 4
  - With discrete GPU, should be faster than implementations 1 and 4, but slower than 2

- Added two new projectors for CT data
  - Projector type 4 is an interpolation-based GPU-only projector
  - Projector type 5 is the branchless distance-driven projector (GPU only)
  - Only supported by implementations 2 and 5
  
- Added one new projector for PET data
  - Projector type 4 is an interpolation-based GPU-only projector
  - Only supported by implementations 2 and 5
  
- Added support for SPECT reconstruction
  - Supports also GATE projection data
  - Supports parallel hole collimators
  - Built-in support for hexagonal and circular holes for the detector-response function
  - Forward/backward projection class object is CPU (implementation 4) only in MATLAB/Octave, OpenCL and CUDA only in Python
  - Supported only by implementations 2 and 4
  - Rotation-based projector
  
- Added support for direction-vector based reconstruction in CT
  - Added support for pitch/roll/yaw/tilt of the detector panel for cone-beam CT
  - Useful for cases where the detector panel gets (intentionally or unintentionally) rotated in all three dimensions
  
- Added support for multi-resolution reconstruction
  - Certain regions can be reconstructed with reduced accuracy (larger voxel size)
  - Recommended for CT only, but should work with PET data as well
  - Resolution can be chosen freely with `options.multiResolutionScale`, for example the default value of 1/4 means that the multi-resolution region has 4 times larger voxel size
  
- Structural/anatomical reference images can now be set as variables rather than as mat-files
  - Previously the reference image had to be stored in a mat-file
  - Mat-file support still exists

- Added subset selection based on projection images/sinograms
  - I.e. every Nth projection image/sinogram is selected
  - Alternatively the projections/sinograms can be selected randomly, by golden angle sampling or with prime factor method
  - Python supports only Nth, random and prime factor, as well as subset types 1-5

- Added support for dual-layer PET reconstruction
  - Supports both with and without detector offset
  - Automatic GATE data import (ROOT only), with crystal offset, the outermost layer should be set as a submodule and innermost as a crystal (both as crystal when without offset)
  - Triple or more layers are not supported inherently
  - Triple or more can be used with custom detector coordinates, i.e. "listmode" format (see below for index-based one)
  - Still experimental feature
  - Supports only span of 1
  - Computing of normalization coefficients are not supported
  - "Listmode" format is still highly recommended
  
- Added an index-based reconstruction method
  - Basically a "listmode" reconstruction method
  - The user inputs the transaxial detector coordinate indices (`options.trIndex`) and axial coordinate indices (`options.axIndex`)
  - The indices need to be for each measurement and two per measurement (source and detector, or detector 1 and detector 2 with PET)
  - Mainly intended for PET data, especially dual/multi-layer, but should work with other types of data too
  - The index should correspond to a coordinate, for transaxial stored in `options.x` and for axial stored in `options.z`
  - For example trIndex values [2,7] would use the `options.x` coordinate values from indices 2 and 7 (3 and 8 in MATLAB/Octave).
  - Zero-based indexing!
  - Built-in support for GATE (ROOT only!) and Inveon
  - With symmetric cases should use 66% less memory than the coordinate-based (list-mode) reconstruction
  - Does not support dynamic data yet!
  - Index numbers cannot be larger than 65535
  - Implementation 2 only!
  
- Added support for 32-bit float (single precision) computations with implementation 4
  - This is now the default setting as well
  - Slightly more memory efficient and also faster
  - Double precision can be selected by setting `options.useSingles = false` (implementation 4 only)
  - Implementations 2, 3 and 5 are single precision only, while implementation 1 is double precision only
  
- Added ROOT support for Windows
  - This uses "legacy" load that causes crashes on Linux systems with Matlab, but works on Windows
  - Requires 64-bit version of ROOT and thus at least Visual Studio 2022
  
- RDP and total variation are now faster to compute with implementation 2 and use much less memory
  
- NLM is now faster to compute with implementation 2
  
- Orthogonal distance-based and volume of intersection based projectors are faster to compute
  - Speed-ups of 6x or more are possible

- Removed MLEM as a separate algorithm
  - MLEM can still be computed by selecting OSEM with 1 subset
  - Any subset-supporting algorithm can now be run without subsets
  - MRAMLA and RAMLA are still separated from MBSREM and BSREM though
  
- Only one algorithm/prior can be selected at a time with any implementation
  - Previously implementation 2 allowed multiple different algorithms/priors
  
- Added seven different image-based preconditioners and two measurement-based preconditioners
  - For image-based these include diagonal, EM, IEM, momentum, normalized gradient, filtering, and curvature based preconditioners
  - For measurement-based, only diagonal and filtering-based preconditioners are available
  - Support is limited to only a few algorithms (MBSREM, PKMA, PDHG (and its variations), FISTA)
  
- Added several new algorithms
  - These include primal-dual hybrid gradient (PDHG), FISTA, LSQR, CGLS, primal-dual Davis-Yin (PDDY), SART
  - PDHG contains several different variations with and without regularization
  - Works with any input data (PET, CT, SPECT)
  - For CT, these use linearized data, with linearization done automatically (if the data is already linearized, put `options.usingLinearizedData = true` before reconstruction)
  
- Added power method to automatically compute the largest singular/eigenvalue
  - Can be used to automatically compute the step-size parameters in PDHG/FISTA
  - Available in both Matlab/Octave and Python
  
- Added an adaptive primal/dual step-size computation method for PDHG

- Added generalized Gaussian random Markov field (GGMRF) prior

- Added several non-local variations of previous priors
  - Non-local RDP
  - Non-local GGMRF
  - Non-local Lange prior

- Added proximal total generalized variation (TGV)
  - Replaced the previous iterative TGV
  - Only supported by PDHG and its variants, and PKMA

- Added proximal TV
  - Separate from the "gradient"-based TV
  - Enable with `options.proxTV = true`
  - Only supported by PDHG and its variants, and PKMA

- Added weighted TV
  - Enabled by using `options.TVtype = 6` along with `options.TV = true`

- Added modified hyperbolic prior
  - Implementation 2 only
  
- Renamed TVtype 4 to modified Lange prior
  - Setting `options.TV = true` and `options.TVtype = 4` uses modified Lange prior instead of TV

- Added new visualization function `volume3Dviewer`
  - Thanks to Nargiza Djurabekova for the first version
  - Functionally very similar to the MATLAB function `sliceViewer` (which requires image processing/medical imaging toolbox)
  - Works on Octave too, though scrolling doesn't update the figure immediately
  
- Added ramp-filtering with several different windowing functions
  - Window functions include Hamming, Hann, Blackman, Nuttal, Parzen, cosine, Gaussian, and Shepp-Logan
  
- Added support for experimental FDK/FBP
  - Values in the reconstructed image are not scaled optimally, especially when using PET data
  - CT, by default, uses FDK weights (`options.useFDKWeights = true`), but they can be turned off as well
  - Should be used mainly for testing purposes
  - Works with PET and CT data
  - Thanks to Hannu Siikonen for help on the FDK implementation
  
- Added the ability to use mask images to limit LORs/measurements and/or voxels to take into account during reconstructions
  - It is possible, for example, to only reconstruct a cylindrical region instead of the whole rectangular volume by inputing a cylindrical mask
  - Alternatively, it is possible to take into account only measurements from a certain region
  - Can improve computation speed
  - Masks should be uint8 2D images
  
- Added support for hybrid projectors
  - Not all combinations are tested (such as projector_type 3 and 5, i.e. 35)
  - Forward projector can thus be different from backprojection (this is already default with projector type 4 when using CT data)
  - The first value is the forward projector, while the second one is backprojection, i.e. 41 uses type 4 for forward and type 1 for backward
  - Combinations involving 1, 4 and 5 should be safe (note that 5 is CT only!)

- For non-Poisson-based algorithms, positivity can be enforced with `options.enforcePositivity = true`

- Added support for object offset
  - I.e. FOV does not need to be in the origin
  
- Combined all regularization parameters into one variable
  - Rather than have separate regularization parameters for ALL function and prior combinations, only one `options.beta` is now used
  - The change, however, is backwards compatible, i.e. you can keep using the old ones as well
  - All new algorithms, however, only use the `options.beta` value
  - Should reduce clutter in the main-files
  
- The above change was also done for relaxation parameter lambda (lambdaN in Python)

- Removed support for PKMA sigma value
  - Code still remains, although commented
  
- Implementations 2, 3 and 4 can now be used in custom C++ code
  - Can be compiled into shared/dynamic libraries
  - Python uses implementation 2 through a dynamic library
  
- Almost complete code overhaul
  - Number of code lines has been reduced significantly
  - main-files continue to function as before
  - Certain cases can be slightly slower to compute than before while others faster than before
  
- Complete rework of the example-files
  - Also added several Python examples
  
- Changed coordinate system to be centered on origin
  - Previously the built-in coordinates were always in the positive x/y/z-axis
  - Any type of coordinates can still be used when using custom coordinates and object offset
  
- Listmode reconstruction is now slightly more memory efficient

- The user can now select which iterations to save
  - Previously only all or last iteration could be saved
  - Now any iteration can be saved
  - E.g. `options.saveNIter = [9;19]` stores iterations 10 and 20 (0 is the first computed iteration, i.e. iteration 1) as well as the last one
  
- Forward projections can be now saved by setting `options.storeFP = true`
  - Stores all subiterations/iterations
  - Not affected by the setting described above
  
- Implementation 2 should be more memory efficient, i.e. use less memory

- Projector type 1 might be slightly slower than before

- In CT, the variables for detector and source offset should be clearer now

- Added support for gaps in PET rings
  - Instead of pseudo rings, you can now add manual gap after certain rings
  - The gap is always the size of one detector pitch in the axial direction
  - `options.ringGaps` should include the ring(s) that are BEFORE the gap(s)
  - E.g. `options.ringGaps = [5;10];` has gaps after rings 5 and 10
  - One-based indexing, i.e. `options.ringGaps = 1;` would refer to the first ring
  
- Added support for high-dimensional CT
  - Even dozens of gigabytes large µCT can be reconstructed on GPUs
  - Supports also very large image sizes
  - Only a subset of the data and image are sent and computed on the GPU at time
  - Depends on the number of subsets, i.e. the data is divided into NumberOfSubsets chunks
  - Limited support of features
  - Supports FDK, PDHG and PKMA only
  - Image-based preconditioners are not supported
  - Activate by setting `options.largeDim = true;`
  - Supported also in Python
  - Implementation 2 only
  
- Added support for subset-based measurement splitting
  - Only the measurement data of the current subset is transfered to the GPU
  - Can help reconstruct large datasets which don't have large image volumes but have large measurement dimensions, such as PET TOF data
  - Activate by setting `options.loadTOF = false;`
  - Supported also in Python
  - Unlike above, only affects measurement data
  - Supported by all algorithms
  - Implementation 2 only

### Bug fixes and enhancements

- Fixed normalization coefficient calculations when using DOI

## OMEGA v1.2.1

### Bug fixes and enhancements

- Fix CT projection data load when using binning values higher than 1

- Compilation fixes for Octave on Windows

- Implementation 2 (OpenCL) can now be used on Octave on Windows
  - This still requires manual building of ArrayFire with MinGW

- Fix implementation 1 and 4 in Octave when using PET data

- Some fixes for older MATLAB versions

- Fix attenuation correction when using GATE MuMap actor
  - Fix errors caused by the above fix

- Several fixes and enhancements for Voxelized_source/phantom_handles
  - The user can now manually select the row and column indices to crop
  - Lesion load when using PNG/TIFF/BMP images or DICOM data is now fixed
  - Errors regarding missing lesion data are now fixed
  - Cropping can now be disabled with source data as well
  - Both functions now output the original row/column/slice indices that are included in the cropped image
  
- SaveInterfile/MetaImage functions now accept custom pixel/voxel size

- Fix raw data load when store_raw_data was false, but use_raw_data was true

- Allow the use of more, or less, image slices than twice the number of crystal rings - 1

- Allow the use of mini blocks for GATE data
  - `options.axial_multip` can now be used to specify axial repetition with modules/submodules when R-sectors are already axially repeated
  
- Implementation 4 is now usable when verbosity is set to 0
  
- Fix crashes in certain cases, for example with certain Biograph data, when using implementation 4

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
