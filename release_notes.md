# Release notes

- SPECT ellipse normalization uses one shared physical region for all resolution volumes, including shifted EFOVs. Explicit radii are preserved.
- The SPECT Gaussian projector retains its collimator origin for blur calculations and uses the clipped ellipse chord for normalization.
- Fixed zero-sized peripheral buffers in MATLAB multi-resolution projection operators.

- Dynamic reconstruction improvements
  - Timeframes can now contain different amount of projection images
  - Measurement-domain correction data is now correctly indexed in dynamic reconstruction
  - Dynamic reconstruction output now uses spatial dimensions followed by saved iteration and timeframe dimensions consistently. With `options.storeMultiResolution`, each cell/list item is one volume with this same dimension order.
- The main FOV in multi-resolution reconstruction can be shifted off-center
- Multi-resolution reconstruction for implementation 5
- options.storeMultiResolution replaces the CELL macro

Metal features:
- ArrayFire support for Metal
- Projector type 4 supports texture input with Metal
- `options.useMAD` selects the Metal math mode

New examples:
- SPECT_main_DIP_PyTorch.py
  - Deep image prior regularized SPECT reconstruction example
  - Highlights OMEGA interoperability with PyTorch framework
  - Tested with Metal (M4) and CUDA (RTX 4060)
- PET_main_gateDualLayerExample.m and gate_PET_dualLayer.py
  - Dual-layer PET GATE example (ROOT import, sinogram and index-based reconstruction)
- Elekta_XVI_CBCT_main.m and Elekta_XVI_CBCT_main.py
  - Elekta XVI CBCT example using the exported projection images (```.his``` files with ```_Frames.xml```)
  - PDHG with RDP as default, FDK with Parker weights for the short scan as alternative
  - Uses the nominal XVI geometry

New features:
- Initial Metal support for PyTorch projection operators
  - SPECT and CBCT
  - See PyTorch examples
- Metal support for reconstructions with ArrayFire backend (both Python and MATLAB)
- ```options.useMAD``` controls HIP unsafe atomics in Python custom operators; notable performance improvement with backprojection (gfx90a only!)
- PET: ```options.transaxialGaps``` adds extra spacing (mm) between crystals inside a block (cryst_per_block - 1 values, e.g. between the stacks of a module)
- PET: ```options.firstBlockAngle``` sets the angle (degrees, counter-clockwise from +x, the direction of the image's first dimension (Nx), default 0) of the first detector block; applied before ```options.offangle```, which only shifts the detector numbering
- PET: Siemens Biograph Vision normalization files (the component-based ```.n```/```.ptd``` calibration file) can now be used with any sinogram size (Ndist, Nang, span, ring difference, segment table, ```ndist_side```, 760 or 798 detectors per ring) in both MATLAB (```loadNormalizationBiograph```) and Python (```omegatomo.fileio.loadNormalizationBiograph```, also selectable in the normalization file dialog)
  - The efficiency of every LOR (crystal efficiencies, crystal interference, geometric and axial effects of the span-1 plane) is summed into the bins of the OMEGA sinogram with the same binning as in ```saveSinogram```, so the span and mashing multiplicity is included; the mean LOR efficiency is one. The coefficients multiply the measured data, ```1 ./ coefficients``` is the efficiency (zero coefficients, e.g. bins with only pseudo detector LORs, give zero efficiency, not Inf)
  - The dead-time component is ignored (count rate dependent), use ```options.global_correction_factor``` for a global dead-time factor. LORs outside of the 520 radial bins of the scanner (large ```Ndist``` with 760 crystals) use the nearest edge bin and a warning is shown
  - The efficiency of every LOR is also multiplied by D^(-n), D being the transaxial distance between the two detectors of the LOR (from the OMEGA detector coordinates) and n = ```options.visionLORLengthExponent``` (default 2, 0 = off). The factor is empirical and was validated against the Siemens reconstruction of measured Vision data, where the radial RMS error decreased from 4.66 to 2.15 % (both MATLAB and Python)
  - ```quantificationFactorBiograph``` (MATLAB; ```omegatomo.fileio.quantificationFactorBiograph``` in Python) returns the factor that converts a Vision image reconstructed with this normalization to Bq/mL, computed from the PET_CALIBRATION and PET_EM_SINO ```.ptd``` files: the quantification and cross calibration factors, the decay during the frame and the branching factor of the isotope, the frame duration, the dead-time factor (from the block singles rates and the dead-time table of the normalization file) and the mean efficiency of the normalization, which ```loadNormalizationBiograph``` now returns as a third output (```returnScale=True``` in Python; the normalization coefficients are unchanged). The result multiplies the final image and uses the Siemens decay reference of the start of the frame. A wrong isotope in the Siemens header can be replaced with ```halfLife```, ```branchingFactor``` and ```frameDuration```; the usage is shown in ```Biograph_Vision_main.m``` (```quantificationFile```)
  - The factor contains an empirical constant (2.19e-3) that was fitted to the Siemens reconstructions in Bq/mL of three measured datasets (a Derenzo phantom, a torso phantom and a pig); it has only been validated for non-TOF reconstructions (the factor is the same for the projector types 1, 2 and 3) and only for the Biograph Vision
  - The old mCT code path is unchanged; the previous Vision code path and the debug file ```normalization_test.n``` written by the MATLAB loader have been removed
  - ```Biograph_Vision_main.m``` now uses the measured-data detector geometry with Vision data (```options.diameter``` 820 mm, the physical face-to-face diameter, ```options.DOI``` 8.6 mm, the fitted mean interaction depth, i.e. an effective ring diameter of 837.2 mm, and ```options.firstBlockAngle``` 5.17 degrees), with which OMEGA reconstructions of measured Vision data match the Siemens reconstructions in size and orientation (previously about 7 % too small transaxially and rotated by 5 degrees)
  - ```options.normalization_datafile``` (both MATLAB and Python) selects the user normalization file (with ```options.use_user_normalization``` in MATLAB, or when ```options.normalization``` is empty in Python) without a file dialog; if empty, the previous behaviour (default file lookup and/or dialog) is used. A missing file gives an error
  - ```Biograph_Vision_main.m``` uses no attenuation image rotation (```options.rotateAttImage = 0```) for measured Vision data with the Siemens CT images (1 is kept for the GATE data)
- Added HIP support, and widened Metal support
  - HIP support includes both implementation 2 and CuPy/PyTorch support
  - Implementation 2 requires ArrayFire to be built for HIP, for which you need at the moment this specific branch: https://github.com/arrayfire/arrayfire/tree/experimental/hip
  - CuPy and PyTorch support in custom reconstruction only projector types 1, 2, 3, and 4, and for 4 only the backprojection, i.e. 14
  - ```install_mex``` (MATLAB) and ```compile.py``` (Python) try to build HIP first (```HIP_PATH```/```ROCM_PATH```) and fall back to CUDA if ROCm is not found or the HIP build fails; OpenCL is always built
- Added support for custom sensitivity image coordinates
  - You can now add both transaxial (`options.xSens`) and axial (`options.zSens`) coordinates specifically for sensitivity image computation
  - All possible combinations are computed
  - xSens should thus be 2xM, where M is the total number of (valid) detectors, and zSens 1xM
- Custom reconstructions in Python no longer require addProjector or initProj
  - You can now ignore both and go straight to calling the operators
  - It's still recommended to at least initialize the projector with initProj
- Added Python ports of ```create_atten_matrix_CT``` and ```attenuationCT_to_511```: ```omegatomo.fileio.createAttenMatrixCT``` and ```omegatomo.fileio.attenuationCTTo511```
  - ```createAttenMatrixCT(options)``` reads a DICOM CT folder (```options.attenuation_datafile```) and returns the attenuation image in 1/mm that can be input to ```options.vaimennus```
  - The water attenuation value follows ```options.scatter_correction```: 0.096 1/cm with scatter correction and 0.093 1/cm without, as in MATLAB
  - Requires the pydicom package
- The workgroup/block size can now be specified by the user with `options.local_size`. The values can, and should be, defined for all three dimensions
- Regularization can now be applied at every Nth iteration with `options.regEveryIter`. For example `options.regEveryIter = 2` computes regularization every second iteration (the first and last iterations are always regularized)
- Add support for hyperbolic prior for CPU (implementation 2 CPU)
- Projector type 2 (orthogonal distance) now supports TOF, both with sinograms and list-mode data, also as the forward or backward projector of a combined projector type (2, 12, 21, 22, 23, 32, 24 and 42)
- Index-based (list-mode) reconstruction now works with the normalization correction (```options.normalization_correction```) in MATLAB and Python
  - The normalization can be input either as a sinogram (```Ndist * Nang * NSinos``` elements, or ```Ndist * Nang * TotSinos``` elements of which the first ```NSinos``` planes are used) or as one value for each event. A sinogram is converted to one value for each event with the same detector pair/ring pair to sinogram element mapping as in the sinogram loaders (new function ```sinogramPairLookup```); events outside the sinogram (radial FOV, ring difference) get zero normalization. Normalization files (```.nrm```, ```.n```/```.ptd``` and mat-files) are read and checked with the sinogram dimensions
  - A sinogram normalization is divided by the number of LORs in each bin (span > 1, mashing), since a sinogram bin is the sum of its LORs and each event is a single LOR. A normalization with one value for each event is used as it is, so it has to be given for a single LOR
  - The sensitivity image (```options.compute_sensitivity_image```) uses the normalization through lookup tables from the detector pair/ring pair to the sinogram element, so all detector pairs are normalized; this requires the sinogram form. LORs outside the sinogram (radial FOV, ring difference) are not counted in the sensitivity image when a sinogram normalization is used
  - Not supported (an error is raised): a sinogram normalization with multiple layers or with pseudo detectors/rings (use one value for each event), a normalization with one value for each event together with the sensitivity image (index-based data without a pair table), and, when the sensitivity image of list-mode data is computed, the normalization with custom (not index-based) detector coordinates or with user-input ```options.xSens```/```options.zSens``` without a pair table (see below), with the CPU implementation, as well as a multiplicative scatter correction/correction vector. Dynamic list-mode data still has the normalization disabled
  - New: the sensitivity image of list-mode data (index-based and custom detector coordinates, with or without ```options.xSens```/```options.zSens```) can be normalized with a pair table, ```options.normSens```: one single-precision multiplier for each ordered detector/ring pair of the sensitivity image grid, ```det_per_ring^2 * rings^2 * nLayers^2``` elements (```det_per_ring = numel(xSens)/2/nLayers``` and ```rings = numel(zSens)/nLayers``` with ```xSens```/```zSens```, the scanner values otherwise; index ```a + b*det_per_ring + (ringA*rings + ringB)*det_per_ring^2```, an array of size ```[det_per_ring, det_per_ring, rings, rings]``` indexed ```(a, b, ringB, ringA)```, layer pairs slowest). This is meant for systems that the cylindrical scanner settings cannot describe, where no sinogram exists. Only the pairs of one orientation are read (```a < b```, or ```a == b``` and ```ringA < ringB```); a zero removes a pair (e.g. both detectors on the same panel). The table must be finite and non-negative, and ```options.normalization``` must then be one value for each event (the table is not used to derive it; a sinogram normalization is an error). No kernel changes, the table is read as the per-event normalization is (no sinogram lookup); the sinogram form with lookup tables is unchanged. Not supported with dynamic data or the CPU implementation; the table is a single device buffer (4 bytes per pair), so it is limited by the maximum buffer size of the device (checked, with an error message)
- Index-based and list-mode PET: the scatter (```options.ScatterC``` with ```options.scatter_correction = true``` and ```options.subtract_scatter = true```) is now used as the additive term of each event, in MATLAB (previously ignored) and Python, with and without TOF
  - ```options.ScatterC``` is one value for each event, in the order of the events (```trIndex```/```axIndex```, or the custom coordinates) before the subsets; the unit is the measured counts expected for the LOR of the event, in the TOF bin of the event with TOF data, over the frame (the model is ```n * (A x) + s```), with no division by the normalization, the attenuation, the number of TOF bins or the number of LORs in a bin. The values must be finite and non-negative, and their number must equal the number of events (an error is raised otherwise)
  - Randoms are still not used with list-mode data (```options.randoms_correction``` is switched off): add them into the same vector. The scatter is not part of the sensitivity image. Not supported (an error is raised): dynamic data, precorrection (```corrections_during_reconstruction = false``` or algorithms without ordinary Poisson such as PDHG), largeDim and scatter smoothing. Non-subtracted (multiplicative) scatter is still rejected with TOF and with the sensitivity image, and in MATLAB it is disabled with a warning (the ```Scatter correction ON``` message is no longer printed in this case)
  - With a sinogram scatter estimate the event value has to be formed by the user: the value of the bin of the event, divided by the number of LORs in the bin with span > 1 or mashing (without the division the scatter is too large by that factor, about 12 % in the image sum in the cylindrical PET example at span 3), and with TOF the value of the TOF bin of the event, with the TOF bin of the LORs that are swapped in the sinogram binning reversed as in ```saveSinogram``` (without the reversal about 1.3 % image difference in the tested case); there is no division by the number of TOF bins for a TOF-resolved estimate
  - Validation (cylindrical PET example, 1.85 million events): the same image as a sinogram reconstruction with the same additive term to float noise (relative RMS difference 4e-6 at span 1, CUDA and OpenCL, with and without normalization and TOF, projector type 1; the projector types 2-4 show the same differences between the sinogram and index-based reconstructions as without scatter), the same result in MATLAB and Python (3e-5), and the same result for subset types 0, 1 and 3 with the events permuted by the library or by the user
  - The GPU backends now add the additive term of list-mode TOF data to the single value of each event (previously it was tiled by the number of TOF bins, which does not match list-mode data and stopped the reconstruction)
- Multi-layer PET
  - GATE ROOT and ASCII sinogram import now map every layer into a common block structure: layers with a different number of crystals (transaxially and axially) and pseudo detectors with layers are supported; the detectors/rings that a layer does not have are left empty at the end of each block
  - MATLAB GATE ASCII import now supports multi-layer data (```options.dualLayerSubmodule``` is supported as in the ROOT import)
  - Span > 1: the projector setup and the loaders now use the same number of sinograms (```TotSinos``` = sum of the segment table times the number of layer combinations), so multi-layer span > 1 data can be reconstructed end to end
  - A scalar ```options.cryst_per_block```, ```options.cryst_per_block_axial``` or ```options.crystH``` is now expanded to all layers; a length between 1 and the number of layers raises an error (previously out-of-bounds reads)
- Python sinogram creation (```omegatomo.util.saveSinogram```)
  - Trues, scatter and randoms sinograms can now be created (```trueIndex```, ```scatterIndex``` and ```randomIndex``` arguments)
- Projector type 4 now supports Joseph-style interpolation: set ```options.dL``` to a negative value to sample the ray once per voxel plane along its dominant direction (implementation 2, all GPU backends, MATLAB and Python)
- Non-local Geman-McClure prior: ```options.NLM = true``` together with ```options.NLGM = true```, scaling parameter ```options.GM_delta``` (default 1) (MATLAB and Python, all implementations)
- ```options.NLMaxWeight``` (default false) includes the voxel itself in the non-local neighbourhood with the largest neighbourhood weight, as in the original NLM; works with all non-local variants
- gpuArray inputs for the MATLAB projection operators (implementation 5, ```A * x``` and ```A' * y```)
  - gpuArray inputs are used directly on the GPU with the CUDA build of implementation 5 (```CUDA_matrixfree_multi_gpu```, built by ```install_mex``` when the CUDA toolkit is found); no extra options are needed
  - Without the CUDA build implementation 5 works with OpenCL as before, but gpuArray inputs give an error
- SPECT: the collimator hole length can differ transaxially and axially with ```options.colLxy``` and ```options.colLz``` (MATLAB and Python, ray tracing and rotation projectors); both default to ```options.colL```, so existing scripts are unchanged
- Multi-ray projection (```options.n_rays_transaxial```/```options.n_rays_axial```) now also works with CT (not recommended: slower with little benefit)
- Index-based data can now be used with the MATLAB projection operators and custom reconstructions (implementation 5)
- Helical CT
  - Projector type 5 (BDD) and the hybrid types 45 and 54 are now supported (MATLAB and Python, all GPU backends); see ```helical_CT_mainExample```
  - Helical CT now works with the Python projector operators
  - The backprojection of types 4 and 5 only processes the projections that can reach each axial block of voxels, which speeds up helical reconstructions (not used with list-mode data, subset type 12 or parallel beam)
- Python
  - The same options object can be used for several reconstructions in the same session: ```options.keepOptionsUnchanged``` (default True) makes the reconstruction work on a private copy of the options. Set it to False to save host memory, in which case the options object should not be reused
  - The projector operators accept multi-dimensional C-ordered arrays; with ```options.useRowMajor``` (True by default only with PyTorch) multi-dimensional PyTorch tensors are returned in the same C ordering without manual reordering. Flattened inputs are still assumed to be in Fortran order
  - The projector operators now support list-mode and index-based data with TOF
  - Dynamic and static attenuation data, both image- and measurement-based, now work with the operators and the reconstructions; attenuation data with the size of the image is treated as an image also when ```options.CT_attenuation = False``` (as in MATLAB)
  - L-filter (```options.L```) and FMH (```options.FMH```) priors (implementation 2)
  - Subset types 6, 7, 10 and 12, as in MATLAB
  - Sinogram gap filling (```options.fill_sinogram_gaps```) with the ```fillmissing``` method (```options.gap_filling_method = 'fillmissing'```, interpolation method from ```options.interpolation_method_fillmissing```); ```inpaint_nans``` and ```fillmissing2``` remain MATLAB only
  - Ports of the MATLAB save functions: ```omegatomo.fileio.saveImage```, ```saveInterfile```, ```saveMetaImage``` and ```scaleImages```
  - ```omegatomo.util.powerMethod``` applies the image-based preconditioners (```options.precondTypeImage```; the gradient- and curvature-based ones are not yet available in Python)
  - The filtering-based measurement preconditioner now works with CuPy custom reconstructions (see ```CBCT_custom_algorithm_example_CuPy.py```)
  - Further MATLAB parity: ```options.OSL_COSEM = 2``` (OSL-ACOSEM), ```options.global_correction_factor```, ```options.useHalf```, scatter variance reduction, a scalar ```options.lambda```, ```options.rhoPrecond```/```options.delta1Precond```, and ```options.sigma_x``` computed from ```options.TOF_FWHM``` when not given
  - The MATLAB option names ```subset_type```, ```use_CUDA```, ```use_CPU``` and ```use_device``` are accepted as aliases of ```subsetType```, ```useCUDA```, ```useCPU``` and ```deviceNum```
  - Normalization, randoms and scatter files are read as in MATLAB (e.g. the ```normalization``` variable of a mat-file, randoms from ```raw_SinDelayed```, ```SinDelayed``` or ```delayed_coincidences```)
  - Multi-bed (step-and-shoot) CT: ```options.bedOffset``` is added to each bed position (the number of projections must be divisible by the number of bed positions)
- ```attenuation_correction_factors(options, file)``` (MATLAB) accepts the Inveon attenuation file as a second input, so no file dialog is needed
- CT: normalization correction (```options.normalization_correction``` with ```options.normalization```) is now supported (MATLAB and Python, implementation 2 and the projection operators; previously it was always switched off for CT)
  - The normalization is input as an array, either one value for each detector pixel (```nRowsD x nColsD```, used for every projection) or one value for each measurement (```nRowsD x nColsD x nProjections```); it is not loaded from a file for CT
  - With intensity data (```options.usingLinearizedData = false```) the normalization is a detector gain on the blank: the expected measurement is ```flat * n * exp(-A x)```. Count-based algorithms (OSEM, ROSEM, OSL-OSEM, ROSEM-MAP, RAMLA, BSREM, RBI, OSL-RBI, DRAMA, MBSREM, MRAMLA, SPS, PKMA, SAGA) with corrections during reconstruction apply the gain in the reconstruction. With precorrection, FDK, or algorithms for which OMEGA linearizes the data (e.g. LSQR, CGLS, PDHG), the measurements are divided by the normalization before the linearization, so all values must then be positive (exclude dead pixels with ```options.maskFP```)
  - With linearized data (```options.usingLinearizedData = true```) and the projection operators the model is ```n * (A x)```
  - When ```options.flat``` is not given, it defaults to the maximum of the normalized data, ```max(y ./ n)```
  - Not supported (an error is raised): largeDim, extended FOV, implementations other than 2 (and ```options.useCPU```), dynamic data, ```options.normalize_scatter```, count-based algorithms with linearized data, and backprojection type 5 (projector types 5, 15, 25, 35, 45 and 55) with linearized data
  - Validation: precorrected reconstructions (LSQR, OSEM, FDK) equal the reconstructions of the same data without the gain to 2-5e-7; with the gain applied in the reconstruction, OSEM and PKMA keep the true image as a fixed point and equal the gain-free reconstruction for a constant normalization (1.2e-6)
- CT: ```options.useLayeredBP``` (default true) stores the projections read by the backprojection of projector type 4 in a layered 2D texture instead of a 3D texture (CUDA, HIP and OpenCL; MATLAB and Python). The result is identical and the backprojection is about 2-5 % faster
  - It is switched off automatically when a subset (or, without subsets, the whole data) has more than 2048 projections, or more than the layer limit of the device, and with helical CT, largeDim and Metal
- Elekta XVI CBCT data loader ```loadElektaXVIData``` (MATLAB and Python ```omegatomo.fileio```): reads the ```.his``` projection images and ```_Frames.xml``` (gantry angles and detector offsets of each projection, exposed frames only), supports binning and sets the nominal XVI geometry and a global flat (air) value (```options.elektaFlatPercentile```, default 99.9th percentile of the data). ```options.elektaUCentreScale``` and ```options.elektaVCentreScale``` (default 1) scale the detector offsets of the XML file

Changed behaviour:
- SPECT ray tracer projectors now normalize emission probability based on a pre-determined superellipse
  - Defaults to the superellipse corresponding to FOV boundary
  - Traced rays are clipped to the inside of the superellipse
  - Currently only ```options.ellipsePower``` values of 2 (elliptic cylinder, an ellipse in the transaxial
    plane) and ```Inf``` (box, the default) are supported; other superellipse powers raise an error
- Memory optimization for SPECT
  - These arrays can now be input as size nRowsD * nColsD * nHeads (one image for each detector head):
    - FP mask
    - BP mask
    - Normalization image / detector uniformity image
  - These arrays have to be input as size 2 * n_rays_transaxial * n_rays_axial x nRowsD x nColsD x nHeads (one image for each detector head):
    - rayShiftsSource
    - rayShiftsDetector
- SPECT multi-ray projector now uses ```options.n_rays_transaxial``` and ```options.n_rays_axial``` as used in sinogram reconstructions with other modalities
  - The old ```options.nRays``` is still supported, but not recommended to use.
  - ```options.nRays=n^2``` converts silently into values of ```options.n_rays_transaxial=n``` and ```options.n_rays_axial=n```
- SPECT rotation projector with Python custom operators now has better emission probability normalization when detectors are inside FOV
- PET transaxial detector coordinates (both MATLAB and Python) now form an exact regular polygon of ```blocks_per_ring * transaxial_multip``` blocks for any block count, with the block faces at ```options.diameter / 2```
  - Previously a quarter of the ring was mirrored with separate x/y gap corrections, which was a regular polygon only when the number of blocks was divisible by 4
  - Results change for scanners where ```blocks_per_ring * transaxial_multip``` is not divisible by 4 (e.g. the Biograph Vision example, 38 blocks, crystals move up to ~1.3 mm; 42 blocks up to ~0.4 mm); unchanged (< 0.001 mm) otherwise
  - ```options.diameter``` is the distance between opposite block faces
- TV-type priors (TV, weighted TV, SATV, anatomical TV types 1-2, APLS) now use the same clamp-to-edge (Neumann) boundary in all cases
  - Previously buffer mode (```options.useImages = false```), the CPU implementation and the MATLAB-based implementations (1, 4 and 5) zero-padded or mirrored the image at the volume boundary, giving different results than the default GPU image mode
  - Results change only at the volume boundary and only for those cases; the default GPU (image) mode is unchanged
- NLM, RDP, GGMRF and hyperbolic priors now clamp the image (and reference image) at the volume boundary also when ```options.useImages = false```, matching the default image mode (previously zero-padded)
- With ```options.scatter_correction = true```, ```create_atten_matrix_CT``` and ```attenuationCT_to_511``` now use the narrow-beam water attenuation 0.096 cm^-1 instead of 0.093 cm^-1
  - CT-based attenuation coefficients are 3.2 % larger and attenuation-corrected images increase by a few percent, more at the center of large objects
  - Without scatter correction nothing changes (0.093 cm^-1); ```attenuationCT_to_511``` has a new optional ```waterMu``` input
- The default workgroup/block size of list-mode and index-based data is now 32 instead of 128 (on wave64 AMD GPUs, i.e. GCN/CDNA, projector type 1 uses 512, limited by the device maximum, and the other projector types keep their earlier defaults: 128 for types 2 and 3 and 64 for type 4), as 128 was a poor size for unsorted events. A user-supplied `options.local_size` is used as before, and the sinogram, projection and other defaults are unchanged
- PET sinogram data with projector type 1 (both forward and backward projection) and the 1D launch (subset types 1-5) now defaults to a workgroup/block size of 256 instead of 128, limited by the device maximum. This is up to about 30 % faster with attenuation correction and never more than about 3 % slower (RTX 4090, RTX 3090, TITAN RTX, RX 7900 XT). NVIDIA OpenCL no longer uses 32 there. Projector types 2 and 3 keep 128
- OpenCL CPUs now default to useImages = false
- PSF padding now uses af::pad and requires AF 3.8 or newer, so OMEGA now requires AF 3.8 minimum
- PET: Projector type 2 is now normalised to the line integral: the forward projection of a constant image is the same as with Siddon (type 1) (median ratio 0.9996 in the NEMA example), independent of the tube width and the voxel size
  - The weights 1 - d/R (R = tube width, ```options.tube_width_z```/```options.tube_width_xy```) are divided by the integral of the cone, pi R^2 / 3 per unit length, expressed in voxels (divided by the voxel volume); previously the weights were unnormalised (sums of "voxel counts")
  - Type 2 reconstructed images are now larger by the factor 3 dx dy dz / (pi R^2) per unit LOR length, i.e. about three orders of magnitude (about 900 times in the NEMA example, 128 x 128 x 63 voxels, R = 2.4 mm) compared with earlier versions; the type 2 image scale no longer changes with the tube width or the number of voxels
  - Regularization parameters (beta) tuned for projector type 2 with earlier versions need to be retuned (the images are now on the same scale as with types 1 and 3)
  - Fixes the +2-3 % excess of the total number of counts of the OSEM-type reconstructions with type 2 (caused by the absolute positivity floor ```options.epps``` on the small image values); now the same as with types 1 and 3 (+0.5 %)
- Projector type 3 (volume of intersection): the scale no longer depends on the voxel size, tube radius or voxel radius, the forward projection of a constant image is now the same as with Siddon (type 1)
  - The tube volume is now scaled by Vmax / (dx dy dz), where Vmax is the volume of the voxel sphere (previously the forward projection of a constant image was too large by the factor Vmax / (dx dy dz), e.g. 1.55 in the NEMA example with 128 x 128 voxels and 2.9 with 64 x 64 voxels)
  - Type 3 images are larger by Vmax / (dx dy dz) (for example 1.55 times with 128 x 128 voxels and 2.9 times with 64 x 64 voxels in the NEMA example; compared with v2.2.1, which also had the pi^2 error, about 15 times larger with 128 x 128 voxels), regularization parameters tuned for type 3 may need retuning
- TOF weights (projector types 1-4, GPU backends) are now the exact probability that an emission in the voxel crossing is recorded in each TOF bin (the Gaussian integrated over the bin and averaged over the crossing), replacing the sampling of the Gaussian at the bin centre
  - The weights are no longer renormalised to sum to one: the sum over the bins is the fraction of the TOF kernel inside the TOF window. This fixes the edge-bin artefacts when the TOF window does not cover the object, the TOF direction of projector type 4 for rays entering past the LOR midpoint, and the TOF bin sums of projector type 3 (now equal to the non-TOF value; list-mode and sinogram TOF agree)
  - Validated against a float64 reference (per-weight error ```<= 2e-5```) and against v2.2.1, whose TOF forward projection of type 1 had 15-22 % errors at the 99th percentile in three of the four tested cases (now ```< 1e-6```)
  - TOF images change slightly (typically < 1 % in the forward projections when the TOF window covers the object); with a TOF window narrower than the object the change is larger (reconstructed activity +11 % in the GATE TOF example)
  - Speed with sinograms is about the same as before; list-mode TOF is considerably faster (one accumulator per event and clipping of the ray to the event's TOF window): without attenuation at 50-100 ps the forward projection is 3.3-6.4 times and the backprojection 1.8-3.1 times faster (about 2.1 and 1.4 times with attenuation)
  - Even numbers of TOF bins now work when loading GATE data (previously a crash or every event shifted by half a bin), a negative ```options.TOF_offset``` is now applied and user-supplied ```TOFCenter``` values are checked to be uniformly spaced
  - The CPU backend (implementation 4) still uses the old TOF model
- Index-based and list-mode reconstructions: image values are about 2 times those of v2.2.1. The sensitivity image previously counted every LOR twice (once in each detector order) and now counts it once (GPU and CPU implementations)
- Axial coordinates and segment table
  - The default ```options.segment_table``` now uses 2 x rings - 1 instead of ```options.Nz``` (MATLAB and Python); unchanged when ```options.Nz``` = 2 x rings - 1
  - The lower limit of the default segment table also uses 2 x rings - 1 instead of ```options.Nz``` (affects only ```options.ring_difference``` < rings - 1)
  - When ```options.rings``` / ```options.det_per_ring``` are not set, they are now derived from the largest axial / transaxial crystal count (and ```options.transaxial_multip```) in both languages; Python no longer fails with ```options.rings = 0```
- The CPU backend (implementation 4, and the implementation 1 MEX-files) now supports only the improved Siddon projector (projector types 1/11) and the SPECT rotation projector (6); projector types 2 and 3 and their hybrids raise an error with implementation 4
- Dynamic reconstructions with saved intermediate images (```options.save_iter```/```options.saveNIter```): the output is now ordered Nx x Ny x Nz x Nt x saved images in MATLAB and Python (previously the last two dimensions were swapped, and the images of the time steps after the first were saved in the wrong place)
- ```options.TOF_noise_FWHM``` is now the FWHM of the noise added to the time differences when loading GATE ROOT data (previously ```options.TOF_FWHM``` was used); the default -1 uses ```options.TOF_FWHM``` and 0 adds no noise (MATLAB and Python)
- Python extended FOV and extrapolation (```CTEFOVCorrection```) now follow MATLAB: the axial extended FOV length is derived from the cone-beam geometry unless ```options.eFOVLengthAxial``` is set (never negative), the default transaxial extrapolation length is 0.25 (was 0.3), ```options.extrapLengthAxial```/```options.extrapLengthTransaxial``` and ```options.eFOVLengthAxial```/```options.eFOVLengthTransaxial``` can be set separately, and ```options.useInpaint``` is MATLAB only
- largeDim now gives an error in unsupported cases: APLS, anatomical reference images with RDP, TV and NLM, and 3D backprojection/prior masks (use a 2D transaxial mask)
- More input checks with list-mode and index-based data: the index-based inputs are checked for presence and format, subset types other than 0, 1 and 3 are rejected, and dynamic TOF index-based data needs dynamic TOF indices
- ```Inveon_PET_main_listmode_example``` (MATLAB and Python) uses the geometry of the Siemens header (```options.cr_p = 1.63```, ```options.cr_pz = 1.592```, ```options.diameter = 161.08```, ```options.firstBlockAngle = 10.79```, ```options.DOI = 4.584```); results differ slightly from earlier versions
- ```GATE_OMEGA_reconstruction.py``` now requires GATE (opengate) 10.1.1 or newer (online coincidence sorter)
- ```options.ProxRDP``` and ```options.ProxNLM``` have been removed (MATLAB and Python). They never worked: the proximal kernels did not exist, so the reconstruction ran without regularization. Setting either option now gives an error; use ```options.RDP``` or ```options.NLM``` instead
- Projector type 4 (all GPU backends): the samples along the ray are now at the midpoints of the sampling intervals (previously one tenth of a step after the ray entered the image), which removes a bias of about +0.4 times the interpolation length in the forward projection (CT and PET). The CT forward projection and the voxel-based CT backprojection are now adjoint to about 1e-5 (previously 4.5e-3)
- Projector type 4 TOF: the TOF weight of each sample is now centred on the sample. Previously it was shifted by half a step towards the second detector of the LOR (e.g. 0.75 mm with a 1.5 mm interpolation length, about 1.7 mm with Joseph interpolation)
- CT with tilted detector panels (```options.pitchRoll``` or a six-row ```options.uV```) and non-square pixels: the detector direction vectors are now scaled as vectors (the row direction by ```options.dPitchX```, the column direction by ```options.dPitchY```) instead of per x/y/z component, which made the pixel grid skewed with wrong pixel sizes (MATLAB and Python, including ```CTDetectorCoordinatesFull```)
  - Results change only for tilted panels with ```dPitchX ~= dPitchY```; e.g. with a 20 degree roll and 1.5 x 1.0 mm pixels the forward projection error of projector type 1 against the exact line integrals dropped from 6 % to 1.5e-6
- COSEM, ECOSEM, ACOSEM and OSL-COSEM now give an error with CT data (MATLAB and Python), as they have no transmission tomography formulation
- The fast PDHG (```options.fastPDHG```) is not used with a 3D backprojection mask and backprojection type 4; the regular path is used instead
- Subset type 12 now gives an error without offset correction (```options.offsetCorrection = true```), MATLAB and Python
- Parallel beam CT (```options.useParallelBeam```) now gives an error with projector type 5 (also as the forward or backward projector of a combined type, e.g. 15, 45, 51 or 54) and with the CPU implementation (implementation 4 and ```options.useCPU```/```options.use_CPU```), MATLAB and Python. Neither had parallel beam support and both silently used the cone beam geometry
- Parallel beam FDK/FBP now uses the same zero-padded ramp filter and DC correction as cone beam FDK (MATLAB and Python)

Other changes:
- Unified SPECT rotation projector logic for all three Python operator backends
- The voxel bounds check in the Siddon projector (projector types 1-3) is now cheaper, which makes the forward projection about 7 % faster on full ordered sinograms. The results are unchanged
- Combined OpenCL, CUDA, Metal, and HIP code into one
  - Functionality between the three should now be identical
  - Doesn't require anything from the user
- Optimized the code in general
  - Projector type 5 and the hybrid x5 should see 10-20% improvement in speeds
  - Otherwise the improvement should be 5-10%
  - Biggest effect for (CB)CT, less so for emission tomography
- Reduced memory usage of implementation 2
  - PDHG (and its variants L1 and KL), PKMA, MBSREM, and BSREM now use less memory as the backprojection is never saved to global memory
  - The above is automatic and doesn't work with all possible combinations, such as when using anatomical references images with priors
  - All other methods should also see a reduction of one volume in memory usage
  - ```options.fastPDHG = false``` disables this for the PDHG variants (MATLAB and Python)
- Commented unused code in CPU RDP (implementation 2 CPU)
- ArrayFire is no longer mandatory with the custom reconstruction in Python
- Projector types 2 and 3: TOF forward projection is faster with CUDA (12-18 %, fewer registers) and HIP (about 9 %); OpenCL is unchanged
- GATE ROOT loading is much faster (only the needed branches are read, multithreaded; MATLAB, Octave and Python), and Inveon list-mode loading is faster and uses less memory in sinogram mode
- Implementation 2: the backprojections of multi-resolution reconstructions run concurrently, and unnecessary synchronizations were removed
- Projector type 3: the volume table is now in global memory instead of constant memory
- Python: the power method computes the Lipschitz values of multi-resolution reconstructions separately for the main and the extra volumes
- ```compile.py``` locates Visual Studio with vswhere (Visual Studio 2017-2026, including Build Tools), finds the OpenCL library also in ```OCL_ROOT/lib``` and no longer fails when the OpenCL path is not set; ```install_mex``` uses the included OpenCL headers for implementations 3 and 5 on Linux
- Python ```volume3Dviewer```: arrow buttons step the slice and time sliders one step at a time, as in the MATLAB viewer
- ```SPECT_custom_algorithm_example_PyTorchCUDA.py``` was renamed to ```SPECT_custom_algorithm_example_PyTorch.py``` (also runs with Metal)
- ```install_mex``` also looks for the OpenCL library in ```OCL_ROOT/lib``` when ```OCL_ROOT/lib/x86_64``` does not exist (as ```compile.py```)
- CT backprojection type 4 with a normalization now interpolates the normalized measurements, making it the adjoint of the normalized forward projection (adjointness 2.6e-4 to the 7e-5 level of the unnormalized pair in the tested case). The MATLAB and Python front ends still disable the normalization correction for CT
- The projector type 4 documentation now states that the ray-based backprojection of non-CT data deposits each sample to the nearest voxel and is therefore not the exact adjoint of the trilinear forward projection (the row sums match)

Bug fixes:
- OSL-OSEM denominator (sensitivity image + prior gradient) is now limited to at least 10 % of the sensitivity image (and at least ```epps```), preventing sign changes, division by zero and blow-up with strong priors (implementation 2 and Python implementation 4)
- Multi-resolution reconstruction now gives an error with the Python CPU implementation instead of crashing
- OSL-OSEM (implementations 2 and 4 in Python, implementation 2 in MATLAB) without a prior now equals OSEM: the backprojection term was missing ```epps```, which set voxels outside the measured region to zero, and the extra volumes of multi-resolution reconstructions were never updated
- Python: the NLM filter parameter set in the examples (```options.sigma```) was ignored; the examples now use ```options.NLMsigma```, and ```options.sigma``` is used as a fallback when ```options.NLMsigma``` is not set
- Python: the anatomical TV weight now uses ```options.T``` as documented (previously the undocumented ```options.B```, default 0.01, was used and ```options.T``` was ignored)
- Fixed texture input for standalone Metal projector calls
- Fixed the missing Metal definition for projector type 4
- APLS prior gradient was incorrect with implementations 2 and 4 (GPU and CPU); it now matches the original article (Ehrhardt et al. 2016) and no longer produces NaN values
  - The reference image is now clamped at the volume boundary and ```epps``` is no longer added to the reference image gradient (also in MATLAB implementation 1/4/5 path)
- APLS reference image can now be input as a variable in MATLAB (previously only mat-files worked) and as a NumPy array in Python on Windows
- APLS reference image size is now checked (wrong size could previously cause out-of-bounds reads)
- MATLAB-based TV type 1 (anatomical weights) gradient was incorrect on the last slices and everywhere in 2D (implementations 1, 4 and 5)
- Anatomical TV type 1 gradient was incorrect with implementations 2 and 4 (GPU and CPU): wrong sign of the voxel's own term and the neighbor terms used the weights of the wrong voxel
- Anatomical TV type 2 gradient was incorrect with implementations 2 and 4: the neighbor terms were missing the reference image term
- Anatomical TV type 2 now uses ```options.T``` as the reference image weight with implementations 2 and 4, as with the MATLAB-based implementations (previously ```options.C```)
- Hyperbolic prior gradient was incorrect with implementation 2 when the image size in x or y was not a multiple of the local work size (16)
- Projector type 3 (volume of intersection) normalized the tube volume incorrectly with the GPU backends (OpenCL, CUDA, HIP, Metal) in versions 2.0–2.2.1: the volume was computed as L·R²/π instead of π·L·R² (the CPU backend was correct)
  - Forward projections were π² (about 9.9) times too large and reconstructed and sensitivity images π² times too small compared with the other projectors and the CPU backend
  - Type 3 reconstructions now have about 9.9 times larger absolute values than with version 2.2.1; regularization parameters (beta) tuned for type 3 with version 2.2.1 may need retuning, and calibration factors derived with version 2.2.1 type 3 must be rescaled by π²
- The Khronos OpenCL C headers (```CL/cl.h```, ```cl_ext.h``` etc.) matching the included C++ bindings (```opencl.hpp```) are now included in ```source/cpp/CL```, and they are used instead of the headers of the installed OpenCL SDK
  - Fixes OpenCL build failures when the system OpenCL headers are older or newer than the included ```opencl.hpp```, e.g. with the OpenCL headers of CUDA 12.x or Ubuntu 22.04 (errors about ```CL_KHR_COMMAND_BUFFER_EXTENSION_VERSION``` or ```cl_command_buffer_update_type_khr```)
  - Only the OpenCL library (```OpenCL.lib```/```libOpenCL.so```) is now needed from the OpenCL SDK
- ```create_atten_matrix_CT``` failed or cropped the CT incorrectly depending on the CT width and the transaxial FOV (the cropped image was one pixel short for some FOVs, the window was off-center and the scale was rounded). It now uses exact area-weighted resampling onto the image grid, supports a FOV larger than the CT (padded with air) and non-square ```options.FOVa_x```/```options.FOVa_y``` with ```options.Nx```/```options.Ny```
  - The CT values are now converted to HU with the DICOM ```RescaleSlope```/```RescaleIntercept```. Previously a fixed offset of 1000 was subtracted, so for e.g. Siemens CT (intercept -1024) all values were 24 HU too high, i.e. the soft tissue/water attenuation coefficients were about 1 % too high at 100-120 kVp
- ```create_atten_matrix_CT``` used ```exist(..., 'folder')```, which older MATLAB releases reject; it now uses ```'dir'```
- Python now raises a ```RuntimeError``` when the native reconstruction fails (the native library previously always returned success, so a failed reconstruction only printed "Reconstruction failed!" and returned a zero image)
- The default MATLAB initial value ```options.x0``` is now 1 (1e-4 for CT), as documented and as in Python, instead of 1e-5; the old value could put the OSEM/MLEM forward projections below ```epps``` and scale the image incorrectly
- The backprojection of the projector types 1-3 again skips the measurements whose input is exactly zero (when the sensitivity image is not being computed), which speeds up reconstructions of low-count sinograms. This check was lost in v2.0
- PDHG, PDHG-KL, PDHG-L1, PDDY and CV with subsets converged to a biased image in versions 2.0–2.2.1: the backprojected dual variable was initialized to the initial image ```options.x0``` instead of zero. This also affected the CT fast PDHG and all largeDim PDHG reconstructions (implementation 2 and the MATLAB-based implementations)
  - The reconstruction effectively had an extra hidden penalty (the sum of ```options.x0``` times the image), which pulled the image towards zero, most in low-sensitivity regions. The bias did not decrease with more iterations and the result depended on ```options.x0```
  - In a PET test (```options.x0 = 1```, RDP prior) the relative distance to the true solution after 100 iterations dropped from 36 % to 12 %
- PDHG, PDHG-KL, PDHG-L1 and CV without subsets evaluated the prior gradient at the extrapolated image, which can be negative, and could diverge (e.g. NaN values with RDP). The prior gradient is now evaluated at the primal image. ProxTV, ProxTGV and TGV, whose dual updates use the extrapolated image, are unchanged
- With ```options.useImages = true``` (default) and other forward projectors than type 5, the NLM, RDP, GGMRF, hyperbolic and TV priors with the GPU backends (OpenCL, CUDA, HIP) used the last forward projected image instead of the current estimate
  - Dynamic reconstructions: the priors of every time step used the image of the last time step
  - With PSF (```options.use_psf```): the priors used the PSF-blurred image
  - BSREM, ROSEM-MAP, SART/POCS and PDDY: the priors used the image from before the update
  - With largeDim, these priors crashed
  - The results now match ```options.useImages = false```
- RDP with a reference image (```options.RDP_use_anatomical```) fixes
  - The precomputation of the reference weights overwrote the sensitivity image (OpenCL, wrong results even with zero beta) or crashed (CUDA/HIP)
  - CUDA/HIP used only one slice of the reference image
  - CT reconstructions with a reference image always failed
  - The fast PDHG ignored the reference image (the prior was zero)
  - largeDim with a reference image now gives an error (it was never supported)
- Fixed a GPU memory leak with largeDim and the CUDA backend
- Proximal TV gradient is now correctly computed when using CUDA (previously it used Z/X/Y gradient order instead of the correct X/Y/Z)
- GGMRF now works with masking
- GGMRF works correctly when using largeDim
- RDP works correctly in OpenCL when using largeDim
- Projector type 4 is now closer to a "normal" case when using largeDim
- Dynamic listmode (coordinate-based) should now work
- Offset correction can now be used with subset type 0
- Masks should work more reliably across different setups
- Fix for offset correction when using non-square detector pixels
- Correct NLM for CPU, previously the results weren't correctly saved (implementation 2 CPU)
- GGMRF previously could have used wrong weight values in very rare cases (implementation 2 CPU)
- NL type 6 (NLGGMRF) now works on the CPU side as well (implementation 2 CPU)
- Small loop fix to the CPU projector type 1 when computing perpendicular rays, one index was previously omitted (implementations 1, 2 CPU, and 4)
- Arc correction in MATLAB should now work for dynamic cases
- MBSREM, MRAMLA and SPS now work in MATLAB when using randoms correction (implementations 1 and 4)
- Gradient computation for the gradient-based preconditioner was previously wrong in MATLAB (implementations 1 and 4)
- DRAMA should now work in MATLAB (implementations 1 and 4)
- Image-based filtering preconditioner should now work in MATLAB (implementations 1 and 4)
- FISTA with L1 regularization now works in MATLAB  (implementations 1 and 4)
- BSREM and ROSEM MAP should correctly regularize now in MATLAB  (implementations 1 and 4)
- Median filter should properly work when using rectangular (non-square) neighborhoods  (implementations 1 and 4)
- Cell measurement inputs should work more reliably (MATLAB)
- Power method should work now if the input image is a cell, and use the input number of iterations and filtering if selected (implementations 1 and 4)
- Weight error now works correctly if the weights vector is too long (MATLAB)
- TV prior should now work with reference image (implementations 1 and 4)
- Fix small indexing bug in arc correction in Python
- TOF indices are now correctly imported when loading ROOT files (GATE 9 and earlier)
- Fixed compilation when platform included several different devices
- Correctly release memory with CUDA when using index-based reconstruction
- Fix RDP with corners when using largeDim, previously one slice was always missing
- Fix TOF indices when computing sensitivity image, previously there was a risk for out of bounds case
- Fixed PSF padding. Now the padding is properly reflected in the boundary
- Fix parallel beam reconstruction when using CUDA
- Inveon SCN-files can be loaded again
  - ```options.fpath``` pointing to an existing .scn file is used without a file dialog
  - With dynamic SCN files ```options.partitions``` must match the number of frames in the file (an error is raised otherwise)
- Multi-ray TOF should now work, though in general it's not recommended
- TGV works when using masks and not using images i.e. useImages = false
- Implementation 2 with CPU no longer gives errors in MATLAB with projector type 11
- Any algorithm that didn't use sensitivity image now works when useImages = false
- Fixed forward projection masks when useImages = false
- Fixed MRP (implementation 2 CPU)
- PET attenuation now works when useImages = false
- With TOF, non-TOF scatter data (```options.ScatterC```) was added to every TOF bin at full value in the reconstruction (corrections during reconstruction) and was not accepted by the precorrection (MATLAB); it is now divided evenly among the TOF bins (as the randoms are) and a message is printed (MATLAB and Python). TOF-resolved scatter data (```Ndist x Nang x NSinos x TOF_bins```) can be used with the precorrection, and now also with corrections during reconstruction (MATLAB and Python): the randoms are divided evenly among the TOF bins and added to the TOF-resolved scatter, which is used as it is in each bin (```options.randomsTOF```, set automatically). It is supported with implementation 2 (the Python library), static data, no increased sampling, and no scatter smoothing or scatter variance reduction; other combinations raise an error, as does any scatter data size that is neither a non-TOF nor a TOF sinogram. The additive term is then TOF-sized, i.e. it needs as much (host and device) memory as the TOF measurements (twice the TOF data memory in total). Non-TOF reconstructions, non-subtracted scatter (```options.subtract_scatter = false```), raw data, list-mode and index-based data are unchanged
- Fixed randoms and scatter loading with dynamic or TOF data. The randoms are divided by the TOF bins also when the scatter is not subtracted (MATLAB and Python; the non-subtracted scatter itself is unchanged), scatter-only data in a cell array keeps all time frames, and a single scatter sinogram can be shared by all the time frames (MATLAB and Python) and is normalized and smoothed only once, and a single randoms sinogram can be combined with scatter data that has its own time frames (MATLAB). TOF-resolved scatter data is still supported only with one time frame
- Added the missing ```useSingles``` default (True, as in MATLAB) to the Python options. The SPECT DEW/TEW scatter correction no longer raises an AttributeError
- ```options.use_user_normalization``` now works as in MATLAB in Python (PET). If it is False, the default normalization file is used and an error is raised if it is not found, no file dialog is shown. If it is True, ```options.normalization_datafile``` is used, or a file dialog is shown if it is empty. A ```normalization_datafile``` set while ```use_user_normalization``` is False is an error instead of being silently ignored
- The normalization and attenuation file dialogs now give a clear error naming the option to set if no display is available (e.g. over ssh) instead of failing or waiting forever
- Projector types 2 and 3: the tube of the ray now ends at the detectors; previously voxels beyond the end points were included when the image was larger than the detector ring (e.g. forward projection of a constant image 5 % too large with a 226 mm and 36 % with a 302 mm field of view on a 260 mm ring, and sensitivity/backprojection images with the same excess). Images within the ring are unchanged
  - The tube was an infinite cylinder around the line through the two detectors, limited only by the image box, so it extended past the end points whenever the image box (side larger than the ring diameter / sqrt(2), 184 mm for a 260 mm ring) contained them; the foot point of each voxel on the ray is now required to lie between the end points (not applied to SPECT, which has its own rule; projector types 1, 4 and 5 are unaffected)
  - Forward projection of a constant image divided by the Siddon one is now 0.998 (type 2) and 0.999 (type 3) for every field of view tested (151-302 mm)
  - The end-point check is enabled automatically (kernel define ORTHCLIP, MATLAB/Python/C++ alike) only when a detector lies inside the image box widened by the tube reach, so images inside the detector ring are unchanged (bitwise identical to before on CUDA and OpenCL) and not slower. On a 260 mm ring it is off for fields of view up to about 151 mm and on from about 184 mm (CT: always on)
  - With TOF the excess was partly hidden by the TOF window, but it is removed in the same way (e.g. a window covering the whole ring and a 226 mm field of view: TOF bin sum 4.6 % above type 1 before, 0.1-0.2 % below after)
  - When the check is on (a detector inside the image), forward projection of types 2 and 3 on CUDA is about 8-22 % slower and backprojection 0-3 % slower (one case, type 3 backprojection with the narrow TOF window, 45 %), because of the extra test per voxel
- CT: projector type 2 is now rejected also as the forward or backward projector of a combined projector type (e.g. 12, 21, 23, 24, 42) in both MATLAB and Python (previously only plain type 2 was rejected)
- Normalization correction (sinogram and index-based data)
  - With ```options.normalize_scatter```, a zero normalization no longer divides the scatter, the scatter keeps its value where the normalization is zero (sinogram data in MATLAB and Python, index-based data in Python). Index-based data needs the scatter as one value for each event, otherwise an error is raised
  - Zero coefficients in a normalization that is loaded from a file (```.nrm```, mat-files, including the default mat-file) are handled the same way in MATLAB and Python: a zero coefficient is a bin without efficiency, so the inverted normalization is zero there (no Inf; the default mat-file previously gave one in MATLAB and Inf in Python), the scatter keeps its value and the precorrection gives zero. SPECT normalization files are unchanged
  - ```options.normalize_scatter``` now works with a normalization that is loaded from a file in MATLAB (in the reconstruction): the normalization is loaded before the scatter is normalized, as in Python. With ```NSinos``` < ```TotSinos```, the scatter can have either ```NSinos``` or ```TotSinos``` sinograms in both languages (the first ```NSinos``` are used)
- Multi-layer PET
  - Fixed: delayed coincidences of mixed-layer events in the MATLAB ASCII import were placed in the wrong layer block
  - Pseudo rings now work with multiple layers (positions from the axial block size, previously the transaxial one); an error is raised for pseudo rings with layers that have a different number of axial crystals
  - Fixed: with multiple layers, ```options.offangle``` is now applied to each layer separately (with pseudo detectors or layers of unequal size the layer split was wrong and the other layer was not rotated correctly)
  - Fixed: the second layer was placed 2 x ```options.DOI``` too far out (the DOI was added twice); it is now at the inner diameter plus 2 x ```options.crystH(1)``` plus the DOI, as the first layer (MATLAB and Python)
  - Fixed: ```options.flip_image``` reversed the whole detector array, which swapped the layers; each layer is now flipped separately
  - MATLAB: with multiple layers the first two outputs of ```detector_coordinates``` (Inf padding removed) are now also flipped/rotated (previously only the padded outputs were rotated with ```options.offangle```)
- Fixed: ```options.flip_image``` with pseudo detectors (or with layers of unequal size) moved the pseudo detectors/padding from the end of each block to its start, so that real detectors got the coordinates of the gaps; the ring is now mirrored with the block structure kept (MATLAB and Python). Without pseudo detectors the coordinates are unchanged
- MATLAB GATE ASCII data import now follows the ROOT import, and dynamic ASCII data works
  - The time window (```options.start```/```options.end```) is applied as in the ROOT import (events with ```start <= t <= end```); previously the upper limit was applied incorrectly
  - The time frame of each event is now determined from the time of the event, independent of the order of the events (previously from the first/last event of each frame), and is passed to the sinogram creation as a zero-based index
  - The sinogram creation uses ```createSinogramASCIICPP``` when it has been built (MATLAB R2018a and later), otherwise ```createSinogramASCII```; the input types are now explicit
- Axial coordinates and segment table
  - Fixed a crash in the span > 1 axial coordinates when ```options.ring_difference``` < rings - 1: the planes beyond the ring difference get no coordinates (Inf) and the boundary segment uses only the accepted ring pairs
  - Fixed the raw-data axial coordinates with pseudo rings (MATLAB counted the pseudo-ring positions instead of their number, Python deleted the wrong rings)
- Precorrection (MATLAB sinogram formation)
  - Fixed TOF randoms precorrection subtracting the full delayed sinogram from every TOF bin; the randoms are now divided evenly between the TOF bins
  - Scatter variance reduction now uses ```options.scatter_variance_reduction``` (previously ```options.variance_reduction```)
  - Plain randoms precorrection is now recorded in ```appliedCorrections```; fixed an undefined variable with user-supplied randoms
  - ```appliedCorrections``` no longer loses the scatter variance reduction entry when scatter smoothing is applied
- Precorrection and corrections (Python)
  - TOF randoms precorrection now divides the randoms evenly between the TOF bins (as in MATLAB)
  - The PET normalization precorrection is now applied after the randoms and scatter subtraction (as in MATLAB), and is applied per LOR to every TOF bin
  - Trues and scatter reconstructions no longer load or apply randoms/scatter, but do apply the normalization correction if selected
  - Measurement files with an ```appliedCorrections``` struct (data stored by OMEGA) are now handled as in MATLAB: the stored precorrected sinogram is used only if it has been precorrected with the selected corrections and options (normalization, randoms and scatter, variance reduction, smoothing); otherwise the raw sinogram is used and precorrected again, and the corrections already in the stored sinogram are not applied twice
- Python sinogram creation (```omegatomo.util.saveSinogram```)
  - Fixed the segment table input: the cumulative segment table (```np.cumsum(options.segment_table)```) is now used, also in the default table; a non-cumulative table raises an error
  - The size variables are converted to integers and the time frames are passed as integers (zero-based), which fixes dynamic sinogram creation; the input sizes, time frames, TOF bins and layers are now checked
- FDK and Parker weights (CT, MATLAB and Python, all GPU backends)
  - The FDK weights now include the magnification of each projection, in the same way for projector types 4 and 5; absolute values of FDK images change, especially with large magnification
  - The FDK distance weight is computed from the actual source position of each projection instead of the projection angle, so custom geometries, clockwise scans, arbitrary start angles and other orbit radii work
  - The FDK ramp filter is zero-padded (FFT length at least 2 x ```options.nRowsD``` - 1) and keeps the DC term; previously detectors with 2^n rows had no padding, which biased the low frequencies by several percent
  - Parker weights (```options.useParkerWeights```) use the actual source and detector coordinates of each projection, work with clockwise scans and angle jumps (e.g. 360 to 0 degrees), include the angular quadrature, and support negative and per-projection ```options.detOffsetRow``` (MATLAB previously ignored ```options.detOffsetRow``` because of a typo)
  - FDK with largeDim uses the total number of projections and the global projection angle of each chunk, the offset correction works with largeDim, and the last projection is also filtered
  - FDK is now correctly weighted with the BDD backprojection (projector type 5)
  - The FDK filtering is done in chunks of projections when the FFT buffer would exceed the ArrayFire limit of 2^31 elements
- Arc correction (sinogram data, MATLAB and Python) was rewritten: the interpolated sinogram is no longer rotated relative to the arc-corrected coordinates, the detector radius includes ```options.DOI```, ```options.flip_image``` is handled, and the normalization, randoms/scatter and ```corrVector``` applied during the reconstruction are interpolated with the same weights. ```options.arc_interpolation``` 'linear' and 'nearest' use fast sparse interpolation weights
- Helical CT: the backprojection used only about half of the projection data (curved detector, MATLAB and Python); the DICOM CT-PD loaders (```loadDICOMCTPD```, MATLAB and Python) now use the correct detector geometry (arc centre and channel offsets); helical CT with the CPU implementation and FDK with helical CT now give an error
- The additive term (scatter and/or randoms) used wrong values from the second subset on with projection data (CT/SPECT) and subsets (implementation 2)
- Projector type 5 (BDD)
  - The normalization weighting (e.g. detector uniformity/normalization images) was wrong in the forward and backprojection (all GPU backends)
  - With ```options.meanFP```/```options.meanBP``` the mean add-back used an unclipped footprint area, giving errors where the footprints extend past the volume or the detector
  - MATLAB projection operators (implementation 5): ```options.meanFP```/```options.meanBP``` were not passed to the kernels, the multi-volume forward projection read wrong offsets and the multi-volume backprojection with 3 or more volumes split the output at wrong offsets
- Projector types 1-3
  - Perpendicular rays were indexed off by one, and TOF did not work, when the source/detector is inside the FOV
  - The multi-ray coordinate division was correct only for perpendicular rays; the ray angle is now taken into account
  - The first voxel and the crossing times of a ray are computed from the same (integer) plane indices and nearly axis-parallel rays use the axis-parallel branch, which removes round-off errors in edge cases
  - Projector type 3: the volume table could contain NaN values at the boundary (MATLAB and Python), giving slightly wrong values and differences between OpenCL and CUDA
- Proximal TV and TGV (ProxTV, ProxTGV, TGV) with CUDA and HIP did not update the last voxels of the image when the number of voxels was not a multiple of 64
- CUDA and HIP kernels now use 64-bit global indices as OpenCL already did; very large volumes could previously overflow silently
- Projector type 6 (SPECT): the attenuation image is now correctly translated and rotated
- Multi-resolution/extended FOV initial values: Python could give non-deterministic or too small values for the extra volumes and MATLAB too large volumes with a user-input ```options.x0```; an undersized ```options.x0``` now raises an error in Python
- largeDim: CT data was linearized a second time with ```options.usingLinearizedData```, and the power method now normalizes over the whole volume instead of each subvolume
- Dynamic reconstructions (Nt > 1)
  - The sensitivity image now follows the time step when the attenuation image or the multiplicative ```corrVector``` differs between time steps (previously the first time step was used for all; implementation 2)
  - Measurement, randoms, attenuation and ```corrVector``` data given per time frame are stacked frame by frame (MATLAB and Python); with one subset or largeDim the frames were previously interleaved
  - A per-frame ```corrVector``` is ordered into subsets frame by frame, and the scatter data of time steps > 0 no longer reads past the end when only one frame was given
  - Python: static normalization, attenuation, ```corrVector``` and sinogram corrections are sliced per frame correctly; dynamic list-mode data with normalization, randoms or scatter correction disables the correction with a warning (as in MATLAB)
  - Dynamic SPECT: the attenuation image of each time step is used with the rotation projector, the detector head vector covers all time frames (MATLAB, previously frames after the first read past the end of the array) and a single normalization is shared by all frames (Python)
- MATLAB implementations 2 and 3: non-subtracted scatter (```options.subtract_scatter = false```) was silently ignored or crashed; it is now merged into ```options.corrVector``` and applied as a multiplicative correction
- Subsets
  - Subset type 4 gave wrong sinogram indices when the number of angles was not divisible by the number of subsets (MATLAB and Python)
  - Python: subset type 3 dropped the last measurement of the last subset
  - Subset type 0 now works with sinogram and projection data (MATLAB and Python)
- SPECT
  - MATLAB TEW scatter precorrection (```options.corrections_during_reconstruction = false```) used the wrong sign for the upper window term
  - With axial flipping the scatter projections (DEW/TEW windows, ```options.ScatterC```) are flipped together with the measurements
- PET geometry and data loading
  - ```options.ringGaps``` were too small with span 1 and too large with span > 1 (MATLAB and Python); index-based data now includes the ring gaps
  - Python: the segment table was one segment short and some spans failed in ```sinogramCoordinates3D```
  - Odd numbers of radial bins (```options.Ndist```) are loaded correctly with Inveon and GATE data, and ```options.ndist_side = 0``` is treated as 1
  - GATE ROOT import: dynamic data (also dynamic TOF with delayed coincidences and index-based dynamic data) loads correctly, the start/end time window is applied also to the delayed coincidences, missing delayed coincidences no longer crash, and the default MATLAB ROOT loader, the legacy loader and the Octave dynamic load work again
  - GATE ROOT sinogram import: events whose sinogram index is outside the sinogram are skipped with a message instead of writing past the end of the sinogram (including the multi-layer bound)
  - ```options.reconstruct_trues``` now works with index-based reconstruction (GATE ROOT data): only the true coincidences are stored and no randoms are used
  - Multi-layer sinograms: the layer pair of each sinogram is stored with the coordinates and ```saveSinogram``` maps it to the correct block
  - Multi-ray PET used the transaxial and axial crystal sizes the wrong way round (MATLAB and Python)
  - Inveon ```.atn```/CT UMAP attenuation (```attenuation_correction_factors```): the attenuation coefficients were half of the correct ones (extra factor 1/2), so attenuation-corrected Inveon images change
- Precorrection
  - Normalization is no longer switched on when precorrected data is used without the normalization correction
  - Randoms/scatter stored as integers are converted to single/double correctly
  - Python: a user-supplied PET normalization (```options.normalization```) is now multiplied with the measurements in the precorrection, as in MATLAB (previously divided); normalization from a file and SPECT are unchanged
- Python
  - ```options.flipImageX/Y/Z``` flipped the attenuation image and the measurement data along the wrong axes (and crashed for Z)
  - Per-projection forward projection masks (```options.maskFP```) with subsets were read at wrong offsets with CuPy/PyOpenCL
  - Saved intermediate images (```options.saveNIter```/```options.save_iter```) and stored forward projections are reshaped correctly (also with several subsets and ```subsetType >= 8```); if a reshape still fails a warning is given and the flat arrays are returned
  - C-ordered coordinate arrays are no longer flattened in Fortran order in the projector (fixes 2D data)
  - TV type 2 without an anatomical reference did not fall back to type 1
  - The 2.5D orthogonal projector used the 3D tube width, and the extended FOV prior region was off by one voxel; the prior mask (```options.usePriorMask```) is no longer applied in the backprojection
  - BSREM is treated as a MAP algorithm in the prior setup
  - SPECT: the DEW/TEW scatter is subtracted before the normalization and the detector head order of ```DetectorVector``` was fixed
  - The randoms file dialog failed (```options.fpath```), ```pitchRoll``` used the wrong detector pitch with non-square pixels, and the offset correction limit was off by one
  - On Windows the native libraries are found also with newer Python versions, which no longer use PATH for DLLs; OMEGA registers the DLL directories from ```AF_PATH```, ```CUDA_PATH``` and ```HIP_PATH```/```ROCM_PATH``` itself
  - Fixed a crash/heap corruption at exit when PyTorch is imported before an implementation 2 reconstruction with OpenCL (the ArrayFire FFT plan cache is now released after the reconstruction)
- MATLAB
  - The orthogonal and volume-of-intersection tube flags of the hybrid projector types (12, 21, 13, 31, 23, 32, 24, 42, 34, 43) are set as in Python, and type 41 is scaled as type 4
  - ```options.gradV2``` never got its default value 2.5 (affects the gradient-based image preconditioner)
  - The maximum number of sinogram angles error is no longer raised for list-mode data with custom coordinates or index-based data
  - SPECT sinogram reconstruction with one subset no longer requires subset types 8-11, and the orthogonal ray tracer with several rays now gives an error (previously the multi-ray branch was silently used)
  - The CPU MEX-files convert ```options.vaimennus``` to the expected type (previously garbage attenuation values could be read); implementation 1 forward projection with several subsets and ```forwardProjection``` with an empty input work
  - Custom list-mode detector coordinates given as two rows were ordered incorrectly, and ```index_maker``` accepts dynamic data with a cell ```options.x```
  - The automatic axial extended FOV length can no longer become negative, and an extended FOV smaller than one multi-resolution voxel is disabled with a warning instead of crashing (also in Python)
  - ```saveImage```, ```saveInterfile```, ```saveMetaImage``` and ```scaleImages``` work with the current ```image_properties``` fields (```saveImage``` used ```Ny``` for the axial voxel size and ```saveMetaImage``` opened its output file read-only); ```saveInterfile``` and ```saveMetaImage``` handle dynamic (4D) data
- MATLAB: ```install_mex``` built the OpenCL implementations 2 and 3 on Linux against the system OpenCL headers instead of the included ones (the source folder was never added to the include path), which failed when the system headers are older than the included ```opencl.hpp``` (e.g. ```CL_MAKE_VERSION``` errors)
- Projector type 4
  - PET/SPECT backprojection with attenuation and ```options.useTotLength = false``` was normalized by the full LOR length instead of the length of the ray inside the image, unlike the forward projection (factors of 0.23-0.89 per LOR in the tested case)
  - CT backprojection with a 3D backprojection mask: each block of 8 z-planes used the mask value of its first plane
  - CT backprojection: the sensitivity image (used e.g. by SART/POCS, RBI and DRAMA) only counted the rays with non-zero measurements; it now counts every ray that hits the detector
  - CT backprojection now computes its scaling from the actual pixel area and source-to-detector distance of each projection (from the detector vectors), which fixes the scale with custom or tilted detector panels and with a varying source-to-detector distance (adjointness 0.42 to 2e-4 in a tilted test). Standard flat-panel geometries are unchanged
  - The FDK arguments of the CT backprojection were passed based on ```options.sourceToCRot > 0``` instead of the FDK settings that compile them in
  - Projector types 44 and 55 (aliases of 4 and 5) silently used projector type 1 in MATLAB and with the Python library, were rejected by the Python reconstructions, failed with the Python projector operators and were not handled by the MATLAB projection operators; they now work everywhere as 4 and 5. Atomics are now also disabled for all CT backprojection types 4 and 5 (also 24, 25, 34, 35, 44 and 55)
  - Python projector operators: Joseph interpolation (negative ```options.dL```) was never selected, so negative values gave invalid projections
  - Python projector operators: the PET/SPECT backprojection with measurement-based attenuation (```options.CT_attenuation = False```) crashed with CuPy and failed with PyOpenCL (missing attenuation argument)
  - Python projector operators: the CT backprojection with a normalization was missing the normalization argument
  - Parallel beam CT backprojection was rewritten. It did not compile without FDK, and with FDK it looped over every detector pixel for every voxel and used the cone beam distance weight (an FDK reconstruction of a cylinder was 1.92 times too high). Each voxel is now projected along the ray direction of the projection onto the detector, with the weight voxel volume / pixel area perpendicular to the ray. FDK reconstructs a cylinder to within 0.2 % for 180 and 360 degree scans, the backprojection is adjoint to the forward projection at the same level as cone beam, and it is about 15 % faster than the cone beam backprojection. Offset correction, subset type 12, the CT normalization and ```options.useLayeredBP``` work as with cone beam
  - Python projector operators: parallel beam (```options.useParallelBeam```) was ignored, so the operators silently used the cone beam geometry
- Python: CT with one detector offset per projection (```options.detOffsetRow``` as a vector) failed when the coordinates were computed (```CTDetSource```); MATLAB was not affected

## OMEGA v2.2.0

### New features
- Major changes to dynamic reconstruction
  - Loop order changed from timestep -> iteration -> sub-iteration to iteration -> sub-iteration -> timestep
  - Allows for regularization through temporal dimension (as well as spatial one)
  - Implemented temporal total variation (`options.temporalTV`) and first-order derivative smoothness (`options.temporal_smoothness`) priors
  - SPECT and PET only, though might work with CT data as well
  - Dynamic reconstruction is supported with (at least) OSEM, RAMLA, LSQR, PDHG, FISTA, PKMA, SART, ASD-POCS and BB
  - Added an example script for PET with both spatial and temporal regularization, and supporting both sinogram and list-mode reconstructions
  - Implementation 2 only!
	
- Multi-resolution support for SPECT
  - Added option to normalize emission probability by either ray-FOV (ray-mask) intersection length or full ray length in SPECT reconstruction
  - This is controlled by `options.useTotLength`, where false is the first case and true the latter
  - Note that the PET functionality of `options.useTotLength` remains the same, i.e. if true, uses the length of the ray between the two detectors, otherwise the length in the image only
  
- Improved FDK reconstruction for CBCT
  - Added (optional) Parker weighting
  - Should be used with any <2pi measurement
  - Activate by setting `options.useParkerWeights` as true
  - The default "weight" value is 0.25, but using 1 gives the "original" Parker weights
  - Adjust the above value with `options.ParkerWeight`
  - Implementation 2 only!
  
- Improved volume3Dviewer
  - Works now with 4D data as well
  - With 4D data, there's a second slider for the 4th dimension
  - Python version now works in interactive sessions and also with 4D data
  
- Improved atomic operations with PET and SPECT
  - What this primarily means is that OpenCL should now be equally fast, and most of the time faster, than CUDA
  - PET and SPECT only!
  - The 32/64-bit atomics are no longer recommended and have been removed from almost all examples
  - Atomics are now computed with floats in practically every case (certain specific hardware might still benefit from the old 32/64-bit atomics, such as integrated GPUs)
  - OpenCL is now recommended for all hardware
  
### Bug fixes and enhancements

- Projector type 4 now works with Metal and CBCT example (MATLAB only!)

- If you have attenuation correction data in 1/cm, you can use those without scaling by setting `options.attIncm` to true

- Related to above, fixed Inveon attenuation in Python

- Small improvements to many examples

- Interactive Python sessions should now clear most GPU memory when using built-in reconstruction

- Added Barzilai-Borwain (BB) algorithm for Python as well (use `options.BB = True` to enable)

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
  - This means that testing will be limited for the other implementations, but potential issues will be fixed
  
### New features

- OMEGA can now be installed in Python through pip using `pip install omegatomo`
  - Examples are not included in this case and must be manually obtained from GitHub
  - Includes precompiled library files for Windows and Linux

- Added support for orthogonal distance ray-tracing projector with SPECT
  - Voxels are weighted by a Gaussian distribution sampled at the orthogonal distance, the variance of which is determined by the parallel distance between detector and voxel center

- Projector type 1 now has built-in support for SPECT parallel-hole and pinhole collimators
  - Variables for focal length: `options.colFxy` and `options.colFz`
  - Both take the value Inf for parallel-hole collimators
  - Both take the value 0 for pinhole collimators
  - Tested with SIMIND v8.0 built-in examples
  
- Added preliminary support for curved helical CT data
  - Functionality is more limited when compared to flat panel (CB)CT
  - Only projector type 4 is supported
  - Rotation of the image in the detector space (`options.offangle`) is NOT supported
  - Reconstruction speed is suboptimal as the code is not optimized
  - In general requires the XYZ-coordinates of both the source and detector as well as the rotation angles and the radius of the circle formed by the detector (`options.helicalRadius`, which should also be the distance from detector to focal spot)
  - Enable by setting `options.useHelical` to true
  - Only built-in reconstruction at the moment (projector operators are NOT supported)
  - See `helical_CT_mainExample` for an example

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
  
- FISTA acceleration can now be used with any algorithm
  - Enable with `options.FISTA_acceleration`
  - Can lead to divergence!
  - Implementation 2 only!

- Added another FISTA acceleration scheme
  - Select with `options.FISTAType`, use either 0 or 1
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
