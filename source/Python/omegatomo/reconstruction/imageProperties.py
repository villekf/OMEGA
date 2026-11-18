# -*- coding: utf-8 -*-
"""
Python port of MATLAB's save_image_properties.m: collects the reconstruction settings
into a dictionary (the second output of MATLAB's reconstructions_main).
"""

import numpy as np

# Port of recNames(0) (the algorithm flags), recNames(11) (full algorithm names, display only), recNames(1) (the prior flags)
# and recNames(10) (full prior names, display only). The lists are kept in the same order as in MATLAB's recNames.m; the full
# names are paired with the flags by position (zip), so both lists must have the same order.
_OS = ['OSEM', 'MRAMLA', 'RAMLA', 'ROSEM', 'RBI', 'DRAMA', 'COSEM', 'ECOSEM', 'ACOSEM', 'LSQR', 'CGLS', 'FISTA', 'FISTAL1', 'FDK']
_MAP = ['OSL_OSEM', 'BSREM', 'MBSREM', 'ROSEM_MAP', 'OSL_RBI', 'OSL_COSEM', 'PKMA', 'SPS', 'PDHG', 'PDHGKL', 'PDHGL1', 'CV', 'PDDY',
        'SART', 'ASD_POCS', 'SAGA', 'BB']
_ALGORITHM_FLAGS = _OS + _MAP
_ALGORITHM_NAMES = [
    'Ordered Subsets Expectation Maximization (OSEM)', 'Modified Row-action Maximum-Likelihood Algorithm (MRAMLA)',
    'Row-action Maximum-Likelihood Algorithm (RAMLA)', 'Relaxed Ordered Subsets Expectation Maximization (ROSEM)',
    'Rescaled Block-iterative algorithm (RBI)', 'Dynamic Row-action Maximum-Likelihood Algorithm (DRAMA)',
    'Complete-data Ordered Subsets Expectation Maximization (COSEM)',
    'Enhanced Complete-data Ordered Subsets Expectation Maximization (ECOSEM)',
    'Accelerated Complete-data Ordered Subsets Expectation Maximization (ACOSEM)', 'LSQR',
    'Conjugate-gradient least-squares (CGLS)', 'Fast Iterative Shrinkage Thresholding Algorithm (FISTA)',
    'Fast Iterative Shrinkage Thresholding Algorithm with L1 regularization (FISTAL1)', 'FDK/FBP algorithm',
    'One-step-late OSEM (OSL_OSEM)',
    'Block-sequential regularized expectation maximization (BSREM)',
    'Modified block-sequential regularized expectation maximization (MBSREM)', 'ROSEM MAP (ROSEM_MAP)', 'One-step-late RBI (OSL_RBI)',
    'One-step-late COSEM (OSL_COSEM)', 'Preconditioned Krasnoselski-Mann algorithm (PKMA)', 'Separable paraboloidal surrogates (SPS)',
    'Primal-dual hybrid gradient (PDHG)', 'Primal-dual hybrid gradient with Kullback-Leibler divergence (PDHGKL)',
    'Primal-dual hybrid gradient with L1 error (PDHGL1)', 'Condat-Vu algorithm', 'Primal-dual Davis-Yin algorithm',
    'Simultaneous Algebraic Reconstruction Technique', 'ASD-POCS', 'SAGA',
    'Barzilia-Borwein']
_PRIOR_FLAGS = ['MRP', 'quad', 'Huber', 'L', 'FMH', 'weighted_mean', 'TV', 'hyperbolic', 'AD', 'APLS', 'TGV', 'NLM', 'RDP', 'GGMRF',
                'ProxTV', 'custom']
_PRIOR_NAMES = [
    'Median Root', 'Quadratic', 'Huber', 'L-filter', 'FIR Median Hybrid', 'Weighted mean', 'Total Variation', 'Hyperbolic function',
    'Anisotropic Diffusion', 'Asymmetric Parallel Level Sets', 'Total Generalized Variation', 'Non-Local Means', 'Relative Difference',
    'General Gaussian Markov Random Field', 'Proximal Total variation', 'Custom']

_PRECOND_IMAGE = ['Diagonal', 'EM', 'IEM', 'Momentum', 'Normalized gradient', 'Filtering', 'Curvature']
_PRECOND_MEAS = ['Diagonal', 'Filtering']


def _py(value):
    """Converts NumPy scalars, 0-d arrays and arrays with one element to plain Python scalars. Other arrays are kept as NumPy
    arrays. Lists are converted element-wise."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return value.reshape(-1)[0].item()
        return value
    if isinstance(value, (list, tuple)):
        return [_py(v) for v in value]
    return value


def _get(options, name, default=None):
    return getattr(options, name, default)


def _flag(options, name):
    """Truth value of an options flag (False if it does not exist or is empty)."""
    value = getattr(options, name, False)
    if isinstance(value, np.ndarray):
        return bool(value.size > 0 and np.all(value))
    return bool(value)


def _num(value):
    """Plain Python float of a scalar or the first element of an array."""
    return float(np.asarray(value).reshape(-1)[0])


def _num2str(value):
    """MATLAB's num2str for a real scalar."""
    value = _num(value)
    if value == np.floor(value):
        return '%d' % value
    return '%.*g' % (min(max(int(np.floor(np.log10(abs(value)))) + 5, 5), 16), value)


def _projectorString(options):
    pt = int(_num(options.projector_type))
    if pt == 0:
        s = "Original Siddon's ray tracer"
    elif pt in (1, 11):
        if _flag(options, 'precompute_lor'):
            s = "Improved Siddon's ray tracer with 1 ray"
        else:
            nT = _num(options.n_rays_transaxial)
            nA = _num(options.n_rays_axial)
            ray = 'rays' if nT > 1 else 'ray'
            aray = 'rays' if nA > 1 else 'ray'
            s = ("Improved Siddon's algorithm used with " + _num2str(nT) + ' transaxial ' + ray + ' and ' + _num2str(nA) + ' axial ' +
                 aray + '.')
    elif pt in (2, 22):
        if _num(options.tube_width_z) == 0:
            s = 'Orthogonal distance-based ray tracer in 2D'
        else:
            s = 'Orthogonal distance-based ray tracer in 3D'
    elif pt in (3, 33):
        s = 'Volume-based ray tracer'
    else:
        names = {
            13: 'Improved Siddon for forward projection, volume-based ray tracer for backprojection',
            23: 'Orthogonal distance-based ray tracer for forward projection, volume-based ray tracer for backprojection',
            12: 'Improved Siddon for forward projection, orthogonal distance-based ray tracer for backprojection',
            32: 'Volume-based ray tracer for forward projection, orthogonal distance-based ray tracer for backprojection',
            31: 'Volume-based ray tracer for forward projection, improved Siddon for backprojection',
            21: 'Orthogonal distance-based ray tracer for forward projection, improved Siddon for backprojection',
            43: 'Interpolation-based ray tracer for forward projection, volume-based ray tracer for backprojection',
            42: 'Interpolation-based ray tracer for forward projection, orthogonal distance-based ray tracer for backprojection',
            34: 'Volume-based ray tracer for forward projection, interpolation-based ray tracer for backprojection',
            24: 'Orthogonal distance-based ray tracer for forward projection, interpolation-based ray tracerfor backprojection',
            4: 'Interpolation-based ray tracer',
            44: 'Interpolation-based ray tracer',
            14: 'Improved Siddon for forward projection, interpolation-based for backprojection',
            41: 'Improved Siddon for backprojection, interpolation-based for forward projection',
            15: 'Improved Siddon for forward projection, branchless distance-driven for backprojection',
            45: 'Interpolation-based for forward projection, branchless distance-driven for backprojection',
            54: 'Branchless distance-driven for forward projection, interpolation-based for backprojection',
            51: 'Branchless distance-driven for forward projection, improved Siddon for backprojection',
            5: 'Branchless distance-driven projector',
            55: 'Branchless distance-driven projector',
            6: 'SPECT rotation-based projector'}
        if pt in names:
            s = names[pt]
        else:
            s = 'Proj types ' + str(pt // 10) + ' and ' + str(pt - (pt // 10) * 10)
    return s + ' (' + str(pt) + ')'


def _implementationString(options):
    # The Python implementation is always the matrix-free ArrayFire one (implementation 2)
    if _flag(options, 'useCUDA'):
        return 'Matrix-free ArrayFire CUDA (2)'
    if _flag(options, 'useCPU'):
        return 'Matrix-free ArrayFire CPU (2)'
    return 'Matrix-free ArrayFire OpenCL (2)'


def imageProperties(options):
    """
    Collects the reconstruction settings into a dictionary (port of MATLAB's save_image_properties).

    The field names are the same as in MATLAB's image_properties struct. reconstructions_main stores the result in
    options.image_properties after the reconstruction. The values are those used in the reconstruction, i.e. they are read from the
    processed options.

    Parameters
    ----------
    options : class object
        OMEGA class object, after the options have been processed (as at the end of reconstructions_main).

    Returns
    -------
    ip : dict
        The image properties. Can be given to saveInterfile/saveMetaImage/saveImage as the options input.
    """
    o = options
    ip = {}
    isCT = _flag(o, 'CT')
    ip['Implementation'] = _implementationString(o)
    ip['deviceNumber'] = _py(_get(o, 'deviceNum', 0))
    ip['Projector'] = _projectorString(o)
    pt = int(_num(o.projector_type))
    if pt in (2, 22) and _num(o.tube_width_z) != 0:
        ip['orthogonalDistanceZ'] = _py(o.tube_width_z)
    if pt in (3, 33):
        ip['tubeRadius'] = _py(o.tube_radius)
        ip['voxelRadius'] = _py(o.voxel_radius)
    if _flag(o, 'use_psf'):
        if _flag(o, 'deblurring'):
            ip['PSF'] = 'PSF enabled with deblurring'
            ip['deblur_iterations'] = _py(o.deblur_iterations)
        else:
            ip['PSF'] = 'PSF enabled'
        ip['PSF_FWHM'] = _py(o.FWHM)
    else:
        ip['PSF'] = 'No PSF'
    ip['nIter'] = _py(o.Niter)
    ip['nSubsets'] = _py(o.subsets)
    if _num(_get(o, 'TOF_bins', 1)) > 1:
        ip['TOF'] = 'TOF enabled'
        ip['N_TOF_bins'] = _py(o.TOF_bins)
        ip['TOF_bin_width'] = _py(o.TOF_width)
        ip['TOF_FWHM'] = _py(o.TOF_FWHM)
    Nx = np.atleast_1d(np.asarray(o.Nx)).reshape(-1)
    Ny = np.atleast_1d(np.asarray(o.Ny)).reshape(-1)
    Nz = np.atleast_1d(np.asarray(o.Nz)).reshape(-1)
    FOVax = np.atleast_1d(np.asarray(o.FOVa_x)).reshape(-1)
    FOVay = np.atleast_1d(np.asarray(o.FOVa_y)).reshape(-1)
    axFOV = np.atleast_1d(np.asarray(o.axial_fov)).reshape(-1)
    # The original (non-extended) sizes of Python's options are only meaningful with the extended FOV (the defaults are 1/0)
    useEFOV = _flag(o, 'useEFOV')
    NxOrig = _num(o.NxOrig) if useEFOV and _flag(o, 'transaxialEFOV') else None
    NyOrig = _num(o.NyOrig) if useEFOV and _flag(o, 'transaxialEFOV') else None
    NzOrig = _num(o.NzOrig) if useEFOV and _flag(o, 'axialEFOV') else None
    FOVxOrig = _num(o.FOVxOrig) if _num(_get(o, 'FOVxOrig', 0)) > 0 else None
    FOVyOrig = _num(o.FOVyOrig) if _num(_get(o, 'FOVyOrig', 0)) > 0 else None
    axialFOVOrig = _num(o.axialFOVOrig) if _num(_get(o, 'axialFOVOrig', 0)) > 0 else None
    if not _flag(o, 'useMultiResolutionVolumes'):
        if NxOrig is not None and NxOrig < Nx[0]:
            ip['Nx'] = _py(NxOrig)
            ip['NxExtended'] = _py(Nx[0])
        else:
            ip['Nx'] = _py(Nx[0])
        if NyOrig is not None and NyOrig < Ny[0]:
            ip['Ny'] = _py(NyOrig)
            ip['NyExtended'] = _py(Ny[0])
        else:
            ip['Ny'] = _py(Ny[0])
        if NzOrig is not None and NzOrig < Nz[0]:
            ip['Nz'] = _py(NzOrig)
            ip['NzExtended'] = _py(Nz[0])
        else:
            ip['Nz'] = _py(Nz[0])
        if FOVxOrig is not None and FOVxOrig < FOVax[0]:
            ip['FOVx'] = _py(FOVxOrig)
            ip['FOVxExtended'] = _py(FOVax[0])
        else:
            ip['FOVx'] = _py(FOVax[0])
        if FOVyOrig is not None and FOVyOrig < FOVay[0]:
            ip['FOVy'] = _py(FOVyOrig)
            ip['FOVyExtended'] = _py(FOVay[0])
        else:
            ip['FOVy'] = _py(FOVay[0])
        if axialFOVOrig is not None and axialFOVOrig < axFOV[0]:
            ip['axialFOV'] = _py(axialFOVOrig)
            ip['axialFOVExtended'] = _py(axFOV[0])
        else:
            ip['axialFOV'] = _py(axFOV[0])
        ip['voxelSizeX'] = float(ip['FOVx']) / float(ip['Nx'])
        ip['voxelSizeY'] = float(ip['FOVy']) / float(ip['Ny'])
        ip['voxelSizeZ'] = float(ip['axialFOV']) / float(ip['Nz'])
    else:
        ip['Nx'] = _py(Nx[0])
        ip['Ny'] = _py(Ny[0])
        ip['Nz'] = _py(Nz[0])
        ip['NxExtended'] = Nx[1:]
        ip['NyExtended'] = Ny[1:]
        ip['NzExtended'] = Nz[1:]
        if FOVxOrig is not None and FOVxOrig < FOVax[0]:
            ip['FOVx'] = _py(FOVxOrig)
        else:
            ip['FOVx'] = _py(FOVax[0])
        ip['FOVxExtended'] = FOVax[1:]
        if FOVyOrig is not None and FOVyOrig < FOVay[0]:
            ip['FOVy'] = _py(FOVyOrig)
        else:
            ip['FOVy'] = _py(FOVay[0])
        ip['FOVyExtended'] = FOVay[1:]
        if axialFOVOrig is not None and axialFOVOrig < axFOV[0]:
            ip['axialFOV'] = _py(axialFOVOrig)
        else:
            ip['axialFOV'] = _py(axFOV[0])
        ip['axialFOVExtended'] = axFOV[1:]
        ip['voxelSizeX'] = float(ip['FOVx']) / float(Nx[0])
        ip['voxelSizeY'] = float(ip['FOVy']) / float(Ny[0])
        ip['voxelSizeZ'] = float(ip['axialFOV']) / float(Nz[0])
        ip['multiResolutionUsed'] = True
        ip['nMultiResolutionVolumes'] = _py(o.nMultiVolumes)
        ip['multiResolutionScalingValue'] = _py(o.multiResolutionScale)
        ip['voxelSizeXExtended'] = ip['FOVxExtended'] / Nx[1:].astype(np.float64)
        ip['voxelSizeYExtended'] = ip['FOVyExtended'] / Ny[1:].astype(np.float64)
        ip['voxelSizeZExtended'] = ip['axialFOVExtended'] / Nz[1:].astype(np.float64)
    ip['name'] = _get(o, 'name', '')
    ip['scannerName'] = _get(o, 'machine_name', '')
    ip['nTimeSteps'] = _py(_get(o, 'Nt', 1))
    ip['offsetAngle'] = _py(o.offangle)
    ip['flippedImage'] = _py(o.flip_image)
    ip['scatterCorrection'] = _py(o.scatter_correction or getattr(o, '_scatterPrecorrected', False))
    if isCT:
        ip['Binning'] = _py(o.binning)
        ip['detectorPitchX'] = _py(o.dPitchX)
        ip['detectorPitchY'] = _py(o.dPitchY)
        ip['nProjections'] = _py(o.nProjections)
        # Python has no extrapolated detector size (nColsDOrig/nRowsDOrig), the detector size is the one used
        ip['nColumnsProjections'] = _py(o.nColsD)
        ip['nRowsProjections'] = _py(o.nRowsD)
        ip['sourceOffsetRow'] = _py(o.sourceOffsetRow)
        ip['sourceOffsetCol'] = _py(o.sourceOffsetCol)
        ip['sourceToDetectorDistance'] = _py(o.sourceToDetector)
        ip['sourceToCenterRotation'] = _py(o.sourceToCRot)
        ip['nBedPositions'] = _py(o.nBed)
        ip['bedOffset'] = _py(o.bedOffset)
        if hasattr(o, 'flatFieldScaling'):
            ip['flatFieldScaling'] = _py(o.flatFieldScaling)
    else:
        ip['nSinogramAngles'] = _py(o.Nang)
        ip['nSinogramRadial'] = _py(o.Ndist)
        ip['nSinograms'] = _py(o.NSinos)
        ip['rawDataUsed'] = _py(o.use_raw_data)
        ip['startTime'] = _py(o.start)
        ip['endTime'] = _py(o.end)
        ip['nRings'] = _py(o.rings)
        if _flag(o, 'use_raw_data'):
            ip['ringDifference'] = _py(o.ring_difference_raw)
            if _num(o.sampling_raw) > 1:
                ip['samplingIncrease'] = _py(o.sampling_raw)
        else:
            ip['ringDifference'] = _py(o.ring_difference)
            if _num(o.sampling) > 1:
                ip['samplingIncrease'] = _py(o.sampling)
        ip['detectorsPerRing'] = _py(o.det_w_pseudo)
        ip['boreDiameter'] = _py(o.diameter)
        ip['blocksPerRing'] = _py(o.blocks_per_ring)
        ip['axialBlocks'] = _py(o.linear_multip)
        ip['transaxialBlocks'] = _py(o.transaxial_multip)
        ip['crystalsPerBlock'] = _py(o.cryst_per_block)
        ip['crystalsPerBlockAxial'] = _py(o.cryst_per_block_axial)
        ip['attenuationCorrection'] = _py(o.attenuation_correction)
        ip['arcCorrection'] = _py(o.arc_correction)
        ip['totalTime'] = _py(o.tot_time)
        ip['normalizationCorrection'] = _py(o.normalization_correction)
        ip['randomsCorrection'] = _py(o.randoms_correction or getattr(o, '_randomsPrecorrected', False))
        if _flag(o, 'randoms_correction') or _flag(o, '_randomsPrecorrected'):
            ip['randomsVarianceReduction'] = _py(o.variance_reduction)
            ip['randomsSmoothing'] = _py(o.randoms_smoothing)
        ip['gapFilling'] = _py(o.fill_sinogram_gaps)
        if _flag(o, 'reconstruct_trues'):
            ip['reconstructTrues'] = _py(o.reconstruct_trues)
        if _flag(o, 'reconstruct_scatter'):
            ip['reconstructScatter'] = _py(o.reconstruct_scatter)
        if hasattr(o, 'DOI'):
            ip['DOI'] = _py(o.DOI)
        if _flag(o, 'use_ASCII'):
            ip['dataType'] = 'Used ASCII GATE data'
        elif _flag(o, 'use_LMF'):
            ip['dataType'] = 'Used LMF GATE data'
        elif _flag(o, 'use_root'):
            ip['dataType'] = 'Used ROOT GATE data'
        elif _num(o.use_machine) == 1:
            ip['dataType'] = 'Used data obtained from a list-mode file'
        elif _num(o.use_machine) == 2:
            ip['dataType'] = 'Used scanner created sinogram data'
        elif _num(o.use_machine) == 3:
            ip['dataType'] = 'Used 32-bit list-mode data'
    if _num(o.subsets) > 1:
        ip['subsetType'] = _py(o.subsetType)
    ip['ordinaryPoissonReconstruction'] = bool(_flag(o, 'corrections_during_reconstruction') and
                                               (_flag(o, 'scatter_correction') or _flag(o, 'randoms_correction')))
    # Algorithm-specific parameters
    if _flag(o, 'ACOSEM') or _num(_get(o, 'OSL_COSEM', 0)) == 1:
        ip['hACOSEM'] = _py(o.h_ACOSEM)
    if _flag(o, 'RAMLA') or _flag(o, 'BSREM'):
        ip['relaxationParameter'] = _py(o.lambdaN)
    if _flag(o, 'MRAMLA') or _flag(o, 'MBSREM'):
        ip['relaxationParameter'] = _py(o.lambdaN)
    if _flag(o, 'ROSEM') or _flag(o, 'ROSEM_MAP'):
        ip['relaxationParameter'] = _py(o.lambdaN)
    if _flag(o, 'DRAMA'):
        ip['relaxationParameter'] = _py(o.lam_drama)
    if _flag(o, 'MRAMLA') or _flag(o, 'MBSREM'):
        ip['U'] = _py(o.U)
    if _flag(o, 'PKMA'):
        ip['rho_PKMA'] = _py(o.rho_PKMA)
        ip['delta_PKMA'] = _py(o.delta_PKMA)
        ip['relaxationParameter'] = _py(o.lambdaN)
        ip['MomentumParameterPKMA'] = _py(o.alpha_PKMA)
    for flagName, fullName in zip(_ALGORITHM_FLAGS, _ALGORITHM_NAMES):
        if _flag(o, flagName):
            ip['algorithm'] = fullName
    for flagName, fullName in zip(_PRIOR_FLAGS, _PRIOR_NAMES):
        if _flag(o, flagName):
            ip['prior'] = fullName
    ip['regularizationParameter'] = _py(o.beta)
    if _flag(o, 'quad'):
        ip['weights'] = _py(o.weights)
        ip['quadratic_prior_weights'] = _py(o.weights_quad)
    if _flag(o, 'Huber'):
        ip['weights_huber'] = _py(o.weights_huber)
        ip['huber_delta'] = _py(o.huber_delta)
    if _flag(o, 'L'):
        ip['L_weights'] = _py(o.a_L)
    if any(_flag(o, n) for n in ('L', 'quad', 'weighted_mean', 'FMH', 'MRP', 'NLM', 'GGMRF')):
        ip['NeighborhoodSizeX'] = _py(o.Ndx)
        ip['NeighborhoodSizeY'] = _py(o.Ndy)
        ip['NeighborhoodSizeZ'] = _py(o.Ndz)
    if _flag(o, 'FMH'):
        ip['fmh_center_weight'] = _py(o.fmh_center_weight)
        ip['fmh_weights'] = _py(o.fmh_weights)
    if _flag(o, 'weighted_mean'):
        ip['weighted_mean_center_weight'] = _py(o.weighted_center_weight)
        ip['weighted_mean_weights'] = _py(o.weighted_weights)
        ip['mean_type'] = _py(o.mean_type)
    if _flag(o, 'TV'):
        ip['TVsmoothing_beta'] = _py(o.TVsmoothing)
        ip['TV_use_anatomical'] = _py(o.TV_use_anatomical)
        tvType = int(_num(o.TVtype))
        if _flag(o, 'TV_use_anatomical'):
            ip['TV_reference_image'] = _py(o.TV_referenceImage)
            # options.T is a method of the options class (the transpose operator) unless the user has set the weight
            ip['T_weighting_parameter'] = _py(o.__dict__.get('T', 0.5))
            if tvType == 3:
                ip['C_original_image_weight'] = _py(o.C)
            if tvType == 1:
                ip['TVType'] = 'Anatomically weighted TV (TV type 1)'
            elif tvType == 2:
                ip['TVType'] = 'Anatomically weighted joint TV (TV type 2)'
            elif tvType == 3:
                ip['TVType'] = 'Anatomically weighted joint prior (TV type 3)'
        else:
            if tvType == 1:
                ip['TVType'] = 'Regular TV (TV type 1)'
            elif tvType == 2:
                ip['TVType'] = 'Regular TV (TV type 2)'
            elif tvType == 3:
                ip['TVType'] = 'Joint prior (hyperbolic function) (TV type 3)'
                ip['HyperbolicWeight'] = _py(o.SATVPhi)
            elif tvType == 4:
                ip['TVType'] = 'Smoothed anisotropic TV (Lange function) (TV type 4)'
                ip['LangePhi'] = _py(o.SATVPhi)
            elif tvType == 6:
                ip['TVType'] = 'Weighted TV (TV type 6)'
    if _flag(o, 'AD'):
        ip['TimeStepAD'] = _py(o.TimeStepAD)
        ip['KAD'] = _py(o.KAD)
        ip['NiterAD'] = _py(o.NiterAD)
        if _num(o.FluxType) == 1:
            ip['FluxTypeAD'] = 'exponential'
        else:
            ip['FluxTypeAD'] = 'quadratic'
        ip['DiffusionTypeAD'] = _py(o.DiffusionType)
    if _flag(o, 'APLS'):
        ip['APLSsmoothing_beta'] = _py(o.APLSsmoothing)
        ip['APLS_reference_image'] = _py(o.APLS_ref_image)
        ip['eta'] = _py(o.eta)
    if _flag(o, 'TGV'):
        ip['alpha0TGV'] = _py(o.alpha0TGV)
        ip['alpha1TGV'] = _py(o.alpha1TGV)
    if _flag(o, 'NLM'):
        ip['NLPatchSizeX'] = _py(o.Nlx)
        ip['NLPatchSizeY'] = _py(o.Nly)
        ip['NLPatchSizeZ'] = _py(o.Nlz)
        ip['NLFilterParameter'] = _py(o.NLMsigma)
        ip['NLM_gaussSTD'] = _py(o.NLM_gauss)
        ip['NLM_use_anatomical'] = _py(o.NLM_use_anatomical)
        if _flag(o, 'NLM_MRP'):
            ip['NLMType'] = 'NLM with NLM filtered image subtracted from the input image (NLM_MRP)'
        elif _flag(o, 'NLTV'):
            ip['NLMType'] = 'Non-local TV (NLTV)'
        elif _flag(o, 'NLRD'):
            ip['NLMType'] = 'Non-local relative difference prior (NLRDP)'
            ip['RDP_gamma'] = _py(o.RDP_gamma)
        elif _flag(o, 'NLLange'):
            ip['NLMType'] = 'Non-local Lange function prior (NLL)'
            ip['LangePhi'] = _py(o.SATVPhi)
        elif _flag(o, 'NLGM'):
            ip['NLMType'] = 'Non-local Geman-McClure prior (NLGM)'
            ip['GM_delta'] = _py(o.GM_delta)
        else:
            ip['NLMType'] = 'MRF-type NLM prior'
    precondMeas = np.asarray(o.precondTypeMeas).reshape(-1).astype(bool)
    precondImage = np.asarray(o.precondTypeImage).reshape(-1).astype(bool)
    if any(_flag(o, n) for n in ('PDHG', 'PDHGKL', 'PDHGL1', 'FISTA', 'FISTAL1')):
        ip['tauCP'] = _py(o.tauCP)
        ip['sigmaCP'] = _py(o.sigmaCP)
        ip['thetaCP'] = _py(o.thetaCP)
        if precondMeas[1]:
            ip['tauCPFilteredIterations'] = _py(o.tauCPFilt)
    if _flag(o, 'RDP'):
        ip['RDP_gamma'] = _py(o.RDP_gamma)
        ip['RDPIncludeCorners'] = _py(o.RDPIncludeCorners)
    if _flag(o, 'GGMRF'):
        ip['GGMRF_p'] = _py(o.GGMRF_p)
        ip['GGMRF_q'] = _py(o.GGMRF_q)
        ip['GGMRF_c'] = _py(o.GGMRF_c)
        ip['GGMRFweights'] = _py(o.weights_quad)
    if any(_flag(o, n) for n in ('PDHG', 'PDHGKL', 'PDHGL1')) and (_flag(o, 'ProxTV') or _flag(o, 'TGV')):
        ip['useL2Ball'] = _py(o.useL2Ball)
    if precondImage.any():
        ip['imageBasedPreconditioners'] = [_PRECOND_IMAGE[k] for k in range(min(len(precondImage), len(_PRECOND_IMAGE))) if precondImage[k]]
    if precondMeas[1]:
        ip['filteringIterations'] = _py(o.filteringIterations)
        ip['filterWindow'] = _py(o.filterWindow)
        if any(_flag(o, n) for n in ('PKMA', 'MBSREM', 'MRAMLA', 'SPS')):
            ip['relaxationParameterFiltering'] = _py(o.lambdaFiltered)
    if len(precondImage) > 3 and precondImage[3]:
        ip['Precond4gradV1'] = _py(o.gradV1)
        ip['Precond4gradV2'] = _py(o.gradV2)
        ip['Precond4gradInitIter'] = _py(o.gradInitIter)
    if precondMeas.any():
        ip['measurementBasedPreconditioners'] = [_PRECOND_MEAS[k] for k in range(min(len(precondMeas), len(_PRECOND_MEAS))) if precondMeas[k]]
    if _num(o.oOffsetZ) != 0 or _num(o.oOffsetX) != 0 or _num(o.oOffsetY) != 0:
        ip['objectOffsetsXYZ'] = np.array([_num(o.oOffsetX), _num(o.oOffsetY), _num(o.oOffsetZ)])
    ip['enforcePositivity'] = _py(o.enforcePositivity)
    ip['savedIterationNumbers'] = np.asarray(o.saveNIter).reshape(-1).astype(np.int64) + 1
    ip['ImagesWereUsed'] = _py(o.useImages)
    ip['MADWasUsed'] = _py(o.useMAD)
    if pt in (4, 44, 45):
        ip['interpolationLength'] = _py(o.dL)
    return ip
