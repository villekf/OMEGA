# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 17:40:13 2024

@author: Ville-Veikko Wettenhovi
"""

import numpy as np
import ctypes
import math
import os
from .coordinates import computePixelCenters
from .coordinates import computePixelSize
from .coordinates import computeProjectorScalingValues
from .coordinates import computeVoxelVolumes
from .indices import indexMaker
from .indices import formSubsetIndices
try:
    import arrayfire as af
except ModuleNotFoundError:
    print('ArrayFire package not found! ArrayFire features are not supported. You can install ArrayFire package with "pip install arrayfire".')
# from SimpleITK import GetSpacing

class projectorClass:
    # These parameters are either NumPy arrays or variables that are not needed in the C++ code
    x = np.empty(0, dtype = np.float32)
    y = np.empty(0, dtype = np.float32)
    z = np.empty(0, dtype = np.float32)
    x0 = np.empty(0, dtype = np.float32)
    Nx = 1
    Ny = 1
    Nz = 1
    index = np.empty(0, dtype = np.uint32)
    x_center = np.empty(0, dtype = np.float32)
    y_center = np.empty(0, dtype = np.float32)
    z_center = np.empty(0, dtype = np.float32)
    N = np.empty(0, dtype = np.uint64)
    saveNIter = np.empty(0, dtype = np.uint32)
    dx = np.empty(1, dtype = np.float32)
    dy = np.empty(1, dtype = np.float32)
    dz = np.empty(1, dtype = np.float32)
    dScaleX4 = np.empty(0, dtype = np.float32)
    dScaleY4 = np.empty(0, dtype = np.float32)
    dScaleZ4 = np.empty(0, dtype = np.float32)
    dSizeX = np.empty(0, dtype = np.float32)
    dSizeY = np.empty(0, dtype = np.float32)
    dScaleX = np.empty(0, dtype = np.float32)
    dScaleY = np.empty(0, dtype = np.float32)
    dScaleZ = np.empty(0, dtype = np.float32)
    kerroin = np.empty(0, dtype = np.float32)
    angles = np.empty(0, dtype = np.float32)
    blurPlanes = np.empty(0, dtype = np.float32)
    radiusPerProj = np.empty(0, dtype = np.float32)
    gFilter = np.empty(0, dtype = np.float32)
    filterIm = np.empty(0, dtype = np.float32)
    filter0 = np.empty(0, dtype = np.float32)
    filter2 = np.empty(0, dtype = np.float32)
    Ffilter = np.empty(0, dtype = np.float32)
    FWHM = np.zeros(3, dtype = np.float32)
    sourceOffsetCol = np.zeros(1, dtype = np.float32)
    sourceOffsetRow = np.zeros(1, dtype = np.float32)
    bedOffset = np.empty(0, dtype = np.float32)
    detOffsetRow = np.zeros(1, dtype = np.float32)
    detOffsetCol = np.zeros(1, dtype = np.float32)
    pitchRoll = np.empty(0, dtype = np.float32)
    OffsetLimit = np.empty(0, dtype = np.float32)
    uV = np.empty(0, dtype = np.float32)
    TOFCenter = np.empty(0, dtype = np.float32)
    precondTypeImage = np.full((7, 1), False)
    precondTypeMeas = np.full((2, 1), False)
    ScatterC = np.empty(0, dtype = np.float32)
    vaimennus = np.empty(0, dtype = np.float32)
    SinM = np.empty(0, dtype = np.float32)
    SinDelayed = np.empty(0, dtype = np.float32)
    corrVector = np.empty(0, dtype = np.float32)
    normalization = np.empty(0, dtype = np.float32)
    scatter_components = np.array([True, True, False, False])
    normalization_options = np.array([True, True, True, False])
    pseudot = np.empty(0, dtype = np.uint32)
    weights = np.empty(0, dtype = np.float32)
    weights_huber = np.empty(0, dtype = np.float32)
    a_L = np.empty(0, dtype = np.float32)
    fmh_weights = np.empty(0, dtype = np.float32)
    weighted_weights = np.empty(0, dtype = np.float32)
    weights_RDP = np.empty(0, dtype = np.float32)
    weights_quad = np.empty(0, dtype = np.float32)
    segment_table = np.empty(0, dtype = np.float32)
    eFOVIndices = np.empty(0, dtype = np.uint8)
    LL = np.empty(0, dtype = np.uint16)
    lambdaN = np.empty(0, dtype = np.float32)
    lambdaFiltered = np.empty(0, dtype = np.float32)
    lam_drama = np.empty(0, dtype = np.float32)
    alpha_PKMA = np.empty(0, dtype = np.float32)
    alphaPrecond = np.empty(0, dtype = np.float32)
    s = np.empty(0, dtype = np.float32)
    maskFP = np.empty(0, dtype = np.uint8)
    maskBP = np.empty(0, dtype = np.uint8)
    maskPrior = np.empty(0, dtype = np.uint8)
    gaussK = np.empty(0, dtype = np.float32)
    xy_index = np.empty(0, dtype=np.uint32)
    z_index = np.empty(0, dtype=np.uint16)
    tauCPFilt = 0.
    sigmaCP = 1.
    tauCP = 0.
    sigma2CP = 1.
    thetaCP = 1.
    normalization_self = np.array([True, True, True, True])
    machine_name = ""
    fpath = ""
    filterWindow = "hamming"
    sampling_interpolation_method = "linear"
    arc_interpolation = "linear"
    scatterPath = ""
    transaxialEFOV = False
    axialEFOV = False
    DOI = 0.
    cryst_per_block = np.empty(0, dtype = np.uint32)
    cryst_per_block_axial = cryst_per_block
    transaxial_multip = 1
    blocks_per_ring = 1
    diameter = 0.
    flip_image = False
    use_machine = 0
    use_ASCII = False
    use_binary = False
    use_LMF = False
    use_root = False
    only_sinos = False
    only_reconstructions = False
    usingLinearizedData = False
    no_data_load = False
    errorChecking = False
    no_data_load = False
    fill_sinogram_gaps = False
    ndist_side = 1
    store_raw_data = False
    sampling_raw = 1
    sampling = 1
    oneD_weights = False
    fmh_center_weight = 4.
    ring_difference_raw = 1
    ring_difference = 1
    obtain_trues = False
    reconstruct_trues = False
    store_scatter = False
    reconstruct_scatter = False
    store_randoms = False
    source = False
    weighted_center_weight = 2.
    offangle = 0.
    binning = 1
    sourceToCRot = 0.
    sourceToDetector = 1.
    nBed = 1
    TOF_FWHM = 0.
    TOF_width = 0.
    TOF_bins_used = 1
    span = 3
    cutoffFrequency = 1.
    normalFilterSigma = 0.25
    useZeroPadding = False
    oOffsetX = 0.
    oOffsetY = 0.
    oOffsetZ = 0.
    FOVa_x = 0.
    FOVa_y = 0.
    axial_fov = 0.
    FOVxOrig = 0.
    FOVyOrig = 0.
    axialFOVOrig = 0.
    tot_time = float("inf")
    start = 0.
    end = tot_time
    compute_normalization = False
    normalization_phantom_radius = float("inf")
    normalization_attenuation = 0.
    normalization_scatter_correction = False
    use_user_normalization = False
    arc_correction = False
    corrections_during_reconstruction = False
    multiResolutionScale = .25
    name = ""
    voxel_radius = 1.
    scatter_variance_reduction = False
    normalize_scatter = False
    scatter_smoothing = False
    variance_reduction = False
    randoms_smoothing = False
    ringGaps = np.empty(0, dtype=np.int32)
    precompute_lor = False
    sigmaZ = -1.
    sigmaXY = -1.
    colL = 0.
    colR = 0.
    colD = 0.
    iR = 0.
    implementation = 2
    rotateAttImage = 0.
    flipAttImageXY = False
    flipAttImageZ = False
    NSinos = 1
    TotSinos = NSinos
    V = np.empty(0, dtype = np.float32)
    uu = 0
    ub = 0
    det_w_pseudo = 1
    startAngle = 0.
    nHeads = 1
    angleIncrement = 0.
    flatFieldScaling = 1.
    CT_attenuation = True
    scatter_correction = 0
    useCUDA = False
    useCPU = False
    NxFull = 1
    NyFull = 1
    NzFull = 1
    beta_drama = 0.
    alpha_drama = 0.
    beta0_drama = 0.
    beta = 0.
    rho_PKMA = .45
    delta_PKMA = 100.
    delta2_PKMA = 100.
    empty_weight = True
    projector_type = 11
    CT = False
    PET = False
    SPECT = False
    saveSens = True
    use_raw_data = 0
    dPitchX = 0.
    dPitchY = 0.
    subsets = 1
    subsetType = 8
    useMaskFP = False
    useMaskBP = False
    offsetCorrection = False
    tube_width_z = 0.
    tube_width_xy = 0.
    use_psf = False
    save_iter = False
    deblurring = False
    use_64bit_atomics = False
    n_rays_transaxial = 1
    n_rays_axial = 1
    RDPIncludeCorners = False
    meanFP = False
    meanBP = False
    Niter = 1
    filteringIterations = 0
    nLayers = 1
    powerIterations = 20
    deviceNum = 0
    platform = 0
    derivativeType = 0
    enforcePositivity = True
    gradV1 = 0.5
    gradV2 = 2.5
    gradInitIter = subsets
    gradLastIter = gradInitIter
    pitch = False
    compute_sensitivity_image = False
    listmode = 0
    nProjections = 1
    Ndist = 1
    Nang = 1
    NSinos = 1
    dL = 0
    epps = 1e-6
    cr_p = 1.
    cr_pz = 1.
    use_32bit_atomics = False
    verbose = 0
    partitions = 1
    Nt = 1
    randoms_correction = 0
    additionalCorrection = 0
    attenuation_correction = 0
    normalization_correction = 0
    global_correction_factor = 1.
    rings = 1
    NxOrig = 1
    NyOrig = 1
    NzOrig = 1
    NxPrior = 1
    NyPrior = 1
    NzPrior = 1
    det_per_ring = Nang * Ndist
    h_ACOSEM = 2.
    U = 10000.
    Ndx = 1
    Ndy = 1
    Ndz = 1
    mean_type = 4
    TVsmoothing = 1e-2
    TV_use_anatomical = False
    TVtype = 1
    tau = 0.
    B = 0.01
    C = 1.
    SATVPhi = 1.
    FluxType = 1
    DiffusionType = 2
    eta = 1e-5
    APLSsmoothing = 1e-5
    Nlx = 2
    Nly = 2
    Nlz = 2
    TOF = False
    TOF_bins = 1
    NiterAD = 1
    sigma_x = 1
    w_sum = 1.
    RDP_gamma = 1.
    NLMsigma = 1.
    NLAdaptiveConstant = 1.
    TimeStepAD = 1.
    KAD = 1.
    huber_delta = 1.
    NLM_use_anatomical = False
    NLAdaptive = False
    NLTV = False
    NLRD = False
    NLLange = False
    NLGGMRF = False
    NLM_MRP = False
    med_no_norm = False
    alpha0TGV = 0.
    alpha1TGV = 0.
    useL2Ball = True
    useMAD = True
    useImages = True
    useEFOV = False
    useExtrapolation = False
    CTAttenuation = True
    flat = 0.
    use2DTGV = False
    hyperbolicDelta = 1.
    useMultiResolutionVolumes = False
    nMultiVolumes = 0
    relaxationScaling = False
    computeRelaxationParameters = False
    PDAdaptiveType = 0
    storeFP = False
    nRowsD = Ndist
    nColsD = Nang
    Nf = 0
    tube_radius = 1.
    OSEM = False
    LSQR = False
    CGLS = False
    SART = False
    FISTA = False
    FISTAL1 = False
    MRAMLA = False
    RAMLA = False
    ROSEM = False
    RBI = False
    DRAMA = False
    COSEM = False
    ECOSEM = False
    ACOSEM = False
    OSL_OSEM = False
    MBSREM = False
    BSREM = False
    ROSEM_MAP = False
    OSL_RBI = False
    OSL_COSEM = 0
    PKMA = False
    SPS = False
    PDHG = False
    PDHGKL = False
    PDHGL1 = False
    PDDY = False
    CV = False
    ASD_POCS = False
    FDK = False
    SAGA = False
    MRP = False
    quad = False
    Huber = False
    L = False
    FMH = False
    weighted_mean = False
    TV = False
    hyperbolic = False
    AD = False
    APLS = False
    TGV = False
    NLM = False
    RDP = False
    GGMRF = False
    ProxTV = False
    ProxRDP = False
    ProxNLM = False
    MAP = False
    custom = False
    Vmax = 0.
    bmax = 0.
    bmin = 0.
    dSizeZBP = 0.
    dSizeXBP = 0.
    inffi = 0
    headerDir = ''
    deblur_iterations = 0
    g_dim_x = 0
    g_dim_y = 0
    g_dim_z = 0
    global_factor = 1.
    GGMRF_p = 0.
    GGMRF_q = 0.
    GGMRF_c = 0.
    orthTransaxial = False
    orthAxial = False
    referenceImage = ''
    APLS_ref_image = ''
    TV_referenceImage = ''
    NLM_referenceImage = ''
    RDP_referenceImage = ''
    storeMultiResolution = False
    extrapLength = 0.2
    axialExtrapolation = False
    transaxialExtrapolation = False
    NLM_gauss = 0.7
    TOF_noise_FWHM = 0.
    TOF_offset = 0.
    gFSize = None
    trans = False
    subset = 0
    largeDim = False
    loadTOF = True
    useAF = False
    useTorch = False
    useCuPy = False
    dualLayerSubmodule = False
    storeResidual = False
    useFDKWeights = True
    useIndexBasedReconstruction = False
    trIndex = np.empty(0, dtype = np.uint16)
    axIndex = np.empty(0, dtype = np.uint16)
    POCS_alpha = 0.2
    POCS_rMax = 0.95
    POCS_alphaRed = 0.95
    POCSepps = 1e-4
    POCS_NgradIter = 20
    FISTA_acceleration = False
    RDP_use_anatomical = False
    nRays = 1
    flipImageX = False
    flipImageY = False
    flipImageZ = False
    crXY = 1
    rayShiftsDetector = np.empty(0, dtype=np.float32)
    rayShiftsSource = np.empty(0, dtype=np.float32)
    CORtoDetectorSurface = 0
    homeAngles = np.empty(0, dtype = np.float32)
    swivelAngles = np.empty(0, dtype = np.float32)
    eFOVLength = 0.4
    FISTAType = 0
    maskFPZ = 1
    maskBPZ = 1
    stochasticSubsetSelection = False
    useTotLength = True
    TOFIndices = np.empty(0, dtype = np.uint8)

    def __init__(self):
        # C-struct
        self.param = self.parameters()
    def addProjector(self):
        if self.OSL_OSEM or self.MBSREM or self.ROSEM_MAP or self.OSL_RBI or self.OSL_COSEM > 0 or self.PKMA or self.SPS or self.PDHG or self.PDHGKL or self.PDHGL1 or self.PDDY or self.CV:
            self.MAP = True
        self.OMEGAErrorCheck()
        if hasattr(self, 'dPitch') and self.dPitch > 0 and self.dPitchX == 0.:
            self.dPitchX = self.dPitch
        if hasattr(self, 'dPitch') and self.dPitch > 0 and self.dPitchY == 0.:
            self.dPitchY = self.dPitch
        if self.CT:
            from .detcoord import setCTCoordinates
            if hasattr(self,'NSinos') and (not hasattr(self,'nProjections') or self.nProjections == 0):
                self.nProjections = self.NSinos
            else:
                self.NSinos = self.nProjections
            self.NSinos = self.nProjections
            self.TotSinos = self.nProjections
            self.span = 3
            if hasattr(self,'Ndist') and (not hasattr(self,'nRowsD') or self.nRowsD == 0):
                self.nRowsD = self.Ndist
            else:
                self.Ndist = self.nRowsD
            if hasattr(self,'Nang') and (not hasattr(self,'nColsD') or self.nColsD == 0):
                self.nColsD = self.Nang
            else:
                self.Nang = self.nColsD
            self.Ndist = self.nRowsD
            self.Nang = self.nColsD
            self.use_raw_data = 0
            self.attenuation_correction = 0
            self.normalization_correction = 0
            self.compute_normalization = False
            self.rings = self.nColsD * self.nBed
            self.det_per_ring = self.nRowsD * self.nProjections
            self.arc_correction = False
            self.det_w_pseudo = self.det_per_ring
            self.cryst_per_block = self.nColsD * self.nRowsD
            if hasattr(self, 'dPitch') and self.dPitch > 0:
                self.cr_pz = self.dPitch
            elif hasattr(self, 'dPitchY'):
                self.cr_pz = self.dPitchY
            elif hasattr(self, 'cr_pz'):
                self.dPitchY = self.cr_pz
            if hasattr(self, 'dPitch') and self.dPitch > 0:
                self.cr_p = self.dPitch
            elif hasattr(self, 'dPitchX'):
                self.cr_p = self.dPitchX
            elif hasattr(self, 'cr_p'):
                self.dPitchX = self.cr_p
            if np.max(np.abs(self.angles.flatten())) > 6. * np.pi:
                self.angles = self.angles * (np.pi / 180.)
            if self.offangle > 0:
                self.angles = self.angles + self.offangle
            setCTCoordinates(self)
        if self.SPECT:
            self.NSinos = self.nProjections
            self.TotSinos = self.nProjections
            self.dPitchX = self.cr_p
            self.dPitchY = self.cr_pz
            self.det_per_ring = self.nRowsD * self.nProjections
            self.arc_correction = False
            self.det_w_pseudo = self.det_per_ring
            self.cryst_per_block = self.nColsD * self.nRowsD
            self.arc_correction = False
            self.span = 3
            self.Ndist = self.nRowsD
            self.Nang = self.nColsD
            self.use_raw_data = 0
        self.TOF = self.TOF_bins > 1 and (self.projector_type == 1 or self.projector_type == 11 or self.projector_type == 3 or self.projector_type == 33 
                                                      or self.projector_type == 13 or self.projector_type == 31 or self.projector_type == 4 or self.projector_type == 14 or self.projector_type == 41 
                                                      or self.projector_type == 44)
        if self.span == 1:
            self.TotSinos = self.rings**2
            self.NSinos = self.TotSinos
        if (self.CT == False and self.SPECT == False) and ((self.subsetType > 7 and self.subsets > 1) or self.subsets == 1):
            self.nProjections = self.NSinos
            self.PET = True
        else:
            self.PET = False
            self.nProjections = self.NSinos
        if self.TOF:
            if self.TOF_bins_used != self.TOF_bins:
                # self.SinM = np.sum(self.SinM,3)
                self.sigma_x = 0.
                # self.TOF = False
            else:
                c = 2.99792458e11
                self.sigma_x = (c*self.TOF_FWHM/2.) / (2. * math.sqrt(2. * math.log(2.)))
                edges_user = np.linspace(-self.TOF_width * self.TOF_bins/2, self.TOF_width * self.TOF_bins / 2, self.TOF_bins + 1, dtype=np.float32)
                edges_user = edges_user[0:-1] + self.TOF_width/2.
                self.TOFCenter = np.zeros(np.shape(edges_user),dtype = np.float32, order='F')
                self.TOFCenter[0] = edges_user[math.floor(np.size(edges_user)/2)]
                self.TOFCenter[1::2] = edges_user[math.floor(np.size(edges_user)/2) + 1:]
                self.TOFCenter[2::2] = edges_user[math.floor(np.size(edges_user)/2) - 1:  : -1]
                if self.TOF_offset > 0:
                    self.TOFCenter = self.TOFCenter + self.TOF_offset
                self.TOFCenter = -self.TOFCenter * c / 2.
        else:
            self.sigma_x = 0.
        if self.maskFP.size > 1 and ((not(self.maskFP.size == (self.nRowsD * self.nColsD)) and not(self.maskFP.size == (self.nRowsD * self.nColsD * self.nProjections)) 
                                      and (self.CT == True or self.SPECT == True)) or (not(self.maskFP.size == (self.Nang * self.Ndist)) and not(self.maskFP.size == (self.Nang * self.Ndist * self.NSinos)) and self.CT == False)):
            if self.CT == True or self.SPECT == True:
                raise ValueError('Incorrect size for the forward projection mask! Must be the size of a single projection image [' + str(self.nRowsD) + ' ' + str(self.nColsD) + ']  or full stack of [' + str(self.nRowsD) + ' ' + str(self.nColsD) + ' ' + str(self.nProjections) + ']')
            else:
                raise ValueError('Incorrect size for the forward projection mask! Must be the size of a single sinogram image [' + str(self.Nang) + ' ' + str(self.Ndist) + '] or 3D stack [' + str(self.Nang) + ' ' + str(self.Ndist) + ' ' + str(self.NSinos) + ']')
        elif self.maskFP.size > 1 and self.maskFP.size == (self.nRowsD * self.nColsD):
            self.useMaskFP = True
            if (self.maskFP.ndim == 3):
                self.maskFPZ = self.maskFP.shape[2]
        else:
            self.useMaskFP = False
        
        if self.maskBP.size > 1 and not(self.maskBP.size == self.Nx * self.Ny) and not(self.maskBP.size == self.Nx * self.Ny * self.Nz):
            raise ValueError('Incorrect size for the backward projection mask! Must be the size of a single image [' + str(self.Nx) + ' ' + str(self.Ny) + '] or 3D stack [' + str(self.Nx) + ' ' + str(self.Ny) + ' ' + str(self.Nz) + ']')
        elif self.maskBP.size == self.Nx * self.Ny:
            self.useMaskBP = True
            if (self.maskBP.ndim == 3):
                self.maskBPZ = self.maskBP.shape[2]
        else:
            self.useMaskBP = False
        list_mode_format = False
        rings = self.rings
        if self.use_raw_data and self.x.size > 1:
            det_per_ring = self.x.size
        else:
            det_per_ring = self.det_per_ring
        
        if self.use_raw_data:
            rings = rings - np.sum(self.pseudot)
            self.rings = rings
            # self.detectors = det_per_ring * rings;
        
        temp = self.pseudot
        if isinstance(temp, int):
            if temp > 0:
                if isinstance(self.cryst_per_block, np.ndarray):
                    self.pseudot = np.array(self.cryst_per_block[0].item() + 1,dtype=np.uint32)
                else:
                    self.pseudot = np.array(self.cryst_per_block + 1,dtype=np.uint32)
            else:
                self.pseudot = np.empty(0, dtype = np.uint32)
        elif temp is None:
            self.pseudot = np.empty(0, dtype = np.uint32)
        elif isinstance(temp, np.ndarray):
            if len(temp) > 0 and np.sum(temp) > 0:
                self.pseudot = np.zeros(temp, dtype=np.uint32)
                for kk in range(1,temp + 1):
                    if isinstance(self.cryst_per_block, np.ndarray):
                        self.pseudot[kk - 1] = np.array(self.cryst_per_block[0].item() + 1,dtype=np.uint32) * kk
                    else:
                        self.pseudot[kk - 1] = np.array(self.cryst_per_block + 1,dtype=np.uint32) * kk
        else:
            self.pseudot = np.array(0,dtype=np.uint32)
        # elif np.sum(temp) == 0 and temp.size > 0:
        #     self.pseudot = np.empty(0, dtype = np.uint32)
        # Whether list-mode or sinogram/raw data is used
        if isinstance(self.x, np.ndarray) and self.x.size > 0 and (self.x.size // 2 == self.SinM.size or self.x.size // 6 == self.SinM.size):
            det_per_ring = self.SinM.size
            self.Nang = 1
            self.Ndist = 1
            self.NSinos = det_per_ring
            self.TotSinos = self.NSinos
            list_mode_format = True
            self.listmode = 1
        elif self.useIndexBasedReconstruction:
            det_per_ring = self.trIndex.size // 2
            self.Nang = 1
            self.Ndist = 1
            self.NSinos = det_per_ring
            self.TotSinos = self.NSinos
            list_mode_format = True
            self.listmode = 1
        else:
            # Compute the indices for the subsets used.
            # For Sinogram data, six different methods to select the subsets are
            # available. For raw data, three methods are available.
            self.listmode = 0
        if self.listmode and self.subsets > 1 and not(self.subsetType == 1) and not(self.subsetType == 3):
            print('Only subset types 0, 1, and 3 are supported with list-mode data! Switching to subset type 0.')
            self.subsetType = 0
        indexMaker(self)
        self.setUpCorrections()
        self.x0 = self.x0.ravel('F')
        
        # Coordinates of the detectors
        if not(self.projector_type == 6):
            if self.listmode == False:
                if self.SPECT:
                    from .detcoord import getCoordinatesSPECT
                    x_det, z_det = getCoordinatesSPECT(self)
                    y = 0
                else:
                    from .detcoord import getCoordinates
                    x_det, y, z_det = getCoordinates(self)
            else:
                if not self.useIndexBasedReconstruction:
                    if self.x.shape[0] == 2:
                        if self.x.flags.f_contiguous:
                            self.x = np.row_stack((self.x[0,:], self.y[0,:], self.z[0,:], self.x[1,:], self.y[1,:], self.z[1,:]))
                        else:
                            self.x = np.asfortranarray(np.row_stack((np.self.x[0,:], self.y[0,:], self.z[0,:], self.x[1,:], self.y[1,:], self.z[1,:])))
                    elif self.x.shape[1] == 2:
                        if self.x.flags.f_contiguous:
                            self.x = np.row_stack((self.x[:,0].T(), self.y[:,0].T(), self.z[:,0].T(), self.x[:,1].T(), self.y[:,1].T(), self.z[:,1].T()))
                        else:
                            self.x = np.asfortranarray(np.row_stack((self.x[:,0].T(), self.y[:,0].T(), self.z[:,0].T(), self.x[:,1].T(), self.y[:,1].T(), self.z[:,1].T())))
                # y = 0
                x_det = 0
                z_det = 0
        else:
            # y = 0
            x_det = 0
            z_det = 0

        if self.use_raw_data == True:
            if list_mode_format == True:
                size_x = self.x.size // 6
            else:
                size_x = x_det.size
        else:
            if list_mode_format == True:
                if not self.useIndexBasedReconstruction:
                    size_x = self.x.size // 6
                else:
                    size_x = self.x.size // 2
            else:
                size_x = self.Ndist
            if self.sampling > 1:
                size_x = size_x * self.sampling
        if self.CT == True or self.projector_type == 6:
            size_x = self.nRowsD
            if self.listmode == True:
                size_x = self.x.size // 6
        else:
            if self.SPECT == True:
                self.dPitch = self.crXY
                self.dPitchY = self.crXY
                self.dPitchX = self.crXY
                self.cr_p = self.crXY
            else:
                self.dPitch = self.cr_p
                self.dPitchY = self.cr_p
                self.dPitchX = self.cr_pz
            self.nProjections = self.NSinos
            self.nRowsD = self.Ndist
            self.nColsD = self.Nang

        # self.size_x = size_x
        # self.totMeas = self.nColsD * self.nRowsD * self.nProjections
        self.nMeas = np.insert(np.cumsum(self.nMeas),0, 0)
        if not isinstance(self.Nx, np.ndarray):
            self.Nx = np.array(self.Nx, dtype=np.uint32, ndmin=1)
        if not isinstance(self.Ny, np.ndarray):
            self.Ny = np.array(self.Ny, dtype=np.uint32, ndmin=1)
        if not isinstance(self.Nz, np.ndarray):
            self.Nz = np.array(self.Nz, dtype=np.uint32, ndmin=1)
        xx, yy, zz = computePixelSize(self)
        formSubsetIndices(self)
        if (self.CT or self.PET or (self.SPECT and not(self.projector_type == 6))) and self.listmode == 0:
            if self.subsetType >= 8 and self.subsets > 1 and not self.FDK:
                if self.CT:
                    x_det = np.reshape(x_det, (self.nProjections, 6))
                    x_det = x_det[self.index,:]
                    x_det = x_det.ravel('C')
                    if self.pitch:
                        z_det = np.reshape(z_det, (self.nProjections, 6))
                        z_det = z_det[self.index,:]
                        z_det = z_det.ravel('C')
                    else:
                        z_det = np.reshape(z_det, (self.nProjections, 2))
                        z_det = z_det[self.index,:]
                        z_det = z_det.ravel('C')
                elif self.SPECT:
                    x_det = x_det[:,self.index]
                    x_det = x_det.ravel('F')
                    z_det = z_det[:,self.index]
                    z_det = z_det.ravel('F')
                else:
                    z_det = np.reshape(z_det, (self.nProjections, -1))
                    z_det = z_det[self.index,:]
                    z_det = z_det.ravel('C')
                if self.CT:
                    self.uV = self.uV[self.index,:]
        if self.listmode == 0:
            self.x = x_det
            self.z = z_det
        computePixelCenters(self, xx, yy, zz)
        computeVoxelVolumes(self)
        if (self.projector_type == 4 or self.projector_type == 5 or self.projector_type == 14 or self.projector_type == 41 
                    or self.projector_type == 15 or self.projector_type == 45 or self.projector_type == 54 or self.projector_type == 51 
                    or self.projector_type == 42 or self.projector_type == 43 or self.projector_type == 24 or self.projector_type == 34):
            computeProjectorScalingValues(self)
        if self.offsetCorrection and self.subsets > 1:
            self.OffsetLimit = self.OffsetLimit[self.index];

        if self.projector_type == 6:
            from .detcoord import SPECTParameters
            SPECTParameters(self)
            if self.subsets > 1 and (self.subsetType == 8 or self.subsetType == 9 or self.subsetType == 10 or self.subsetType == 11):
                self.angles = self.angles[self.index];
                self.blurPlanes = self.blurPlanes[self.index];
                self.gFilter = self.gFilter[:,:,:,self.index];
            # self.gFilter = self.gFilter.ravel('F').astype(dtype=np.float32)
        ## This part is used when the observation matrix is calculated on-the-fly

        if self.nMeas.size == 1:
            self.nMeas = np.insert(self.nMeas, 0, 0)

        if self.subsets > 1:
            self.subset = 0
        if self.subsetType >= 8 or self.subsets == 1:
            kerroin = self.nColsD * self.nRowsD
        else:
            kerroin = 1
        self.nMeasSubset = np.zeros((self.subsets, 1), dtype = np.int64);
        self.nProjSubset = np.zeros((self.subsets, 1), dtype = np.int64);
        self.nTotMeas = self.nMeas * kerroin
        for kk in range(self.subsets):
            self.nMeasSubset[kk] = self.nMeas[kk + 1] * kerroin - self.nMeas[kk] * kerroin
            self.nProjSubset[kk] = self.nMeas[kk + 1] - self.nMeas[kk]
        if self.listmode == 1:
            self.x = self.x.astype(dtype=np.float32)
            if self.x.flags.f_contiguous:
                self.x = self.x.ravel('F')
        # Compute PSF kernel
        self.PSFKernel()
        self.N = self.Nx.astype(np.uint64) * self.Ny.astype(np.uint64) * self.Nz.astype(np.uint64)
        if not isinstance(self.FOVa_x, np.ndarray):
            self.FOVa_x = np.array(self.FOVa_x, dtype=np.float32, ndmin=1)
        if not isinstance(self.FOVa_y, np.ndarray):
            self.FOVa_y = np.array(self.FOVa_y, dtype=np.float32, ndmin=1)
        if not isinstance(self.axial_fov, np.ndarray):
            self.axial_fov = np.array(self.axial_fov, dtype=np.float32, ndmin=1)
        if self.projector_type in [2, 3, 22, 33]:
            if self.projector_type in [3, 33]:
                self.orthTransaxial = True
            elif (self.projector_type in [2, 22]) and (self.tube_width_xy > 0):
                self.orthTransaxial = True
            else:
                self.orthTransaxial = False
        if self.projector_type in [2, 3, 22, 33]:
            if self.projector_type in [3, 33]:
                self.orthAxial = True
            elif (self.projector_type in [2, 22]) and (self.tube_width_z > 0):
                self.orthAxial = True
            else:
                self.orthAxial = False
        if self.use_32bit_atomics and self.use_64bit_atomics:
            self.use_32bit_atomics = False
        self.x0 = self.x0.astype(dtype=np.float32)
        if isinstance(self.x, int):
            self.x = np.zeros(1, dtype=np.float32)
            self.z = np.zeros(1, dtype=np.float32)
        else:
            self.x = self.x.astype(dtype=np.float32)
            self.z = self.z.astype(dtype=np.float32)
        if isinstance(self.partitions, np.ndarray):
            if self.partitions.size >= 1:
                self.Nt = self.partitions.size
            # else:
            #     self.Nt = self.partitions[0].item()
        else:
            self.Nt = self.partitions
        if self.listmode and self.randoms_correction:
            self.randoms_correction = False
        
        
    def OMEGAErrorCheck(self):
        if not self.CT and not self.SPECT and (self.FOVa_x >= self.diameter or self.FOVa_y >= self.diameter) and self.diameter > 0:
            raise ValueError(f"Transaxial FOV is larger than the scanner diameter ({self.diameter})!")
        if not self.CT and not self.SPECT and self.axial_fov < (self.rings * self.cr_pz - self.cr_pz):
            raise ValueError("Axial FOV is too small, crystal ring(s) on the boundary have no slices!")
        
        if not(self.PDHG or self.PDHGKL or self.PDHGL1 or self.PDDY or self.PKMA or self.FISTA or self.FISTAL1 or self.MBSREM or self.SPS or self.MRAMLA) and any(self.precondTypeImage):
            print("Image-based preconditioning selected, but the selected algorithm(s) do not support preconditioning. No preconditioning will be performed.")
            print("Supported algorithms are: MBSREM, MRAMLA, PKMA, SPS, PDHG, PDHGL1, PDHGKL, FISTA, FISTAL1, PDDY")
            self.precondTypeImage = np.array([False, False, False, False, False, False])
        
        if np.sum(self.precondTypeImage[0:3]) > 1:
            raise ValueError("Only one of the first 3 image-based preconditioners can be selected at a time!")
        
        if not(self.PDHG or self.PDHGKL or self.PDHGL1 or self.PDDY or self.PKMA or self.FISTA or self.FISTAL1 or self.MBSREM or self.SPS or self.MRAMLA) and any(self.precondTypeMeas):
            print("Measurement-based preconditioning selected, but the selected algorithm does not support preconditioning. No preconditioning will be performed.")
            print("Supported algorithms are: MBSREM, MRAMLA, PKMA, SPS, PDHG, PDHGL1, PDHGKL, FISTA, FISTAL1, PDDY")
            self.precondTypeMeas = np.array([False, False])        
        
        if not self.CT and not self.SPECT and self.span > self.ring_difference and self.NSinos > 1 and not self.use_raw_data:
            raise ValueError(f"Span value cannot be larger than ring difference ({self.ring_difference})!")
        
        if not self.CT and not self.SPECT and (self.span % 2 == 0 or self.span <= 0) and not self.use_raw_data:
            raise ValueError("Span value has to be odd and positive!")
        
        if not self.CT and not self.SPECT and self.ring_difference >= self.rings and not self.use_raw_data:
            print(f"Ring difference can be at most {self.rings - 1}. Setting it to the maximum possible.")
            self.ring_difference = self.rings - 1
        
        if not self.CT and not self.SPECT and self.ring_difference < 0 and not self.use_raw_data:
            raise ValueError("Ring difference has to be at least 0!")
        
        if not self.CT and not self.SPECT and self.Nang > self.det_w_pseudo / 2 and not self.use_raw_data and self.Nang > 1:
            raise ValueError(f"Number of sinogram angles can be at most the number of detectors per ring divided by two ({self.det_w_pseudo / 2})!")
        
        if not self.CT and not self.SPECT and self.TotSinos < self.NSinos and not self.use_raw_data:
            raise ValueError(f"The numnber of sinograms used ({self.NSinos}) is larger than the total number of sinograms ({self.TotSinos})!")
        
        if not self.CT and not self.SPECT and (self.ndist_side > 1 and self.Ndist % 2 == 0 or self.ndist_side < -1 and self.Ndist % 2 == 0) and not self.use_raw_data:
            raise ValueError("ndist_side can be either 1 or -1!")
        
        if not self.CT and not self.SPECT and self.ndist_side == 0 and self.Ndist % 2 == 0 and not self.use_raw_data:
            raise ValueError("ndist_side cannot be 0 when Ndist is even!")
        if self.useIndexBasedReconstruction and self.projector_type > 3:
            raise ValueError('Index-based recpnstruction only supports projector types 1-3!')
        
        # if not self.CT and not self.SPECT and ((self.sampling % 2 > 0 and self.sampling != 1) or self.sampling < 0) and not self.use_raw_data:
        #     raise ValueError("Sampling rate has to be divisible by two and positive or one!")
        
        # if not self.CT and not self.SPECT and ((self.sampling_raw % 2 > 0 and self.sampling_raw != 1) or self.sampling_raw < 0) and self.use_raw_data:
        #     raise ValueError("Sampling rate has to be divisible by two and positive or one!")
        
        # if (self.sampling > 1 and not self.use_raw_data) or (self.sampling_raw > 1 and self.use_raw_data) and self.implementation == 1:
        #     print("Increased sampling rate is not supported with implementation 1. Using sample rate of 1.")
        #     self.sampling = 1
        
        # if self.SPECT and self.implementation == 1:
        #     raise ValueError("Implementation 1 is not supported with SPECT data.")
        
        # if not self.CT and not self.SPECT and self.arc_correction and self.use_raw_data:
        #     print("Arc correction is not supported for raw data. Disabling arc correction.")
        #     self.arc_correction = False
        
        # if not self.CT and not self.SPECT and self.arc_correction and self.implementation == 1:
        #     print("Arc correction is not supported with implementation 1. Disabling arc correction.")
        #     self.arc_correction = False
        
        if self.Nt < 1:
            print("Number of time steps is less than one. Using one time step.")
            self.partitions = 1
            self.Nt = 1
        
        if self.start > self.end:
            raise ValueError("Start time is later than end time!")
        
        if self.start > self.tot_time:
            raise ValueError("Start time is larger than the total time of the measurement!")
        
        if self.Niter < 1:
            print("Number of iterations is less than one! Using one iteration.")
            self.Niter = 1
        if self.subsets < 1:
            print("Number of subsets is less than one! Using one subset.")
            self.subsets = 1
        
        # if self.implementation == 1 and self.useMultiResolutionVolumes:
        #     raise ValueError("Multi-resolution reconstruction is not supported with implementation 1!")
        
        if self.useMultiResolutionVolumes and not self.useEFOV:
            print("Multi-resolution reconstruction selected, but extended FOV is not selected! Disabling multi-resolution volumes.")
            self.useMultiResolutionVolumes = False
        
        if (self.LSQR or self.CGLS) and self.useMultiResolutionVolumes:
            raise ValueError("Multi-resolution reconstruction is not supported with LSQR or CGLS!")
            
        # if not self.CT and not self.SPECT and self.det_per_ring == self.det_w_pseudo and self.fill_sinogram_gaps:
        #     raise ValueError('Gap filling is only supported with pseudo detectors!')
        if not self.largeDim and not self.x0.any():
            if not self.CT:
                print('Initial value is an empty array, using the default values (1)')
                self.x0 = np.ones((self.Nx, self.Ny, self.Nz), dtype=np.float32, order='F')
            else:
                print('Initial value is an empty array, using the default values (1e-4)')
                self.x0 = np.ones((self.Nx, self.Ny, self.Nz), dtype=np.float32, order='F') * 1e-4
        if not self.largeDim and self.x0.size < self.Nx * self.Ny * self.Nz:
            raise ValueError(f"Initial value has a matrix size smaller ({np.prod(self.x0.shape)}) than the actual image size ({self.Nx*self.Ny*self.Nz})!")
        if not self.largeDim and self.x0.size > self.Nx * self.Ny * self.Nz:
            print(f"Initial value has a matrix size larger ({np.prod(self.x0.shape)}) than the actual image size ({self.Nx*self.Ny*self.Nz})! Attempting automatic resize.")
            try:
                from skimage.transform import resize #scikit-image
                self.x0 = resize(self.x0, (self.Nx, self.Ny, self.Nz))
            except ModuleNotFoundError:
                print('skimage package not found! Unable to perform automatic resize! Install scikit-image package with "pip install scikit-image".')
            # apuVar = np.prod(self.x0.shape) / (self.Nx * self.Ny * self.Nz) ** (1 / 3)
        # if os.name == 'posix':
        #     if len(self.fpath) > 1 and self.fpath[-1] != '/':
        #         self.fpath += '/'
        # elif os.name == 'nt':
        #     if len(self.fpath) > 1 and self.fpath[-1] not in ('\\', '/'):
        #         self.fpath += '\\'
        # else:
        #     if len(self.fpath) > 1 and self.fpath[-1] != '/':
        #         self.fpath += '/'
        # if not self.CT and not self.SPECT and self.use_LMF and self.randoms_correction:
        #     print('Randoms correction is set to true although LMF input is selected. No randoms correction will be performed.')
        #     self.randoms_correction = False
        if self.TV_use_anatomical and self.TV and not os.path.exists(self.TV_referenceImage) and self.MAP and not isinstance(self.TV_referenceImage, np.ndarray):
            raise FileNotFoundError('Anatomical reference image for TV was not found on path!')
        if self.NLM_use_anatomical and self.NLM and not os.path.exists(self.NLM_referenceImage) and self.MAP and not isinstance(self.NLM_referenceImage, np.ndarray):
            raise FileNotFoundError('Anatomical reference image for NLM was not found on path!')
        if self.RDP_use_anatomical and self.RDP and self.RDPIncludeCorners and not os.path.exists(self.RDP_referenceImage) and self.MAP and not isinstance(self.RDP_referenceImage, np.ndarray):
            raise FileNotFoundError('Reference image for RDP was not found on path!')
        if self.precondTypeImage[2] and not os.path.exists(self.referenceImage) and not isinstance(self.referenceImage, np.ndarray):
            raise FileNotFoundError('Reference image for precondititiong was not found on path!')
        if self.RDP_use_anatomical and self.RDP and not self.RDPIncludeCorners:
            raise ValueError('Reference image for RDP is only supported with options.RDPIncludeCorners = True')
        if self.implementation == 2 and self.useCPU and self.RDP and self.RDPIncludeCorners:
            raise ValueError('RDP with include corners is supported only on OpenCL and CUDA!')
        if self.TV and self.TVtype == 2 and not self.TV_use_anatomical:
            raise ValueError('Using TV type = 2, but no anatomical reference set. Use options.TVtype = 1 if anatomical weighting is not used!')
        if self.projector_type not in [1, 2, 3, 4, 5, 6, 11, 14, 12, 13, 21, 22, 31, 32, 33, 41, 51, 15, 44, 45, 54, 55]:
            raise ValueError('The selected projector type is not supported!')
        if self.APLS and not os.path.exists(self.APLS_ref_image) and self.MAP:
            raise FileNotFoundError('APLS selected, but the anatomical reference image was not found on path!')
        if self.epps <= 0:
            print('Epsilon value is zero or less than zero; must be a positive value. Using the default value (1e-6).')
            self.epps = 1e-6
        if self.projector_type in [5, 4, 44, 14, 41, 15, 45, 54, 51, 55] and self.useCPU:
            raise ValueError('The selected projector type is not supported with CPU implementation!')
        if np.sum(self.precondTypeImage) == 0 and (self.PKMA or self.MRAMLA or self.MBSREM):
            print('No image-based preconditioner selected with PKMA/MRAMLA/MBSREM. EM preconditioner is highly recommended!')
        # if len(self.epps) > 1:
        #     print('Epsilon has to be a scalar value! Using the default value (1e-6).')
        #     self.epps = 1e-6
        # if not self.CT and not self.SPECT and self.store_scatter and sum(self.scatter_components) <= 0:
        #     raise ValueError('Store scatter selected, but no scatter components have been selected!')
        # if self.implementation == 1 and not self.precompute_lor and (self.n_rays_transaxial > 1 or self.n_rays_axial > 1):
        #     raise ValueError('Multiray Siddon is not supported with implementation 1!')
        # if self.hyperbolic and self.implementation != 2:
        #     raise ValueError('Hyperbolic prior is only available when using implementation 2!')
        # elif self.implementation == 2 and self.TV and self.TVtype == 3:
        #     self.TV = False
        #     self.hyperbolic = True
        #     if self.TV_use_anatomical:
        #         print('Hyperbolic prior does not support anatomic weighting at the moment when using implementation 2')
        #         self.TV_use_anatomical = False
        # if not self.CT and not self.SPECT and self.use_machine == 2 and self.use_raw_data:
        #     print('Sinogram data cannot be used when raw data is set to true, using raw data instead')
        #     self.use_machine = 1
        if not self.CT and not self.SPECT and self.reconstruct_trues and self.reconstruct_scatter:
            print('Both reconstruct trues and scatter selected, reconstructing only trues.')
            self.reconstruct_scatter = False
        if (self.CGLS or self.LSQR) and self.subsets > 1:
            print('CGLS or LSQR do not support subsets! Setting subsets to 1.')
            self.subsets = 1
        if self.subsets <= 0:
            print('Subsets set to 0 or less than 0. Setting subsets to 1.')
            self.subsets = 1
        if self.TOF_bins > 1 and self.TOF_bins_used == 1 and not self.CT and not self.SPECT:
            print('Summing TOF bins.')
            # self.TOF_bins = self.TOF_bins_used
        
        # if self.implementation == 1 and not self.CT:
        #     if self.projector_type in [2, 3, 31, 32, 21, 12]:
        #         raise ValueError('Orthogonal distance-based/volume-based projector is NOT supported when using implementation 1!')
        #     if self.TOF_bins > 1:
        #         raise ValueError('TOF is not supported with implementation 1!')
        
        # if self.projector_type == 2 and self.implementation == 1 and self.precompute_lor:
        #     warnings.warn('Orthogonal distance-based projector is not recommended with implementation 1!')
        if np.sum(self.precondTypeImage) == 0 and (self.PKMA or self.MRAMLA or self.MBSREM):
            print('No image-based preconditioner selected with PKMA/MRAMLA/MBSREM. EM preconditioner is highly recommended!')
        if self.useCPU and (self.projector_type == 5 or self.projector_type == 45 or self.projector_type == 54 or self.projector_type == 51 or self.projector_type == 15):
            raise ValueError('Selected projector type is not supported with CPU implementation!')
        if self.projector_type == 2 and self.CT:
            raise ValueError('Orthogonal distance-based projector is NOT supported when using CT data!')
        
        if self.projector_type in [5, 15, 51, 45, 54] and not self.CT:
            raise ValueError('Projector type 5 is only supported with CT data!')
        
        if self.projector_type == 6 and not self.SPECT:
            raise ValueError('Projector type 6 is only supported with SPECT data!')
            
        if (not(self.projector_type == 6) and not(self.projector_type == 1) and not(self.projector_type == 11)) and self.SPECT:
            raise ValueError('SPECT only supports projector types 1 and 6!')
        
        if self.projector_type == 6:
            if self.Nx != self.nRowsD:
                raise ValueError('options.Nx has to be same as options.nRowsD when using projector type 6')
            if self.Ny != self.nRowsD:
                raise ValueError('options.Ny has to be same as options.nRowsD when using projector type 6')
            if self.Nz != self.nColsD:
                raise ValueError('options.Nz has to be same as options.nColsD when using projector type 6')
            if self.subsets > 1 and self.subsetType < 8:
                raise ValueError('Subset types 0-7 are not supported with projector type 6!')
        
        if self.FDK and (self.Niter > 1 or self.subsets > 1):
            if self.largeDim:
                self.Niter = 1
            else:
                print('When using FDK/FBP, the number of iterations and subsets must be set as 1. Setting both to 1.')
                self.subsets = 1
                self.Niter = 1
        if self.useCUDA and self.useCPU:
            raise ValueError('Both CUDA and CPU selected! Select only one!')
        
        if self.TOF_bins_used > 1 and (self.projector_type not in [1, 11, 3, 33, 31, 13, 4, 41, 14]) and not self.CT and not self.SPECT:
            raise ValueError('TOF is currently only supported with improved Siddon (projector_type = 1), interpolation-based projector (projector_type = 4) and volume of intersection (projector_type = 3)')
        
        if self.TOF_bins_used > 1 and self.TOF_width <= 0 and not self.CT and not self.SPECT:
            raise ValueError('TOF width (self.TOF_width) must be greater than zero.')
        
        if self.TOF_bins_used > 1 and self.TOF_FWHM == 0 and not self.CT and not self.SPECT:
            raise ValueError('TOF enabled, but the TOF FWHM (self.TOF_FWHM) is zero. FWHM must be nonzero.')
        
        # if self.TOF_bins > 1 and self.use_raw_data and not self.CT and not self.SPECT:
        #     raise ValueError('TOF data is only available with sinogram data. Disable raw data (self.use_raw_data = False).')
        
        if self.corrections_during_reconstruction and (self.scatter_correction or self.randoms_correction) and (self.PDHG or self.PDHGL1 or self.FISTA or self.LSQR or self.CGLS or self.FISTAL1):
            raise ValueError('Randoms/scatter correction cannot be applied during the reconstruction with the selected algorithm!')
            
        
        if self.precondTypeMeas[1] and (self.subsetType < 8 and not(self.subsetType == 4)):
            raise ValueError('Filtering-based preconditioner only works for subset types 4 and 8-11!')
            
        if self.verbose > 0:
            if self.use_ASCII and self.use_machine == 0:
                dispi = 'Using ASCII data'
            # elif self.use_LMF and self.use_machine == 0:
            #     dispi = 'Using LMF data'
            elif self.use_root and self.use_machine == 0:
                dispi = 'Using ROOT data'
            # elif self.use_binary and self.use_machine == 0:
            #     dispi = 'Using BINARY data'
            elif self.use_machine == 1:
                dispi = 'Using data obtained from list-mode file'
            elif self.use_machine == 2:
                dispi = 'Using scanner created sinogram data'
            elif self.use_machine == 3:
                dispi = 'Using 32-bit list-mode data'
            else:
                dispi = None
            
            if self.TOF_bins_used > 1 and self.TOF_FWHM > 0:
                dispi = f"{dispi} with TOF ({self.TOF_bins_used} bins)." if dispi else f"With TOF ({self.TOF_bins_used} bins)."
            else:
                dispi = f"{dispi}." if dispi else "."
            
            if dispi != ".":
                print(dispi)
            
            # if self.only_sinos:
            #     print('Loading only data.')
            if self.useMultiResolutionVolumes:
                if self.transaxialEFOV:
                    Nx = self.NxOrig + round((self.Nx - self.NxOrig) * self.multiResolutionScale)
                    Ny = self.NyOrig + round((self.Ny - self.NyOrig) * self.multiResolutionScale)
                else:
                    Nx = self.Nx
                    Ny = self.Ny
                    
                if self.axialEFOV:
                    Nz = self.NzOrig + round((self.Nz - self.NzOrig) * self.multiResolutionScale)
                else:
                    Nz = self.Nz
            else:
                Nx = self.Nx
                Ny = self.Ny
                Nz = self.Nz
            
            if not (self.compute_normalization or self.only_sinos):
                if self.deviceNum < 0:
                    raise ValueError('Device number has to be positive!')
                try:
                    if not self.useCUDA and af.get_active_backend() != 'opencl':
                        af.set_backend('opencl')
                    dispaus = f"Using implementation {self.implementation} with "
                    info = af.device.info_str()
                    loc = info.find('[' + str(self.deviceNum) + ']')
                    if loc == -1:
                        loc = info.find('-' + str(self.deviceNum) + '-')
                    loc2 = info[loc:].find('(Compute')
                    dispaus += info[loc + 4 : loc + loc2 - 1]
                    print(dispaus)
                except NameError:
                    print('Selected device number is ' + str(self.deviceNum))
                    
                
                if self.projector_type == 1 and not self.precompute_lor:
                    if self.implementation == 1:
                        print("Improved Siddon's algorithm selected with 1 ray.")
                    else:
                        ray = 'rays' if self.n_rays_transaxial > 1 else 'ray'
                        aray = 'rays' if self.n_rays_axial > 1 else 'ray'
                        print(f"Improved Siddon's algorithm selected with {self.n_rays_transaxial} transaxial {ray} and {self.n_rays_axial} axial {aray}.")
                elif self.projector_type == 1 or self.projector_type == 11:
                    print("Improved Siddon's algorithm selected with 1 ray.")
                elif self.projector_type == 2 or self.projector_type == 22:
                    dispi = 'Orthogonal distance-based ray tracer selected'
                    if self.tube_width_z > 0:
                        dispi = f"{dispi} in 3D mode." 
                    else: 
                        dispi = f"{dispi} in 2.5D mode."
                    print(dispi)
                elif self.projector_type == 21:
                    print("Improved Siddon's algorithm selected for forward projection, orthogonal for backprojection.")
                elif self.projector_type == 3 or self.projector_type == 33:
                    print('Volume of intersection based ray tracer selected.');
                elif self.projector_type == 31:
                    print('Improved Siddon''s algorithm selected for forward projection, Volume of intersection based ray tracer for backprojection.')
                elif self.projector_type == 13:
                    print('Volume of intersection based ray tracer selected for forward projection, improved Siddon''s algorithm for backprojection.')
                elif self.projector_type == 4:
                    print('Interpolation-based projector selected.')
                elif self.projector_type == 5:
                    print('Branchless distance-driven based projector selected.')
                elif self.projector_type == 41:
                    print('Interpolation-based projector selected for forward projection, improved Siddon for backprojection.')
                elif self.projector_type == 14:
                    print('Improved Siddon projector selected for forward projection, interpolation-based projector for backprojection.')
                elif self.projector_type == 15:
                    print('Improved Siddon projector selected for forward projection, branchless distance-driven projector for backprojection.')
                elif self.projector_type == 45:
                    print('Interpolation-based projector selected for forward projection, branchless distance-driven projector for backprojection.')
                elif self.projector_type == 54:
                    print('Branchless distance-driven projector selected for forward projection, interpolation-based projector for backprojection.')
                elif self.projector_type == 51:
                    print('Branchless distance-driven projector selected for forward projection, improved Siddon for backprojection.')
                elif self.projector_type == 6:
                    print('Rotation-based projector selected (SPECT).')
        
                if self.use_psf:
                    if self.deblurring:
                        print('PSF ON with deblurring phase.')
                    else:
                        print('PSF ON.')
                if self.attenuation_correction and not self.CT:
                    print('Attenuation correction ON.')
                
                if self.randoms_correction and not self.CT:
                    dispi = 'Randoms correction ON'
                    if self.variance_reduction:
                        dispi += ' with variance reduction'
                        if self.randoms_smoothing:
                            dispi += ' and smoothing'
                    elif self.randoms_smoothing:
                        dispi += ' with smoothing'
                    dispi += '.'
                    print(dispi)
                
                if self.scatter_correction:
                    dispi = 'Scatter correction ON'
                    if self.scatter_smoothing:
                        dispi += ' with smoothing'
                    dispi += '.'
                    print(dispi)
                
                if any(self.precondTypeImage):
                    if self.precondTypeImage[0]:
                        print('Using image-based preconditioning with diagonal preconditioner.')
                    elif self.precondTypeImage[1]:
                        print('Using image-based preconditioning with EM preconditioner.')
                    elif self.precondTypeImage[2]:
                        print('Using image-based preconditioning with IEM preconditioner.')
                
                    if self.precondTypeImage[3]:
                        print('Using image-based preconditioning with momentum preconditioner.')
                
                    if self.precondTypeImage[4]:
                        print('Using image-based preconditioning with normalized gradient preconditioner.')
                
                    if self.precondTypeImage[5]:
                        dispP = 'Using image-based preconditioning with filtering preconditioner'
                        if self.filterWindow == 'hamming':
                            dispP += ' with Hamming window.'
                        elif self.filterWindow == 'hann':
                            dispP += ' with Hann window.'
                        elif self.filterWindow == 'blackman':
                            dispP += ' with Blackman window.'
                        elif self.filterWindow == 'nuttal':
                            dispP += ' with Nuttal window.'
                        elif self.filterWindow == 'gaussian':
                            dispP += ' with Gaussian window.'
                        elif self.filterWindow == 'shepp-logan':
                            dispP += ' with Shepp-Logan window.'
                        elif self.filterWindow == 'cosine':
                            dispP += ' with cosine window.'
                        elif self.filterWindow == 'parzen':
                            dispP += ' with Parzen window.'
                        else:
                            dispP += ' (no windowing).'
                        print(dispP)
                if any(self.precondTypeMeas):
                    if self.precondTypeMeas[0]:
                        print('Using measurement-based preconditioning with diagonal preconditioner.')
                
                    if self.precondTypeMeas[1]:
                        dispP = 'Using measurement-based preconditioning with filtering preconditioner'
                        if self.filterWindow == 'hamming':
                            dispP += ' with Hamming window.'
                        elif self.filterWindow == 'hann':
                            dispP += ' with Hann window.'
                        elif self.filterWindow == 'blackman':
                            dispP += ' with Blackman window.'
                        elif self.filterWindow == 'nuttal':
                            dispP += ' with Nuttal window.'
                        elif self.filterWindow == 'gaussian':
                            dispP += ' with Gaussian window.'
                        elif self.filterWindow == 'shepp-logan':
                            dispP += ' with Shepp-Logan window.'
                        elif self.filterWindow == 'cosine':
                            dispP += ' with cosine window.'
                        elif self.filterWindow == 'parzen':
                            dispP += ' with Parzen window.'
                        else:
                            dispP += ' (no windowing).'
                        print(dispP)
                
                if self.oOffsetZ != 0 or self.oOffsetX != 0 or self.oOffsetY != 0:
                    print(f'Object offset is [{self.oOffsetX}, {self.oOffsetY}, {self.oOffsetZ}] (XYZ).')
                    
            if self.normalization_correction and not self.CT:
                print('Normalization correction ON.')
            elif self.normalization_correction and self.compute_normalization and not self.CT:
                print('Normalization correction cannot be applied when computing normalization coefficients. Disabling normalization correction.')
                self.normalization_correction = False
            elif self.compute_normalization and not self.CT:
                print('Computing normalization coefficients.')
            if not (self.compute_normalization or self.only_sinos):
                if self.corrections_during_reconstruction and (self.normalization_correction or self.randoms_correction or self.scatter_correction):
                    print('Corrections applied during reconstruction (ordinary Poisson).')
                elif not self.corrections_during_reconstruction and (self.normalization_correction or self.randoms_correction or self.scatter_correction):
                    print('Corrections applied to the measurement data.')
            
                if self.arc_correction and not self.use_raw_data and not self.CT:
                    print('Arc correction ON.')
                if self.Nt == 1:
                    if self.CT:
                        dispi = 'Using STATIC projection data'
                    else:
                        dispi = 'Using STATIC sinogram data'
                else:
                    if self.CT:
                        dispi = 'Using DYNAMIC projection data'
                    else:
                        dispi = 'Using DYNAMIC sinogram data'
                
                if self.reconstruct_trues:
                    dispi += ' (trues)'
                elif self.reconstruct_scatter:
                    dispi += ' (scatter).'
                else:
                    dispi += ' (prompts)'
                
                if self.sampling > 1:
                    dispi += f' with {self.sampling}x sampling'
                
                if self.Nt > 1:
                    if self.sampling > 1:
                        dispi += f' and {self.Nt} time steps'
                    else:
                        dispi += f' with {self.Nt} time steps'
                
                dispi += '.'
                print(dispi)
                if self.subsets > 1:
                    if self.subsetType == 1:
                        print(f'Every {self.subsets}th column measurement is taken per subset.')
                    elif self.subsetType == 2:
                        print(f'Every {self.subsets}th row measurement is taken per subset.')
                    elif self.subsetType == 3:
                        print('Using random subset sampling.')
                    elif self.subsetType == 4:
                        print(f'Every {self.subsets}th sinogram column is taken per subset.')
                    elif self.subsetType == 5:
                        print(f'Every {self.subsets}th sinogram row is taken per subset.')
                    elif self.subsetType == 6:
                        print(f'Using angle-based subset sampling with {self.n_angles} angles combined per subset.')
                    elif self.subsetType == 7:
                        print('Using golden angle-based subset sampling.')
                    elif self.subsetType == 8:
                        print(f'Using every {self.subsets}th sinogram/projection image.')
                    elif self.subsetType == 9:
                        print('Using sinograms/projection images in random order.')
                    elif self.subsetType == 10:
                        print('Using golde angle sampling with sinograms/projections.')
                    elif self.subsetType == 11:
                        print('Using prime factor ordering of projections/sinograms into subsets.')
                    elif self.subsetType == 0:
                        print(f'Dividing data into {self.subsets} segments.')
                
                print(f'Using an image (matrix) size of {Nx}x{Ny}x{Nz} with {self.Niter} iterations and {self.subsets} subsets.')
            elif self.CT:
                print(f'Using an image (matrix) size of {Nx}x{Ny}x{Nz} with {self.Niter} iterations and {self.subsets} subsets.')

    
    
    def setUpCorrections(self):
        nx = self.Nx
        ny = self.Ny
        nz = self.Nz
        self.NxFull = nx
        self.NyFull = ny
        self.NzFull = nz
        if self.useEFOV:
            if self.useMultiResolutionVolumes:
                if self.axialEFOV == True and self.transaxialEFOV == True:
                    self.nMultiVolumes = 6;
                    if np.ceil(self.Nx * self.multiResolutionScale) == np.floor(self.NxOrig * self.multiResolutionScale) + np.ceil((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale)*2:
                        self.Nx = np.row_stack((self.NxOrig, np.floor(self.NxOrig * self.multiResolutionScale), np.floor(self.NxOrig * self.multiResolutionScale), 
                            np.ceil((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale), np.ceil((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale), 
                            np.floor(self.NxOrig * self.multiResolutionScale), np.floor(self.NxOrig * self.multiResolutionScale))).astype(np.uint32)
                    elif np.ceil(self.Nx * self.multiResolutionScale) == np.ceil(self.NxOrig * self.multiResolutionScale) + np.floor((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale)*2:
                        self.Nx = np.row_stack((self.NxOrig, np.ceil(self.NxOrig * self.multiResolutionScale), np.ceil(self.NxOrig * self.multiResolutionScale), 
                            np.floor((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale), np.floor((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale), 
                            np.ceil(self.NxOrig * self.multiResolutionScale), np.ceil(self.NxOrig * self.multiResolutionScale))).astype(np.uint32)
                    else:
                        self.Nx = np.row_stack((self.NxOrig, np.ceil(self.NxOrig * self.multiResolutionScale), np.ceil(self.NxOrig * self.multiResolutionScale), 
                            np.ceil((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale), np.ceil((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale), 
                            np.ceil(self.NxOrig * self.multiResolutionScale), np.ceil(self.NxOrig * self.multiResolutionScale))).astype(np.uint32)
                    if np.ceil(self.Ny * self.multiResolutionScale) == np.floor(self.NyOrig * self.multiResolutionScale) + np.ceil((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale)*2:
                        self.Ny = np.row_stack((self.NyOrig, np.floor(self.NyOrig * self.multiResolutionScale), np.floor(self.NyOrig * self.multiResolutionScale), 
                            np.ceil(self.Ny * self.multiResolutionScale), np.ceil(self.Ny * self.multiResolutionScale),
                            np.ceil((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale), np.ceil((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale))).astype(np.uint32)
                    elif np.ceil(self.Ny * self.multiResolutionScale) == np.ceil(self.NyOrig * self.multiResolutionScale) + np.floor((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale)*2:
                        self.Ny = np.row_stack((self.NyOrig, np.ceil(self.NyOrig * self.multiResolutionScale), np.ceil(self.NyOrig * self.multiResolutionScale), 
                            np.ceil(self.Ny * self.multiResolutionScale), np.ceil(self.Ny * self.multiResolutionScale),
                            np.floor((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale), np.floor((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale))).astype(np.uint32)
                    else:
                        self.Ny = np.row_stack((self.NyOrig, np.ceil(self.NyOrig * self.multiResolutionScale), np.ceil(self.NyOrig * self.multiResolutionScale), 
                            np.ceil(self.Ny * self.multiResolutionScale), np.ceil(self.Ny * self.multiResolutionScale),
                            np.ceil((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale), np.ceil((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale))).astype(np.uint32)
                    self.Nz = np.row_stack((self.NzOrig, np.ceil((self.Nz - self.NzOrig) / 2 * self.multiResolutionScale), np.ceil((self.Nz - self.NzOrig) / 2 * self.multiResolutionScale), 
                        np.ceil(self.Nz * self.multiResolutionScale), np.ceil(self.Nz * self.multiResolutionScale), 
                        np.ceil(self.Nz * self.multiResolutionScale), np.ceil(self.Nz * self.multiResolutionScale))).astype(np.uint32)
                    self.FOVa_x = np.row_stack((self.FOVxOrig, self.FOVxOrig, self.FOVxOrig, 
                        (self.FOVa_x - self.FOVxOrig) / 2, (self.FOVa_x - self.FOVxOrig) / 2, 
                        self.FOVxOrig, self.FOVxOrig)).astype(np.float32)
                    self.FOVa_y = np.row_stack((self.FOVyOrig, self.FOVyOrig, self.FOVyOrig, 
                        self.FOVa_y, self.FOVa_y, 
                        (self.FOVa_y - self.FOVyOrig) / 2, (self.FOVa_y - self.FOVyOrig) / 2)).astype(np.float32)
                    self.axial_fov = np.row_stack((self.axialFOVOrig, (self.axial_fov - self.axialFOVOrig) /2, (self.axial_fov - self.axialFOVOrig) /2, 
                        self.axial_fov, self.axial_fov, self.axial_fov, self.axial_fov)).astype(np.float32)
                    if self.x0.shape[0] == nx:
                        try:
                            from skimage.transform import downscale_local_mean #scikit-image
                        except ModuleNotFoundError:
                            print('skimage package not found! Unable to perform automatic resize! Install scikit-image package with "pip install scikit-image".')
                        apu = downscale_local_mean(self.x0, int(1 / self.multiResolutionScale))
                        # apu = zoom(self.x0, self.multiResolutionScale)
                        x1 = apu[self.Nx[3].item() : self.Nx[3].item() + self.Nx[1].item(), self.Ny[5].item() : self.Ny[5].item() + self.Ny[1].item(), 
                            0 : self.Nz[1].item()]
                        x2 = apu[self.Nx[4].item() : self.Nx[4].item() + self.Nx[2].item(), self.Ny[6].item() : self.Ny[6].item() + self.Ny[2].item(), 
                            -self.Nz[1].item():]
                        x3 = apu[0 : self.Nx[3].item(), :, :]
                        # if apu.shape[0] % 2 == 0:
                        #     x4 = apu[self.Nx[4].item() + self.Nx[2].item() - 1: , :, :]
                        # else:
                        x4 = apu[self.Nx[4].item() + self.Nx[2].item() : , :, :]
                        x5 = apu[self.Nx[3].item() : self.Nx[3].item() + self.Nx[1].item(), 0:self.Ny[5].item(), 
                            0 : self.Nz[3].item()]
                        # if apu.shape[1] % 2 == 0:
                        x6 = apu[self.Nx[4].item() : self.Nx[4].item() + self.Nx[2].item(), self.Ny[6].item() + self.Ny[2].item() : , 
                                 0 : self.Nz[4].item()]
                        # else:
                        #     x6 = apu[self.Nx[4].item() : self.Nx[4].item() + self.Nx[2].item(), self.Ny[6].item() + self.Ny[2].item() - 1 : , 
                        #              0 : self.Nz[4].item()]
                        self.x0 = self.x0[(self.x0.shape[0] - self.NxOrig) // 2 : (self.x0.shape[0] - self.NxOrig) // 2 + self.NxOrig, 
                            (self.x0.shape[1] - self.NyOrig) // 2 : (self.x0.shape[1] - self.NyOrig) // 2 + self.NyOrig,
                            (self.x0.shape[2] - self.NzOrig) // 2 : (self.x0.shape[2] - self.NzOrig) // 2 + self.NzOrig]
                        # if self.implementation == 2 or self.implementation == 3 or self.implementation == 5:
                        self.x0 = np.concatenate([self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F'), x3.ravel('F'), x4.ravel('F'), x5.ravel('F'), x6.ravel('F')])
                elif self.transaxialEFOV == True and self.axialEFOV == False:
                    self.nMultiVolumes = 4
                    self.Nx = np.row_stack((self.NxOrig, np.ceil((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale), np.ceil((self.Nx - self.NxOrig) / 2 * self.multiResolutionScale), 
                        np.ceil(self.NxOrig * self.multiResolutionScale), np.ceil(self.NxOrig * self.multiResolutionScale))).astype(np.uint32)
                    self.Ny = np.row_stack((self.NyOrig, np.ceil(self.Ny * self.multiResolutionScale), np.ceil(self.Ny * self.multiResolutionScale),
                        np.ceil((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale), np.ceil((self.Ny - self.NyOrig) / 2 * self.multiResolutionScale))).astype(np.uint32)
                    self.Nz = np.row_stack((self.Nz, np.ceil(self.Nz * self.multiResolutionScale), np.ceil(self.Nz * self.multiResolutionScale), 
                        np.ceil(self.Nz * self.multiResolutionScale), np.ceil(self.Nz * self.multiResolutionScale))).astype(np.uint32)
                    self.FOVa_x = np.row_stack((self.FOVxOrig, (self.FOVa_x - self.FOVxOrig) / 2, (self.FOVa_x - self.FOVxOrig) / 2, 
                        self.FOVxOrig, self.FOVxOrig)).astype(np.float32)
                    self.FOVa_y = np.row_stack((self.FOVyOrig, self.FOVa_y, self.FOVa_y, 
                        (self.FOVa_y - self.FOVyOrig) / 2, (self.FOVa_y - self.FOVyOrig) / 2)).astype(np.float32)
                    self.axial_fov = np.row_stack((self.axial_fov, self.axial_fov, self.axial_fov,
                        self.axial_fov, self.axial_fov)).astype(np.float32)
                    if self.x0.shape[0] == nx:
                        try:
                            from skimage.transform import downscale_local_mean #scikit-image
                        except ModuleNotFoundError:
                            print('skimage package not found! Unable to perform automatic resize! Install scikit-image package with "pip install scikit-image".')
                        apu = downscale_local_mean(self.x0, int(1 / self.multiResolutionScale))
                        # apu = zoom(self.x0, self.multiResolutionScale)
                        x1 = apu[0 : self.Nx[1].item(), :, :]
                        # if apu.shape[0] % 2 == 0:
                        x2 = apu[self.Nx[1].item() + self.Nx[3].item() : , :, :]
                        # else:
                        #     x2 = apu[self.Nx[1].item() + self.Nx[3].item() - 1 : , :, :]
                        x3 = apu[self.Nx[1].item() : self.Nx[1].item() + self.Nx[3].item(), 0:self.Ny[3].item(), :]
                        x4 = apu[self.Nx[1].item() : self.Nx[1].item() + self.Nx[4].item(), self.Ny[1].item() - self.Ny[3].item() : , :]
                        self.x0 = np.single(self.x0[(self.x0.shape[0] - self.NxOrig) // 2 : (self.x0.shape[0] - self.NxOrig) // 2 + self.NxOrig, 
                                                    (self.x0.shape[1] - self.NyOrig) // 2 : (self.x0.shape[1] - self.NyOrig) // 2 + self.NyOrig, :])
                        # if self.implementation == 2 or self.implementation == 3 or self.implementation == 5
                        self.x0 = np.concatenate([self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F'), x3.ravel('F'), x4.ravel('F')])
                else:
                    self.nMultiVolumes = 2;
                    self.Nx = np.row_stack((self.Nx, np.ceil(self.Nx * self.multiResolutionScale), np.ceil(self.Nx * self.multiResolutionScale))).astype(np.uint32)
                    self.Ny = np.row_stack((self.Ny, np.ceil(self.Ny * self.multiResolutionScale), np.ceil(self.Ny * self.multiResolutionScale))).astype(np.uint32)
                    self.Nz = np.row_stack((self.NzOrig, np.ceil((self.Nz - self.NzOrig) // 2 * self.multiResolutionScale), np.ceil((self.Nz - self.NzOrig) / 2 * self.multiResolutionScale))).astype(np.uint32)
                    self.FOVa_x = np.row_stack((self.FOVa_x, self.FOVa_x, self.FOVa_x)).astype(np.float32)
                    self.FOVa_y = np.row_stack((self.FOVa_y, self.FOVa_y, self.FOVa_y)).astype(np.float32)
                    self.axial_fov = np.row_stack((self.axialFOVOrig, (self.axial_fov - self.axialFOVOrig) / 2, (self.axial_fov - self.axialFOVOrig) / 2)).astype(np.float32)
                    if self.x0.shape[0] == nx:
                        try:
                            from skimage.transform import downscale_local_mean #scikit-image
                        except ModuleNotFoundError:
                            print('skimage package not found! Unable to perform automatic resize! Install scikit-image package with "pip install scikit-image".')
                        apu = downscale_local_mean(self.x0, int(1 / self.multiResolutionScale))
                        # apu = zoom(self.x0, self.multiResolutionScale)
                        x1 = apu[:, :, :self.Nz[1].item()]
                        x2 = apu[:, :, -self.Nz[2].item():]
                        self.x0 = np.single(self.x0[:, :, (self.x0.shape[2] - self.NzOrig) // 2 : (self.x0.shape[2] - self.NzOrig) // 2 + self.NzOrig])
                        # if self.implementation == 2 or self.implementation == 3 or self.implementation == 5
                        self.x0 = np.concatenate([self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F')])
                        # self.x0 = np.row_stack((self.x0.flatten('F'),x1.flatten('F'), x2.flatten('F')))
                self.NxPrior = self.Nx[0].item()
                self.NyPrior = self.Ny[0].item()
                self.NzPrior = self.Nz[0].item()
            else:
                if self.eFOVIndices.size < 1:
                    self.eFOVIndices = np.zeros((self.Nz,1), dtype=np.uint8)
                    self.eFOVIndices[(self.Nz - self.NzOrig)//2 : -(self.Nz - self.NzOrig)//2 - 1] = 1
                self.NzPrior = np.sum(self.eFOVIndices, dtype=np.uint32)
                self.maskPrior = np.zeros((self.Nx, self.Ny), dtype=np.uint8)
                self.maskPrior[(self.Nx - self.NxOrig)//2 : -(self.Nx - self.NxOrig)//2 - 1, (self.Ny - self.NyOrig)//2 : -(self.Ny - self.NyOrig)//2 - 1] = 1
                self.NxPrior = np.sum(self.maskPrior[:,int(np.round(self.maskPrior.shape[1]/2))], dtype=np.uint32)
                self.NyPrior = np.sum(self.maskPrior[int(np.round(self.maskPrior.shape[0]/2)),:], dtype=np.uint32)
                if self.useMaskBP:
                    self.maskPrior = self.maskPrior + (1 - self.maskBP)
        if self.offsetCorrection:
            self.OffsetLimit = np.zeros((self.nProjections, 1), dtype=np.float32);
            for kk in range(self.nProjections):
                sx = self.x[kk, 0];
                sy = self.x[kk, 1];
                dx = self.x[kk, 3];
                dy = self.x[kk, 4];
                ii = np.arange(0., self.nRowsD, 0.25, dtype=np.float32) - self.nRowsD / 2
                dx = dx + self.z[kk, 0] * ii + self.z[kk, 3] * ii
                dy = dy + self.z[kk, 1] * ii + self.z[kk, 4] * ii
                dist = np.abs((dx - sx) * (sy) - ((sx) * (dy - sy))) / np.sqrt((dx - sx)**2 + (dy - sy)**2)
                ind = np.argmin(dist)
                self.OffsetLimit[kk] = ind * self.dPitchY/4
            
            
    def PSFKernel(self):
        def gaussianKernel(x, y, z, sigma_x, sigma_y, sigma_z = 0.):
            if sigma_z == 0.:
                gaussK = np.exp(-(x**2 / (2. * sigma_x**2) + y**2 / (2. * sigma_y**2)))
            else:
                gaussK = np.exp(-(x**2 / (2. * sigma_x**2) + y**2 / (2. * sigma_y**2) + z**2 / (2. * sigma_z**2)))
            gaussK = gaussK / np.sum(gaussK)
            gaussK = gaussK.astype(dtype=np.float32)
            return gaussK
        if self.use_psf:
            g_pituus_z = 0
            g_pituus_x = math.ceil(2. * (self.FWHM[0].item() / (2. * math.sqrt(2. * math.log(2.)))) / self.dx[0].item())
            g_pituus_y = math.ceil(2. * (self.FWHM[1].item() / (2. * math.sqrt(2. * math.log(2.)))) / self.dy[0].item())
            if len(self.FWHM) >= 3:
                g_pituus_z = math.ceil(2. * (self.FWHM[2].item() / (2. * math.sqrt(2. * math.log(2.)))) / self.dz[0].item())
            g_x = np.reshape(np.linspace(-g_pituus_x * self.dx[0].item(), g_pituus_x * self.dx[0].item(), 2 * g_pituus_x + 1, dtype=np.float32), (-1, 1, 1), order='F')
            g_y = np.reshape(np.linspace(-g_pituus_y * self.dy[0].item(), g_pituus_y * self.dy[0].item(), 2 * g_pituus_y + 1, dtype=np.float32), (1, -1, 1), order='F')
            if len(self.FWHM) >= 3:
                g_z = np.zeros((1,1,g_pituus_z * 2 + 1), dtype=np.float32, order='F')
                g_z[0,0,:] = np.reshape(np.linspace(-g_pituus_z * self.dz[0], g_pituus_z * self.dz[0], 2*g_pituus_z + 1, dtype=np.float32), (1,1,-1), order='F')
                self.gaussK = gaussianKernel(g_x, g_y, g_z, self.FWHM[0].item() / (2. * math.sqrt(2. * math.log(2.))), self.FWHM[1].item() / (2. * math.sqrt(2. * math.log(2.))), self.FWHM[2].item() / (2. * math.sqrt(2. * math.log(2.))));
            else:
                self.gaussK = gaussianKernel(g_x, g_y, 0, self.FWHM[0].item() / (2. * math.sqrt(2. * math.log(2.))), self.FWHM[1].item() / (2. * math.sqrt(2. * math.log(2.))))
            self.gaussK = np.asfortranarray(self.gaussK)
            self.g_dim_x = g_pituus_x
            self.g_dim_y = g_pituus_y
            self.g_dim_z = g_pituus_z

    def initProj(self):
        from omegatomo.reconstruction.prepass import prepassPhase
        from omegatomo.reconstruction.prepass import parseInputs
        from omegatomo.reconstruction.prepass import loadCorrections
        if self.useAF:
            # import arrayfire as af
            if af.get_active_backend() != 'opencl' and not self.useCUDA:
                af.set_backend('opencl')

            af.device.set_device(self.deviceNum)
        if self.useTorch and self.useAF:
            raise ValueError('Arrayfire and PyTorch cannot be used at the same time! Select only one!')
        if self.useTorch:
            import torch
            torch.cuda.init()
        if self.useTorch and not self.useCUDA:
            raise ValueError('PyTorch does not work with OpenCL! You can still use OpenCL manually with PyTorch, but you have to manually transfer the OpenCL data first to host (NumPy array) and then to Torch (or vice versa, i.e. Torch --> NumPy --> OpenCL)')
        if self.useCuPy and self.useCUDA:
            import cupy as cp
        
        self.NVOXELS = 8
        self.TH = 100000000000.
        self.TH32 = 100000.
        self.NVOXELS5 = 1
        self.NVOXELSFP = 8
        if np.size(self.weights) > 0:
            self.empty_weight = False
        mDataFound = self.SinM.size > 0
        loadCorrections(self)
        parseInputs(self, mDataFound)
        prepassPhase(self)
        if self.listmode > 0 and self.subsets > 1 and self.subsetType > 0:
            if self.useIndexBasedReconstruction:
                self.trIndex = self.trIndex[:,self.index]
                self.axIndex = self.axIndex[:,self.index]
            # else:
            #     self.x = np.reshape(self.x, [6, -1], order='F')
            #     self.x = self.x[:,self.index]
            #     self.x = self.x.ravel('F')
        if self.useIndexBasedReconstruction and self.listmode > 0:
            self.trIndex = self.trIndex.ravel('F')
            self.axIndex = self.axIndex.ravel('F')

        if self.projector_type in [1, 11, 14, 15, 12, 13]:
            self.FPType = 1
        elif self.projector_type in [2, 21, 22, 23, 24, 25]:
            self.FPType = 2
        elif self.projector_type in [3, 31, 32, 33, 34, 35]:
            self.FPType = 3
        elif self.projector_type in [4, 41, 42, 43, 44, 45]:
            self.FPType = 4
        elif self.projector_type in [5, 51, 52, 53, 54, 55]:
            self.FPType = 5
        elif self.projector_type == 6:
            self.FPType = 6
        else:
            raise ValueError('Invalid forward projector!')
        if self.projector_type in [1, 11, 21, 31, 41, 51]:
            self.BPType = 1
        elif self.projector_type in [2, 12, 22, 32, 42, 52]:
            self.BPType = 2
        elif self.projector_type in [3, 13, 23, 33, 43, 53]:
            self.BPType = 3
        elif self.projector_type in [4, 14, 24, 34, 44, 54]:
            self.BPType = 4
        elif self.projector_type in [5, 15, 25, 35, 45, 55]:
            self.BPType = 5
        elif self.projector_type == 6:
            self.BPType = 6
        else:
            raise ValueError('Invalid backprojector!')
        # if self.useAF == False and (self.FPType == 5 or self.BPType == 5):
        #     raise ValueError('Branchless distance-driven (projector type 5) can only be used with Arrayfire!')
        if (self.useAF == False and self.useCuPy == False) and self.projector_type == 6:
            raise ValueError('Projector type 6 can only be used with Arrayfire (OpenCL) or CuPy (CUDA)!')
            
        if not self.projector_type == 6:
            headerDir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'opencl')) + "/"
            with open(headerDir + 'general_opencl_functions.h', encoding="utf8") as f:
                hlines = f.read()
            if self.FPType in [1, 2, 3]:
                with open(headerDir + 'projectorType123.cl', encoding="utf8") as f:
                    linesFP = f.read()
            elif self.FPType in [4]:
                with open(headerDir + 'projectorType4.cl', encoding="utf8") as f:
                    linesFP = f.read()
            elif self.FPType in [5]:
                with open(headerDir + 'projectorType5.cl', encoding="utf8") as f:
                    linesFP = f.read()
            if self.BPType in [1, 2, 3]:
                with open(headerDir + 'projectorType123.cl', encoding="utf8") as f:
                    linesBP = f.read()
            elif self.BPType in [4]:
                with open(headerDir + 'projectorType4.cl', encoding="utf8") as f:
                    linesBP = f.read()
            elif self.BPType in [5]:
                with open(headerDir + 'projectorType5.cl', encoding="utf8") as f:
                    linesBP = f.read()
            globalSize = [None] * self.subsets
            # self.mSize = [None] * self.subsets
            for i in range(self.subsets):
                if (self.FPType == 5):
                    globalSize[i] = (self.nRowsD, (self.nColsD + self.NVOXELSFP - 1) // self.NVOXELSFP, self.nProjSubset[i].item())
                    localSize = (16, 16, 1)
                    erotus = (localSize[0] - (globalSize[i][0] % localSize[0]), localSize[1] - (globalSize[i][1] % localSize[1]), 0)
                    globalSize[i] = (self.nRowsD + erotus[0], (self.nColsD + self.NVOXELSFP - 1) // self.NVOXELSFP + erotus[1], self.nProjSubset[i].item())
                elif ((self.CT or self.SPECT or self.PET) and self.listmode == 0):
                    globalSize[i] = (self.nRowsD, self.nColsD, self.nProjSubset[i].item())
                    localSize = (16, 16, 1)
                    erotus = (localSize[0] - (globalSize[i][0] % localSize[0]), localSize[1] - (globalSize[i][1] % localSize[1]), 0)
                    globalSize[i] = (self.nRowsD + erotus[0], self.nColsD + erotus[1], self.nProjSubset[i].item())
                else:
                    globalSize[i] = (self.nMeasSubset[i].item(), 1, 1)
                    localSize = (128, 1, 1)
                    erotus = (localSize[0] - (globalSize[i][0] % localSize[0]), localSize[1] - (globalSize[i][1] % localSize[1]), 0)
                    globalSize[i] = (self.nMeasSubset[i].item() + erotus[0], 1, 1)
            self.globalSizeFP = globalSize.copy()
            self.localSizeFP = localSize + tuple()
            self.erotusBP = [0] * (self.nMultiVolumes + 1) * 2
            localSize = (16, 16, 1)
            for ii in range(self.nMultiVolumes + 1):
                apu = [self.Nx[ii].item() % localSize[0], self.Ny[ii].item() % localSize[1], 0]
                if apu[0] > 0:
                    self.erotusBP[ii * 2] = localSize[0] - apu[0]
                if apu[1] > 0:
                    self.erotusBP[ii * 2 + 1] = localSize[1] - apu[1]
            
            if self.BPType in [1, 2, 3] or (self.BPType == 4 and not self.CT):
                globalSize = [[None] * (self.nMultiVolumes + 1)] * self.subsets
                for i in range(self.subsets):
                    for ii in range(self.nMultiVolumes + 1):
                        globalSize[i][ii] = self.globalSizeFP[i]
                self.localSizeBP = self.localSizeFP + tuple()
            else:
                globalSize = [[None] * (self.nMultiVolumes + 1)] * self.subsets
                for i in range(self.subsets):
                    for ii in range(self.nMultiVolumes + 1):
                        if self.BPType == 4:
                            globalSize[i][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], (self.Nz[ii].item()  + self.NVOXELS - 1) // self.NVOXELS)
                        elif self.BPType == 5:
                            if self.pitch:
                                globalSize[i][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
                            else:
                                globalSize[i][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], (self.Nz[ii].item()  + self.NVOXELS5 - 1) // self.NVOXELS5)
                        else:
                            globalSize[i][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
                self.localSizeBP = localSize + tuple()
            self.globalSizeBP = globalSize.copy()
                                
            
            self.Nxy = self.Nx[0].item() * self.Ny[0].item()
            if self.useCUDA:
                if self.useImages and not self.useCuPy:
                    self.useImages = False
                if self.use_64bit_atomics or self.use_32bit_atomics:
                    self.use_64bit_atomics = False
                    self.use_32bit_atomics = False
                bOpt = ('-DCUDA','-DPYTHON',)
            else:
                bOpt =('-cl-single-precision-constant -DOPENCL',)
            if self.useMAD:
                if self.useCUDA:
                    bOpt += ('--use_fast_math','-DUSEMAD',)
                else:
                    bOpt += (' -cl-fast-relaxed-math -DUSEMAD',)
            if self.useImages:
                bOpt += ('-DUSEIMAGES',)
            if (self.FPType == 2 or self.BPType == 2 or self.FPType == 3 or self.BPType == 3):
                if self.orthTransaxial:
                    bOpt += ('-DCRYSTXY',)
                if self.orthAxial:
                    bOpt += ('-DCRYSTZ',)
                with open(headerDir + 'opencl_functions_orth3d.h') as f:
                    hlines2 = f.read()
                if self.FPType in [2, 3]:
                    linesFP = hlines + hlines2 + linesFP
                else:
                    linesFP = hlines + linesFP
                if self.BPType in [2, 3]:
                    linesBP = hlines + hlines2 + linesBP
                else:
                    linesBP = hlines + linesBP
            else:
                linesFP = hlines + linesFP
                linesBP = hlines + linesBP
            if self.FPType == 3 or self.BPType == 3:
                bOpt += ('-DVOL',)
            if self.useMaskFP:
                bOpt += ('-DMASKFP',)
            if self.useMaskBP:
                bOpt += ('-DMASKBP',)
            if self.OffsetLimit.size > 0:
                bOpt += ('-DOFFSET',)
            if self.attenuation_correction and self.CTAttenuation:
                bOpt += ('-DATN',)
            elif self.attenuation_correction and not self.CTAttenuation:
                bOpt += ('-DATNM',)
            if self.normalization_correction:
                bOpt += ('-DNORM',)
            if self.additionalCorrection:
                bOpt += ('-DSCATTER',)
            if self.randoms_correction:
                bOpt += ('-DRANDOMS',)
            if self.nLayers > 1:
                if self.useIndexBasedReconstruction:
                    bOpt += ('-DNLAYERS=' + str(self.nLayers),)
                else:
                    bOpt += ('-DNLAYERS=' + str(self.nProjections // (self.nLayers * self.nLayers)),)
            if self.TOF:
                bOpt += ('-DTOF',)
            if self.CT:
                bOpt += ('-DCT',)
            elif self.SPECT:
                bOpt += ('-DSPECT',)
            elif self.PET:
                bOpt += ('-DPET',)

            bOpt += ('-DNBINS=' + str(self.TOF_bins_used),)
            if self.listmode:
                bOpt += ('-DLISTMODE',)
            if self.listmode > 0 and self.useIndexBasedReconstruction:
                bOpt += ('-DINDEXBASED',)
            if (self.FPType == 1 or self.BPType == 1 or self.FPType == 4 or self.BPType == 4) and self.n_rays_transaxial * self.n_rays_axial > 1:
                bOpt += ('-DN_RAYS=' + str(self.n_rays_transaxial * self.n_rays_axial),)
                bOpt += ('-DN_RAYS2D=' + str(self.n_rays_transaxial),)
                bOpt += ('-DN_RAYS3D=' + str(self.n_rays_axial),)
            if self.pitch:
                bOpt += ('-DPITCH',)
            if (((self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7))) and not self.CT and not self.SPECT and not self.PET and self.listmode == 0):
                bOpt += ('-DSUBSETS',)
            if self.subsets > 1 and self.listmode == 0:
                bOpt += ('-DSTYPE=' + str(self.subsetType),'-DNSUBSETS=' + str(self.subsets),)
            
            bOptFP = bOpt + ('-DFP',)
            if self.localSizeFP[1] > 1:
                bOptFP += ('-DLOCAL_SIZE=' + str(self.localSizeFP[0]),'-DLOCAL_SIZE2=' + str(self.localSizeFP[1]),)
            else:
                bOptFP += ('-DLOCAL_SIZE=' + str(self.localSizeFP[0]),'-DLOCAL_SIZE2=' + str(1),)
            if self.FPType in [1, 2, 3]:
                bOptFP += ('-DSIDDON',)
                if self.FPType in [2, 3]:
                    bOptFP += ('-DORTH',)
                if self.FPType == 3:
                    bOptFP += ('-DVOL',)
                if self.use_64bit_atomics:
                    bOptFP += ('-DCAST=long',)
                elif self.use_32bit_atomics:
                    bOptFP += ('-DCAST=int',)
                else:
                    bOptFP += ('-DCAST=float',)
            elif self.FPType == 4:
                bOptFP += ('-DPTYPE4','-DNVOXELS=' + str(self.NVOXELS),)
                if not self.CT:
                    if self.use_64bit_atomics:
                        bOptFP += ('-DCAST=long',)
                    elif self.use_32bit_atomics:
                        bOptFP += ('-DCAST=int',)
                    else:
                        bOptFP += ('-DCAST=float',)
            elif self.FPType == 5:
                bOptFP += ('-DPROJ5','-DNVOXELSFP=' + str(self.NVOXELSFP),)
                if self.meanFP:
                    bOptFP += ('-DMEANDISTANCEFP',)
            
            bOptBP = bOpt + ('-DBP',)
            if self.localSizeBP[1] > 1:
                bOptBP += ('-DLOCAL_SIZE=' + str(self.localSizeBP[0]),'-DLOCAL_SIZE2=' + str(self.localSizeBP[1]),)
            else:
                bOptBP += ('-DLOCAL_SIZE=' + str(self.localSizeBP[0]),'-DLOCAL_SIZE2=' + str(1),)
            if self.BPType in [1, 2, 3]:
                bOptBP += ('-DSIDDON',)
                if self.BPType in [2, 3]:
                    bOptBP += ('-DORTH',)
                if self.BPType == 3:
                    bOptBP += ('-DVOL',)
                bOptBP += ('-DATOMICF',)
                if self.use_64bit_atomics:
                    bOptBP += ('-DATOMIC','-DCAST=long','-DTH=' + str(self.TH),)
                elif self.use_32bit_atomics:
                    bOptBP += (' -DATOMIC32',' -DCAST=int',' -DTH=' + str(self.TH32),)
                else:
                    bOptBP += ('-DCAST=float',)
            elif self.BPType == 4:
                bOptBP += ('-DPTYPE4','-DNVOXELS=' + str(self.NVOXELS),)
                if not self.CT:
                    bOptBP += ('-DATOMICF',)
                    if self.use_64bit_atomics:
                        bOptBP += ('-DATOMIC','-DCAST=long','-DTH=' + str(self.TH),)
                    elif self.use_32bit_atomics:
                        bOptBP += (' -DATOMIC32',' -DCAST=int',' -DTH=' + str(self.TH32),)
                    else:
                        bOptBP += ('-DCAST=float',)
            elif self.BPType == 5:
                bOptBP += ('-DPROJ5','-DNVOXELS5=' + str(self.NVOXELS5),)
                if self.meanBP:
                    bOptBP += ('-DMEANDISTANCEBP',)
        else:
            if self.useCUDA:
                if self.useTorch:
                    self.gFilter = np.ascontiguousarray(self.gFilter)
                    self.gFilter = np.transpose(self.gFilter, (1, 0, 2, 3))
                    self.d_gFilter = torch.tensor(self.gFilter, device='cuda')
                    self.angles = np.degrees(self.angles)
                    self.d_gFilter = self.d_gFilter.permute(2, 0, 1, 3).unsqueeze(1)
            else:
                self.d_gFilter = af.interop.np_to_af_array(self.gFilter)
            self.uu = 0
        
        if self.useCUDA:
            if (self.BPType == 5 or self.FPType == 5 or self.FPType == 4) and not self.useCuPy:
                raise ValueError('Unsupported projector for CUDA! Only projector types 1, 11, 12, 13, 2, 21, 23, 3, 31, 32, 14 and 6 are supported')
            self.no_norm = 1
            self.mSize = self.nRowsD * self.nColsD * self.nProjections
            self.d_d = [None] * (self.nMultiVolumes + 1)
            self.d_b = [None] * (self.nMultiVolumes + 1)
            self.d_bmax = [None] * (self.nMultiVolumes + 1)
            self.d_Nxyz = [None] * (self.nMultiVolumes + 1)
            self.dSize = [None] * (self.nMultiVolumes + 1)
            self.d_Scale = [None] * (self.nMultiVolumes + 1)
            self.d_Scale4 = [None] * (self.nMultiVolumes + 1)
            self.d_x = [None] * self.subsets
            self.d_z = [None] * self.subsets
            if self.projector_type != 6:
                if self.useCuPy:
                    # if self.FPType == 5:
                    #     raise ValueError('Not yet supported')
                    self.d_Sens = cp.empty(shape=(1,1), dtype=cp.float32)
                    if (self.listmode == 0 and not self.CT) or self.useIndexBasedReconstruction:
                        self.d_x[0] = cp.asarray(self.x.ravel())
                    elif self.CT and self.listmode == 0:
                        apu = self.x.ravel()
                        for i in range(self.subsets):
                            self.d_x[i] = cp.asarray(apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                    elif self.listamode > 0 and not self.useIndexBasedReconstruction:
                        apu = self.x.ravel()
                        for i in range(self.subsets):
                            if self.loadTOF:
                                self.d_x[i] = cp.asarray(apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                    if (self.CT and self.listmode == 0):
                        if self.pitch:
                            kerroin = 6
                        else:
                            kerroin = 2
                        apu = self.z.ravel()
                        for i in range(self.subsets):
                            self.d_z[i] = cp.asarray(apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                    else:
                        if (self.PET and self.listmode == 0):
                            if self.nLayers > 1:
                                kerroin = 3
                            else:
                                kerroin = 2
                            apu = self.z.ravel()
                            for i in range(self.subsets):
                                self.d_z[i] = cp.asarray(apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                        elif self.listmode == 0 or (self.listmode > 0 and self.useIndexBasedReconstruction):
                            self.d_z[0] = cp.asarray(self.z.ravel())
                        else:
                            for i in range(self.subsets):
                                self.d_z[i] = cp.asarray(np.zeros(1,dtype=np.float32))
                    if (self.attenuation_correction and not self.CTAttenuation):
                        self.d_atten = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_atten[i] = cp.asarray(self.vaimennus[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                    elif (self.attenuation_correction and self.CTAttenuation):
                        if not self.useImages:
                            self.d_atten = cp.asarray(self.vaimennus)
                        else:
                            chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                            array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item())
                            array.copy_from(self.vaimennus.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                            res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                            self.d_atten = cp.cuda.texture.TextureObject(res, tdes)
                    if self.useMaskFP:
                        if not self.useImages:
                            self.d_maskFP = cp.asarray(self.maskFP)
                        else:
                            chl = cp.cuda.texture.ChannelFormatDescriptor(8,0,0,0, cp.cuda.runtime.cudaChannelFormatKindUnsigned)
                            array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD)
                            array.copy_from(self.maskFP.reshape((self.nColsD, self.nRowsD)))
                            res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                            self.d_maskFP = cp.cuda.texture.TextureObject(res, tdes)
                    if self.useMaskBP:
                        if not self.useImages:
                            self.d_maskBP = cp.asarray(self.maskBP)
                        else:
                            chl = cp.cuda.texture.ChannelFormatDescriptor(8,0,0,0, cp.cuda.runtime.cudaChannelFormatKindUnsigned)
                            array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item())
                            array.copy_from(self.maskFP.reshape((self.Ny[0].item(), self.Nx[0].item())))
                            res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                            self.d_maskFP = cp.cuda.texture.TextureObject(res, tdes)
                    if self.TOF:
                        self.d_TOFCenter = cp.asarray(self.TOFCenter)
                    if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                        self.d_V = cp.asarray(self.V)
                    if (self.normalization_correction):
                        self.d_norm = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_norm[i] = cp.asarray(self.normalization[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                    if (self.additionalCorrection):
                        self.d_corr = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_corr[i] = cp.asarray(self.corrVector[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                    if (self.listmode != 1 and ((not self.CT and not self.SPECT and not self.PET) and (self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7)))):
                        self.d_zindex = [None] * self.subsets
                        self.d_xyindex = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_xyindex[i] = cp.asarray(self.xy_index[self.nMeas[i] : self.nMeas[i + 1]])
                            self.d_zindex[i] = cp.asarray(self.z_index[self.nMeas[i] : self.nMeas[i + 1]])
                    if (self.listmode > 0 and self.useIndexBasedReconstruction):
                        self.d_trIndex = [None] * self.subsets
                        self.d_axIndex = [None] * self.subsets
                        for i in range(self.subsets):
                            if self.loadTOF:
                                self.d_trIndex[i] = cp.asarray(self.trIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                                self.d_axIndex[i] = cp.asarray(self.axIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                    if self.OffsetLimit.size > 0 and ((self.BPType == 4 and self.CT) or self.BPType == 5):
                        self.d_T = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_T[i] = cp.asarray(self.offsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
                    mod = cp.RawModule(code=linesFP, options=bOptFP)
                    # import sys
                    # mod.compile(log_stream=sys.stdout)
                    if self.FPType in [1, 2, 3]:
                        self.knlF = mod.get_function('projectorType123')
                    elif self.FPType == 4:
                        self.knlF = mod.get_function('projectorType4Forward')
                    elif self.FPType == 5:
                        self.knlF = mod.get_function('projectorType5Forward')
                    mod = cp.RawModule(code=linesBP, options=bOptBP)
                    if self.BPType in [1, 2, 3]:
                        self.knlB = mod.get_function('projectorType123')
                    elif self.BPType == 4 and not self.CT:
                        self.knlB = mod.get_function('projectorType4Forward')
                    elif self.BPType == 4 and self.CT:
                        self.knlB = mod.get_function('projectorType4Backward')
                    elif self.BPType == 5:
                        self.knlB = mod.get_function('projectorType5Backward')
                    
                    if self.use_psf:
                        with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
                            lines = f.read()
                        lines = hlines + lines
                        bOpt += ('-DCAST=float','-DPSF','-DLOCAL_SIZE=' + str(localSize[0]),'-DLOCAL_SIZE2=' + str(localSize[1]),)
                        mod = cp.RawModule(code=lines, options=bOpt)
                        self.knlPSF = mod.get_function('Convolution3D_f')
                        self.d_gaussPSF = cp.asarray(self.gaussK.ravel('F'))
                        
                    if self.FPType in [1, 2, 3]:
                        self.kIndF = (cp.float32(self.global_factor), cp.float32(self.epps), cp.uint32(self.nRowsD), cp.uint32(self.det_per_ring), cp.float32(self.sigma_x), cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                    elif self.FPType == 4:
                        self.kIndF = (cp.uint32(self.nRowsD), cp.uint32(self.nColsD), cp.float32(self.dPitchX),cp.float32(self.dPitchY),cp.float32(self.dL),cp.float32(self.global_factor),)
                    elif self.FPType == 5:
                        self.kIndF = (cp.uint32(self.nRowsD), cp.uint32(self.nColsD), cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                    if self.FPType in [2,3]:
                        if self.FPType == 2:
                            self.kIndF += (cp.float32(self.tube_width_z),)
                        else:
                            self.kIndF += (cp.float32(self.tube_radius),)
                        self.kIndF += (cp.float32(self.bmin), cp.float32(self.bmax), cp.float32(self.Vmax),)
                    if self.useMaskFP:
                        self.kIndF += (self.d_maskFP,)
                    if self.FPType in [1, 2, 3]:
                        if self.TOF:
                            self.kIndF += (self.d_TOFCenter, )
                        if self.FPType in [2, 3]:
                            self.kIndF += (self.d_V, )
                        self.kIndF += (cp.uint32(self.nColsD),)
                    if self.FPType == 4 and self.TOF:
                        self.kIndF += (self.d_TOFCenter, )
                        self.kIndF += (cp.float32(self.sigma_x), )
                    if self.attenuation_correction and self.CTAttenuation and self.FPType in [1, 2, 3, 4]:
                        self.kIndF += (self.d_atten,)
                        
                    
                    if self.BPType == 4 or self.BPType == 5:
                        self.kIndB = (cp.uint32(self.nRowsD), cp.uint32(self.nColsD), cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                    if self.BPType == 4 and not self.CT:
                        self.kIndB += (cp.float32(self.dL),)
                        self.kIndB += (cp.float32(self.global_factor),)
                    if self.BPType in [1, 2, 3]:
                        self.kIndB = (cp.float32(self.global_factor), cp.float32(self.epps), cp.uint32(self.nRowsD), cp.uint32(self.det_per_ring), cp.float32(self.sigma_x), cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                        if self.BPType in [2, 3]:
                            if self.BPType == 2:
                                self.kIndB  += (cp.float32(self.tube_width_z),)
                            else:
                                self.kIndB  += (cp.float32(self.tube_radius),)
                            self.kIndB += (cp.float32(self.bmin),)
                            self.kIndB += (cp.float32(self.bmax),)
                            self.kIndB += (cp.float32(self.Vmax),)
                        if self.useMaskFP:
                            self.kIndB += (self.d_maskFP,)
                    if self.useMaskBP:
                        self.kIndB += (self.d_maskBP,)
                    if self.BPType in [1, 2, 3]:
                        if self.TOF:
                            self.kIndB += (self.d_TOFCenter,)
                        if self.BPType in [2, 3]:
                            self.kIndB += (self.d_V,)
                        self.kIndB += (cp.uint32(self.nColsD),)
                    if self.BPType == 4 and not self.CT and self.TOF:
                        self.kIndB += (self.d_TOFCenter,)
                        self.kIndB += (cp.float32(self.sigma_x),)
                    if self.attenuation_correction and self.CTAttenuation and self.BPType in [1, 2, 3, 4] and not self.CT:
                        self.kIndB += (self.d_atten,)
                else:
                    import pycuda as cuda
                    from pycuda.compiler import SourceModule
                    import pycuda.autoinit
                    self.d_Sens = cuda.gpuarray.empty(shape=(1,1), dtype=np.float32)
                    for k in range(self.nMultiVolumes + 1):
                        self.d_d[k] = cuda.gpuarray.vec.make_float3(self.dx[k].item(), self.dy[k].item(), self.dz[k].item())
                        self.d_b[k] = cuda.gpuarray.vec.make_float3(self.bx[k].item(), self.by[k].item(), self.bz[k].item())
                        self.d_bmax[k] = cuda.gpuarray.vec.make_float3(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item(), self.by[k].item() + self.Ny[k].item() * self.dy[k].item(), self.bz[k].item() + self.Nz[k].item() * self.dz[k].item())
                        self.d_Nxyz[k] = cuda.gpuarray.vec.make_uint3(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
                        if (self.FPType == 4 or self.FPType == 5 or self.BPType == 4 or self.BPType == 5):
                            self.d_Scale4[k] = cuda.gpuarray.vec.make_float3(self.dScaleX4[k].item(), self.dScaleY4[k].item(), self.dScaleZ4[k].item())
                            if self.FPType == 5 or self.BPType == 5:
                                self.dSize[k] = cuda.gpuarray.vec.make_float2(self.dSizeX[k].item(), self.dSizeY[k].item())
                                self.d_Scale[k] = cuda.gpuarray.vec.make_float3(self.dScaleX[k].item(), self.dScaleY[k].item(), self.dScaleZ[k].item())
                                if k == 0:
                                    self.dSizeBP = cuda.gpuarray.vec.make_float2(self.dSizeXBP, self.dSizeZBP)
                    self.d_dPitch = cuda.gpuarray.vec.make_float2(self.dPitchX, self.dPitchY)
                    if (self.listmode == 0 and not self.CT) or self.useIndexBasedReconstruction:
                        self.d_x[0] = cuda.gpuarray.to_gpu(self.x.ravel())
                    elif self.CT and self.listmode == 0:
                        apu = self.x.ravel()
                        for i in range(self.subsets):
                            self.d_x[i] = cuda.gpuarray.to_gpu(apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                    elif self.listamode > 0 and not self.useIndexBasedReconstruction:
                        apu = self.x.ravel()
                        for i in range(self.subsets):
                            if (i == 0 or self.loadTOF):
                                self.d_x[i] = cuda.gpuarray.to_gpu(apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                    if (self.CT and self.listmode == 0):
                        if self.pitch:
                            kerroin = 6
                        else:
                            kerroin = 2
                        apu = self.z.ravel()
                        for i in range(self.subsets):
                            self.d_z[i] = cuda.gpuarray.to_gpu(apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                    else:
                        if (self.PET and self.listmode == 0):
                            if self.nLayers > 1:
                                kerroin = 3
                            else:
                                kerroin = 2
                            apu = self.z.ravel()
                            for i in range(self.subsets):
                                self.d_z[i] = cuda.gpuarray.to_gpu(apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                        elif self.listmode == 0:
                            self.d_z[0] = cuda.gpuarray.to_gpu(self.z.ravel())
                    if (self.attenuation_correction and not self.CTAttenuation):
                        self.d_atten = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_atten[i] = cuda.gpuarray.to_gpu(self.vaimennus[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                    elif (self.attenuation_correction and self.CTAttenuation):
                        self.d_atten = cuda.gpuarray.to_gpu(self.vaimennus)
                    if self.useMaskFP:
                        self.d_maskFP = cuda.gpuarray.to_gpu(self.maskFP)
                    if self.useMaskBP:
                        self.d_maskBP = cuda.gpuarray.to_gpu(self.maskBP)
                    if self.TOF:
                        self.d_TOFCenter = cuda.gpuarray.to_gpu(self.TOFCenter)
                    if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                        self.d_V = cuda.gpuarray.to_gpu(self.V)
                    if (self.normalization_correction):
                        self.d_norm = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_norm[i] = cuda.gpuarray.to_gpu(self.normalization[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                    if (self.additionalCorrection):
                        self.d_corr = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_corr[i] = cuda.gpuarray.to_gpu(self.corrVector[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                    if (self.listmode != 1 and ((not self.CT and not self.SPECT and not self.PET) and (self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7)))):
                        self.d_zindex = [None] * self.subsets
                        self.d_xyindex = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_xyindex[i] = cuda.gpuarray.to_gpu(self.xy_index[self.nMeas[i] : self.nMeas[i + 1]])
                            self.d_zindex[i] = cuda.gpuarray.to_gpu(self.z_index[self.nMeas[i] : self.nMeas[i + 1]])
                    if self.OffsetLimit.size > 0 and ((self.BPType == 4 and self.CT) or self.BPType == 5):
                        self.d_T = [None] * self.subsets
                        for i in range(self.subsets):
                            self.d_T[i] = cuda.gpuarray.to_gpu(self.offsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
                    try:
                        mod = SourceModule(linesFP, options=list(bOptFP), no_extern_c=True)
                    except cuda.driver.CompileError as e:
                        print("Compilation error:", e)
                        raise ValueError("Forward projection compilation failed!")
                    if self.FPType in [1, 2, 3]:
                        self.knlF = mod.get_function('projectorType123')
                    elif self.FPType == 4:
                        self.knlF = mod.get_function('projectorType4Forward')
                    elif self.FPType == 5:
                        self.knlF = mod.get_function('projectorType5Forward')
                    try:
                        mod = SourceModule(linesBP, options=list(bOptBP), no_extern_c=True)
                    except cuda.driver.CompileError as e:
                        print("Compilation error:", e)
                        raise ValueError("Backprojection compilation failed!")
                    if self.BPType in [1, 2, 3]:
                        self.knlB = mod.get_function('projectorType123')
                    elif self.BPType == 4 and not self.CT:
                        self.knlB = mod.get_function('projectorType4Forward')
                    elif self.BPType == 4 and self.CT:
                        self.knlB = mod.get_function('projectorType4Backward')
                    elif self.BPType == 5:
                        self.knlB = mod.get_function('projectorType5Backward')
                    self.kIndF = (np.float32(self.global_factor), np.float32(self.epps), np.uint32(self.nRowsD), np.uint32(self.det_per_ring), np.float32(self.sigma_x), np.float32(self.d_dPitch['x'].item()),np.float32(self.d_dPitch['y'].item()),)
                    
                    if self.use_psf:
                        with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
                            lines = f.read()
                        lines = hlines + lines
                        bOpt += ('-DCAST=float','-DPSF','-DLOCAL_SIZE=' + str(localSize[0]),'-DLOCAL_SIZE2=' + str(localSize[1]),)
                        mod = SourceModule(lines, options=list(bOpt), no_extern_c=True)
                        self.knlPSF = mod.get_function('Convolution3D_f')
                        self.d_gaussPSF = cuda.gpuarray.to_gpu(self.gaussK.ravel('F'))
                    
                    if self.FPType in [2,3]:
                        if self.FPTye == 2:
                            self.kIndF += (np.float32(self.tube_width_z),)
                        else:
                            self.kIndF += (np.float32(self.tube_radius),)
                        self.kIndF += (np.float32(self.bmin), np.float32(self.bmax), np.float32(self.Vmax),)
                    if self.useMaskFP:
                        self.kIndF += (self.d_maskFP.gpudata,)
                    if self.FPType in [1, 2, 3]:
                        if self.TOF:
                            self.kIndF += (self.d_TOFCenter.gpudata, )
                        if self.FPType in [2, 3]:
                            self.kIndF += (self.d_V.gpudata, )
                        self.kIndF += (np.uint32(self.nColsD),)
                    if self.attenuation_correction and self.CTAttenuation and self.FPType in [1, 2, 3, 4]:
                        self.kIndF += (self.d_atten.gpudata,)
                        
                    
                    if self.BPType == 4 or self.BPType == 5:
                        self.kIndB = (np.uint32(self.nRowsD), np.uint32(self.nColsD), np.float32(self.d_dPitch['x'].item()),np.float32(self.d_dPitch['y'].item()),)
                    if self.BPType == 4 and not self.CT:
                        self.kIndB += (np.float32(self.dL),)
                        self.kIndB += (np.float32(self.global_factor),)
                    if self.BPType in [1, 2, 3]:
                        self.kIndB = (np.float32(self.global_factor), np.float32(self.epps), np.uint32(self.nRowsD), np.uint32(self.det_per_ring), np.float32(self.sigma_x), np.float32(self.d_dPitch['x'].item()),np.float32(self.d_dPitch['y'].item()),)
                        if self.BPType in [2, 3]:
                            if self.BPType == 2:
                                self.kIndB  += (np.float32(self.tube_width_z),)
                            else:
                                self.kIndB  += (np.float32(self.tube_radius),)
                            self.kIndB += (np.float32(self.bmin),)
                            self.kIndB += (np.float32(self.bmax),)
                            self.kIndB += (np.float32(self.Vmax),)
                        if self.useMaskFP:
                            self.kIndB += (self.d_maskFP.gpudata,)
                    if self.useMaskBP:
                        self.kIndB += (self.d_maskBP.gpudata,)
                    if self.BPType in [1, 2, 3]:
                        if self.TOF:
                            self.kIndB += (self.d_TOFCenter.gpudata,)
                        if self.BPType in [2, 3]:
                            self.kIndB += (self.d_V.gpudata,)
                        self.kIndB += (np.uint32(self.nColsD),)
                    if self.BPType == 4 and not self.CT and self.TOF:
                        self.kIndB += (self.d_TOFCenter.gpudata,)
                        self.kIndB += (np.float32(self.sigma_x),)
                    if self.attenuation_correction and self.CTAttenuation and self.BPType in [1, 2, 3, 4] and not self.CT:
                        self.kIndB += (self.d_atten.gpudata,)
        else:
            if self.projector_type != 6:
                import pyopencl as cl
                
                if self.useAF:
                    ctx = af.opencl.get_context(retain=True)
                    self.clctx = cl.Context.from_int_ptr(ctx)
                    q = af.opencl.get_queue(True)
                    self.queue = cl.CommandQueue.from_int_ptr(q)
                else:
                    platforms = cl.get_platforms()
                    dList = platforms[self.platform].get_devices()
                    dList = [dList[self.deviceNum]]
                    self.clctx = cl.Context(devices=dList)
                    self.queue = cl.CommandQueue(self.clctx)
                
                self.no_norm = 1
                self.mSize = self.nRowsD * self.nColsD * self.nProjections
                # length = [None] * self.subsets
                # for i in range(self.subsets):
                #     if self.subsetType >= 8:
                self.d_Sens = cl.array.empty(self.queue, shape=(1,1), dtype=cl.cltypes.float)
                
                
                self.d_d = [None] * (self.nMultiVolumes + 1)
                self.d_b = [None] * (self.nMultiVolumes + 1)
                self.d_bmax = [None] * (self.nMultiVolumes + 1)
                self.d_Nxyz = [None] * (self.nMultiVolumes + 1)
                self.dSize = [None] * (self.nMultiVolumes + 1)
                self.d_Scale = [None] * (self.nMultiVolumes + 1)
                self.d_Scale4 = [None] * (self.nMultiVolumes + 1)
                for k in range(self.nMultiVolumes + 1):
                    self.d_d[k] = cl.cltypes.make_float3(self.dx[k].item(), self.dy[k].item(), self.dz[k].item())
                    self.d_b[k] = cl.cltypes.make_float3(self.bx[k].item(), self.by[k].item(), self.bz[k].item())
                    self.d_bmax[k] = cl.cltypes.make_float3(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item(), self.by[k].item() + self.Ny[k].item() * self.dy[k].item(), self.bz[k].item() + self.Nz[k].item() * self.dz[k].item())
                    self.d_Nxyz[k] = cl.cltypes.make_uint3(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
                    if (self.FPType == 4 or self.FPType == 5 or self.BPType == 4 or self.BPType == 5):
                        self.d_Scale4[k] = cl.cltypes.make_float3(self.dScaleX4[k].item(), self.dScaleY4[k].item(), self.dScaleZ4[k].item())
                        if self.FPType == 5 or self.BPType == 5:
                            self.dSize[k] = cl.cltypes.make_float2(self.dSizeX[k].item(), self.dSizeY[k].item())
                            self.d_Scale[k] = cl.cltypes.make_float3(self.dScaleX[k].item(), self.dScaleY[k].item(), self.dScaleZ[k].item())
                            if k == 0:
                                self.dSizeBP = cl.cltypes.make_float2(self.dSizeXBP, self.dSizeZBP)
                self.d_dPitch = cl.cltypes.make_float2(self.dPitchX, self.dPitchY)
                self.d_x = [None] * self.subsets
                self.d_z = [None] * self.subsets
                if (self.listmode == 0 and not self.CT) or self.useIndexBasedReconstruction:
                    self.d_x[0] = cl.array.to_device(self.queue, self.x.ravel())
                elif self.CT and self.listmode == 0:
                    apu = self.x.ravel()
                    for i in range(self.subsets):
                        self.d_x[i] = cl.array.to_device(self.queue, apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                elif self.listmode > 0 and not self.useIndexBasedReconstruction:
                    apu = self.x.ravel()
                    for i in range(self.subsets):
                        if self.loadTOF:
                            self.d_x[i] = cl.array.to_device(self.queue, apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                if (self.CT and self.listmode == 0):
                    if self.pitch:
                        kerroin = 6
                    else:
                        kerroin = 2
                    apu = self.z.ravel()
                    for i in range(self.subsets):
                        self.d_z[i] = cl.array.to_device(self.queue, apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                else:
                    if (self.PET and self.listmode == 0):
                        if self.nLayers > 1:
                            kerroin = 3
                        else:
                            kerroin = 2
                        apu = self.z.ravel()
                        for i in range(self.subsets):
                            self.d_z[i] = cl.array.to_device(self.queue, apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                    elif self.listmode == 0 or (self.listmode > 0 and self.useIndexBasedReconstruction):
                        self.d_z[0] = cl.array.to_device(self.queue, self.z.ravel())
                    else:
                        for i in range(self.subsets):
                            self.d_z[i] = cl.array.to_device(self.queue, np.zeros(1,dtype=np.float32))
                if (self.attenuation_correction and not self.CTAttenuation):
                    self.d_atten = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_atten[i] = cl.array.to_device(self.queue, self.vaimennus[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                elif (self.attenuation_correction and self.CTAttenuation):
                    imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
                    self.d_atten = cl.Image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.vaimennus, shape=(self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()))
                    # self.d_atten = cl.image_from_array(self.clctx, np.reshape(self.vaimennus, (self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()), order='F'))
                if self.useMaskFP:
                    self.d_maskFP = cl.image_from_array(self.clctx, self.maskFP)
                if self.useMaskBP:
                    self.d_maskBP = cl.image_from_array(self.clctx, self.maskBP)
                if self.TOF:
                    self.d_TOFCenter = cl.array.to_device(self.queue, self.TOFCenter)
                if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                    self.d_V = cl.array.to_device(self.queue, self.V)
                if (self.normalization_correction):
                    self.d_norm = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_norm[i] = cl.array.to_device(self.queue, self.normalization[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                if (self.additionalCorrection):
                    self.d_corr = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_corr[i] = cl.array.to_device(self.queue, self.corrVector[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                if (self.listmode != 1 and ((not self.CT and not self.SPECT and not self.PET) and (self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7)))):
                    self.d_zindex = [None] * self.subsets
                    self.d_xyindex = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_xyindex[i] = cl.array.to_device(self.queue, self.xy_index[self.nMeas[i] : self.nMeas[i + 1]])
                        self.d_zindex[i] = cl.array.to_device(self.queue, self.z_index[self.nMeas[i] : self.nMeas[i + 1]])
                if (self.listmode > 0 and self.useIndexBasedReconstruction):
                    self.d_trIndex = [None] * self.subsets
                    self.d_axIndex = [None] * self.subsets
                    for i in range(self.subsets):
                        if self.loadTOF:
                            self.d_trIndex[i] = cl.array.to_device(self.queue, self.trIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                            self.d_axIndex[i] = cl.array.to_device(self.queue, self.axIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                if self.OffsetLimit.size > 0 and ((self.BPType == 4 and self.CT) or self.BPType == 5):
                    self.d_T = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_T[i] = cl.array.to_device(self.queue, self.offsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
                # d_Sens = cl.Buffer(clctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=Sens)
                # d_x = cl.Buffer(self.clctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.x)
                # z = cl.Buffer(clctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.z)
                prg = cl.Program(self.clctx, linesFP).build(' '.join(bOptFP))
                if self.FPType in [1, 2, 3]:
                    self.knlF = prg.projectorType123
                elif self.FPType == 4:
                    self.knlF = prg.projectorType4Forward
                elif self.FPType == 5:
                    self.knlF = prg.projectorType5Forward
                prg = cl.Program(self.clctx, linesBP).build(' '.join(bOptBP))
                if self.BPType in [1, 2, 3]:
                    self.knlB = prg.projectorType123
                elif self.BPType == 4 and not self.CT:
                    self.knlB = prg.projectorType4Forward
                elif self.BPType == 4 and self.CT:
                    self.knlB = prg.projectorType4Backward
                elif self.BPType == 5:
                    self.knlB = prg.projectorType5Backward
                
                if self.use_psf:
                    with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
                        lines = f.read()
                    lines = hlines + lines
                    bOpt +=(' -DCAST=float',' -DPSF',' -DLOCAL_SIZE=' + str(localSize[0]), ' -DLOCAL_SIZE2=' + str(localSize[1]),)
                    prg = cl.Program(self.clctx, lines).build(' '.join(bOpt))
                    self.knlPSF = prg.Convolution3D_f
                    self.d_gaussPSF = cl.array.to_device(self.queue, self.gaussK.ravel('F'))
                    
                self.kIndF = 0
                if self.FPType == 4 or self.FPType == 5:
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.nRowsD))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.nColsD))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, self.d_dPitch)
                    self.kIndF += 1
                if self.FPType == 4:
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.dL))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.global_factor))
                    self.kIndF += 1
                if self.FPType in [1, 2, 3]:
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.global_factor))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.epps))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.nRowsD))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.det_per_ring))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.sigma_x))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, self.d_dPitch)
                    self.kIndF += 1
                    if self.FPType in [2, 3]:
                        if self.FPType == 2:
                            self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.tube_width_z))
                            self.kIndF += 1
                        else:
                            self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.tube_radius))
                            self.kIndF += 1
                        self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.bmin))
                        self.kIndF += 1
                        self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.bmax))
                        self.kIndF += 1
                        self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.Vmax))
                        self.kIndF += 1
                if self.useMaskFP:
                    self.knlF.set_arg(self.kIndF, self.d_maskFP)
                    self.kIndF += 1
                if self.FPType in [1, 2, 3]:
                    if self.TOF:
                        self.knlF.set_arg(self.kIndF, self.d_TOFCenter.data)
                        self.kIndF += 1
                    if self.FPType in [2, 3]:
                        self.knlF.set_arg(self.kIndF, self.d_V.data)
                        self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.nColsD))
                    self.kIndF += 1
                if self.FPType == 4 and not self.CT and self.TOF:
                    self.knlF.set_arg(self.kIndF, self.d_TOFCenter.data)
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.sigma_x))
                    self.kIndF += 1
                if self.attenuation_correction and self.CTAttenuation and self.FPType in [1, 2, 3, 4]:
                    self.knlF.set_arg(self.kIndF, self.d_atten)
                    self.kIndF += 1
                        
                    
                
                self.kIndB = 0
                if self.BPType == 4 or self.BPType == 5:
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.nRowsD))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.nColsD))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, self.d_dPitch)
                    self.kIndB += 1
                if self.BPType == 4 and not self.CT:
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.dL))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.global_factor))
                    self.kIndB += 1
                if self.BPType in [1, 2, 3]:
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.global_factor))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.epps))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.nRowsD))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.det_per_ring))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.sigma_x))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, self.d_dPitch)
                    self.kIndB += 1
                    if self.BPType in [2, 3]:
                        if self.BPType == 2:
                            self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.tube_width_z))
                            self.kIndB += 1
                        else:
                            self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.tube_radius))
                            self.kIndB += 1
                        self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.bmin))
                        self.kIndB += 1
                        self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.bmax))
                        self.kIndB += 1
                        self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.Vmax))
                        self.kIndB += 1
                    if self.useMaskFP:
                        self.knlB.set_arg(self.kIndB, self.d_maskFP)
                        self.kIndB += 1
                if self.useMaskBP:
                    self.knlB.set_arg(self.kIndB, self.d_maskBP)
                    self.kIndB += 1
                if self.BPType in [1, 2, 3]:
                    if self.TOF:
                        self.knlB.set_arg(self.kIndB, self.d_TOFCenter.data)
                        self.kIndB += 1
                    if self.BPType in [2, 3]:
                        self.knlB.set_arg(self.kIndB, self.d_V.data)
                        self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.nColsD))
                    self.kIndB += 1
                if self.BPType == 4 and not self.CT and self.TOF:
                    self.knlB.set_arg(self.kIndB, self.d_TOFCenter.data)
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.sigma_x))
                    self.kIndB += 1
                if self.attenuation_correction and self.CTAttenuation and self.BPType in [1, 2, 3, 4] and not self.CT:
                    self.knlB.set_arg(self.kIndB, self.d_atten)
                    self.kIndB += 1
                
                
            
    def computeConvolution(self, f, ii = 0):
        globalSize = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
        kInd = 0
        if self.useCUDA:
            if self.useTorch:
                import torch
            if self.useCuPy:
                import cupy as cp
            else:
                import pycuda as cuda
        else:
            import pyopencl as cl
        if self.useAF:
            import arrayfire as af
            if isinstance(f, af.array.Array):
                ptr = f.raw_ptr()
                f = cl.MemoryObject.from_int_ptr(ptr)
                self.knlPSF.set_arg(kInd, f)
            else:
                self.knlPSF.set_arg(kInd, f.data)
        else:
            if self.useCUDA:
                if self.useTorch:
                    output = torch.zeros(self.N[ii].item(), dtype=torch.float32, device='cuda')
                elif self.useCuPy:
                    output = cp.zeros(self.N[ii].item(), dtype=cp.float32)
                else:
                    output = cuda.gpuarray.zeros(self.N[ii].item(), dtype=np.float32)
            else:
                output = cl.array.zeros(self.queue, self.N[ii].item(), dtype=cl.cltypes.float)
            if not self.useCUDA:
                self.knlPSF.set_arg(kInd, f.data)
        if self.useCUDA:
            if self.useCuPy:
                fD = cp.asarray(f)
                outputD = cp.asarray(output)
                self.knlPSF((globalSize[0] // 16, globalSize[1] // 16, globalSize[2] // 16), (16,16,1),(fD, outputD, self.d_gaussPSF, cp.int32(self.g_dim_x), cp.int32(self.g_dim_y), cp.int32(self.g_dim_z)))
            else:
                if self.useTorch:
                    class Holder(cuda.driver.PointerHolderBase):
                        def __init__(self, t):
                            super(Holder, self).__init__()
                            self.t = t
                            self.gpudata = t.data_ptr()
                        def get_pointer(self):
                            return self.t.data_ptr()
                    fD = Holder(f)
                    outputD = Holder(output)
                    self.knlPSF(fD, outputD, self.d_gaussPSF.gpudata, np.int32(self.g_dim_x), np.int32(self.g_dim_y), np.int32(self.g_dim_z), block=(16,16,1), grid=(globalSize[0], globalSize[1], globalSize[2]))
                else:
                    self.knlPSF(f.gpudata, output.gpudata, self.d_gaussPSF.gpudata, np.int32(self.g_dim_x), np.int32(self.g_dim_y), np.int32(self.g_dim_z), block=(16,16,1), grid=(globalSize[0], globalSize[1], globalSize[2]))
            if self.useTorch:
                torch.cuda.synchronize()
        else:
            kInd += 1
            if self.useAF:
                output = af.data.constant(0., self.N[ii].item(), dtype=af.Dtype.f32)
                outPtr = cl.MemoryObject.from_int_ptr(output.raw_ptr())
                self.knlPSF.set_arg(kInd, outPtr)
            else:
                self.knlPSF.set_arg(kInd, output.data)
            kInd += 1
            self.knlPSF.set_arg(kInd, self.d_gaussPSF.data)
            kInd += 1
            self.knlPSF.set_arg(kInd, (cl.cltypes.int)(self.g_dim_x))
            kInd += 1
            self.knlPSF.set_arg(kInd, (cl.cltypes.int)(self.g_dim_y))
            kInd += 1
            self.knlPSF.set_arg(kInd, (cl.cltypes.int)(self.g_dim_z))
            cl.enqueue_nd_range_kernel(self.queue, self.knlPSF, globalSize, (16, 16, 1))
            self.queue.finish()
        if self.useAF:
            af.device.unlock_array(output)
        return output

    def forwardProject(self, f, subset = -1):
        if subset == -1:
            subset = self.subset
        if self.projector_type == 6:
            if not self.useCUDA:
                import arrayfire as af
                u1 = np.sum(self.nProjSubset[:subset])
                y = af.data.constant(0., self.nRowsD, self.nColsD, self.nProjSubset[subset].item())
                for ii in range(self.nMultiVolumes + 1):
                    if isinstance(f,list):
                        apuArr = af.data.moddims(f[ii], self.Nx[ii].item(), self.Ny[ii].item(), self.Nz[ii].item())
                    else:
                        apuArr = af.data.moddims(f, self.Nx[ii].item(), self.Ny[ii].item(), self.Nz[ii].item())
                    for kk in range(self.nProjSubset[subset].item()):
                        kuvaRot = af.image.rotate(apuArr, -self.angles[u1].item(), method=af.INTERP.BILINEAR) # [128, 128, 96]
                        kuvaRot = af.data.reorder(kuvaRot, 2, 1, 0) # [96, 128, 128]
                        kuvaRot = af.signal.convolve2(kuvaRot, self.d_gFilter[:, :, :, u1]) # [96, 128, 128]
                        kuvaRot = kuvaRot[:, :, self.blurPlanes[u1].item():] # [96, 128, 90]
                        kuvaRot = af.data.reorder(kuvaRot, 2, 1, 0) # [90, 128, 96]
                        kuvaRot = af.data.reorder(af.algorithm.sum(kuvaRot, 0), 1, 2, 0) # [1, 128, 96] -->  [128, 96]
                        y[:, :, kk] += kuvaRot
                        u1 += 1
                y = af.flat(y)
                af.eval(y)
            else:
                if not self.useCuPy:
                    raise ValueError('SPECT is only supported on CUDA with CuPy!')
                else:
                    import torch
                    from torchvision.transforms.functional import rotate
                    from torchvision.transforms import InterpolationMode
                    import torch.nn.functional as F
                    u1 = np.sum(self.nProjSubset[:subset])
                    y = torch.zeros((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD), dtype=torch.float32).cuda()
                    for ii in range(self.nMultiVolumes + 1):
                        if isinstance(f,list):
                            apuArr = torch.reshape(f[ii], (self.Nz[ii].item(), self.Ny[ii].item(), self.Nx[ii].item()))
                        else:
                            apuArr = torch.reshape(f, (self.Nz[ii].item(), self.Ny[ii].item(), self.Nx[ii].item()))
                        for kk in range(self.nProjSubset[subset].item()):
                            kuvaRot = rotate(apuArr, -self.angles[u1].item(), InterpolationMode.BILINEAR)
                            kuvaRot = torch.permute(kuvaRot, (2, 1, 0))
                            kuvaRot = kuvaRot.unsqueeze(0)
                            kuvaRot = F.conv2d(kuvaRot, self.d_gFilter[:, :, :, :, u1], padding=(self.d_gFilter.shape[2] // 2, self.d_gFilter.shape[3] // 2), groups=kuvaRot.shape[1])
                            kuvaRot = kuvaRot.squeeze(0)
                            kuvaRot = kuvaRot[self.blurPlanes[u1].item():, :, :]
                            kuvaRot = torch.sum(kuvaRot, 0)
                            kuvaRot = torch.permute(kuvaRot, (1, 0))
                            y[kk, :, :] += kuvaRot
                            u1 += 1
                    y = y.ravel()
        else:
            if self.useCUDA:
                if self.useCuPy:
                    import cupy as cp
                    if not self.loadTOF:
                        if self.useIndexBasedReconstruction and self.listmode > 0:
                            self.d_trIndex[0] = cp.asarray(self.trIndex[self.nMeas[subset] * 2 : self.nMeas[subset + 1] * 2])
                            self.d_axIndex[0] = cp.asarray(self.axIndex[self.nMeas[subset] * 2 : self.nMeas[subset + 1] * 2])
                        elif self.listmode > 0:
                            apu = self.x.ravel()
                            self.d_x[0] = cp.asarray(apu[self.nMeas[subset] * 6 : self.nMeas[subset + 1] * 6])
                    if self.useTorch:
                        import torch
                        if self.subsetType > 7 or self.subsets == 1:
                            y = torch.zeros(self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=torch.float32, device='cuda')
                        else:
                            y = torch.zeros(self.nMeasSubset[subset].item(), dtype=torch.float32, device='cuda')
                        yD = cp.asarray(y)
                    else:
                        if self.subsetType > 7 or self.subsets == 1:
                            y = cp.zeros(self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=cp.float32)
                        else:
                            y = cp.zeros(self.nMeasSubset[subset].item(), dtype=cp.float32)
                    for k in range(self.nMultiVolumes + 1):
                        if isinstance(f,list):
                            if self.use_psf:
                                f[k] = self.computeConvolution(f[k])
                            if self.useTorch:
                                fD = cp.asarray(f[k])
                        else:
                            if self.use_psf:
                                f = self.computeConvolution(f)
                            if self.useTorch:
                                fD = cp.asarray(f)
                        if self.FPType == 5:
                            intIm = cp.zeros((self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()), dtype=cp.float32, order='F')
                            if self.useTorch:
                                intIm[1:,1:,:] = cp.transpose(fD.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (1, 2, 0))
                            else:
                                if isinstance(f,list):
                                    intIm[1:,1:,:] = cp.transpose(f[k].reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (1, 2, 0))
                                else:
                                    intIm[1:,1:,:] = cp.transpose(f.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (1, 2, 0))
                            intIm = intIm.cumsum(0)
                            intIm = intIm.cumsum(1)
                            intIm = intIm.ravel('F')
                            chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                            array = cp.cuda.texture.CUDAarray(chl, self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item())
                            array.copy_from(intIm.reshape((self.Nx[k].item(), self.Nz[k].item() + 1, self.Ny[k].item() + 1)))
                            res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                            ff2 = cp.cuda.texture.TextureObject(res, tdes)
                            intIm = cp.zeros((self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()), dtype=cp.float32, order='F')
                            if self.useTorch:
                                intIm[1:,1:,:] = cp.transpose(fD.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (0, 2, 1))
                            else:
                                if isinstance(f,list):
                                    intIm[1:,1:,:] = cp.transpose(f[k].reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (0, 2, 1))
                                else:
                                    intIm[1:,1:,:] = cp.transpose(f.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (0, 2, 1))
                            intIm = intIm.cumsum(0)
                            intIm = intIm.cumsum(1)
                            intIm = intIm.ravel('F')
                            array2 = cp.cuda.texture.CUDAarray(chl, self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item())
                            array2.copy_from(intIm.reshape((self.Ny[k].item(), self.Nz[k].item() + 1, self.Nx[k].item() + 1)))
                            res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array2)
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                            ff = cp.cuda.texture.TextureObject(res, tdes)
                        kIndLoc = self.kIndF
                        if self.FPType == 1 or self.FPType == 2 or self.FPType == 3 or self.FPType == 4:
                            if (self.attenuation_correction and not self.CTAttenuation):
                                kIndLoc += (self.d_atten[subset],)
                        if self.FPType == 5 or self.FPType == 4:
                            kIndLoc += (cp.uint32(self.Nx[k].item()),)
                            kIndLoc += (cp.uint32(self.Ny[k].item()),)
                            kIndLoc += (cp.uint32(self.Nz[k].item()),)
                            kIndLoc += (cp.float32(self.bx[k].item()),)
                            kIndLoc += (cp.float32(self.by[k].item()),)
                            kIndLoc += (cp.float32(self.bz[k].item()),)
                            if self.FPType == 5:
                                kIndLoc += (cp.float32(self.dSizeX[k].item()),)
                                kIndLoc += (cp.float32(self.dSizeY[k].item()),)
                                kIndLoc += (cp.float32(self.dx[k].item()),)
                                kIndLoc += (cp.float32(self.dy[k].item()),)
                                kIndLoc += (cp.float32(self.dz[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleX[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleY[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleZ[k].item()),)
                            else:
                                kIndLoc += (cp.float32(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item()),)
                                kIndLoc += (cp.float32(self.by[k].item() + self.Ny[k].item() * self.dy[k].item()),)
                                kIndLoc += (cp.float32(self.bz[k].item() + self.Nz[k].item() * self.dz[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleX4[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleY4[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleZ4[k].item()),)
                        if self.FPType == 4:
                            if isinstance(f,list):
                                chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                array = cp.cuda.texture.CUDAarray(chl, self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
                                if self.useTorch:
                                    array.copy_from(fD.reshape((self.Nz[k].item(), self.Ny[k].item(), self.Nx[k].item())))
                                else:
                                    array.copy_from(f[k].reshape((self.Nz[k].item(), self.Ny[k].item(), self.Nx[k].item())))
                                res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeBorder, cp.cuda.runtime.cudaAddressModeBorder,cp.cuda.runtime.cudaAddressModeBorder), 
                                                                        filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                                ff = cp.cuda.texture.TextureObject(res, tdes)
                                kIndLoc += (ff,)
                            else:
                                chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item())
                                if self.useTorch:
                                    array.copy_from(fD.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                                else:
                                    array.copy_from(f.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                                res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeBorder, cp.cuda.runtime.cudaAddressModeBorder,cp.cuda.runtime.cudaAddressModeBorder), 
                                                                        filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                                ff = cp.cuda.texture.TextureObject(res, tdes)
                                kIndLoc += (ff,)
                            if self.useTorch:
                                kIndLoc += (yD,)
                            else:
                                kIndLoc += (y,)
                            if (self.listmode == 0 and not self.CT):
                                kIndLoc += (self.d_x[0],)
                            else:
                                kIndLoc += (self.d_x[subset], )
                            if (self.CT or self.PET or self.listmode > 0):
                                kIndLoc += (self.d_z[subset],)
                            else:
                                kIndLoc += (self.d_z[0],)
                            kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                kIndLoc += (self.d_xyindex[subset],)
                                kIndLoc += (self.d_zindex[subset],)
                            if (self.normalization_correction):
                                kIndLoc += (self.d_norm[subset],)
                            if (self.additionalCorrection):
                                kIndLoc += (self.d_corr[subset],)
                            kIndLoc += (cp.uint8(self.no_norm),)
                            kIndLoc += (cp.uint64(self.nMeasSubset[subset].item()),)
                            kIndLoc += (cp.uint32(subset),)
                            kIndLoc += (cp.int32(k),)
                        elif self.FPType == 5:
                            kIndLoc += (self.d_x[subset], )
                            kIndLoc += (self.d_z[subset],)
                            # self.knlF.set_arg(kIndLoc, d_im)
                            # kIndLoc += 1
                            # self.knlF.set_arg(kIndLoc, d_imInt)
                            kIndLoc += (ff,)
                            kIndLoc += (ff2,)
                            if self.useTorch:
                                kIndLoc += (yD,)
                            else:
                                kIndLoc += (y,)
                            kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                            # if self.meanFP:
                        elif self.FPType in [1, 2, 3]:
                            if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                                kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                            if (((self.listmode == 0 and not self.CT) or self.useIndexBasedReconstruction)) or (not self.loadTOF and self.listmode > 0):
                                kIndLoc += (self.d_x[0],)
                            else:
                                kIndLoc += (self.d_x[subset], )
                            if (self.CT or self.PET or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                                kIndLoc += (self.d_z[subset],)
                            else:
                                kIndLoc += (self.d_z[0],)
                            if (self.normalization_correction):
                                kIndLoc += (self.d_norm[subset],)
                            if (self.additionalCorrection):
                                kIndLoc += (self.d_corr[subset],)
                            kIndLoc += (self.d_Sens,)
                            kIndLoc += (cp.uint32(self.Nx[k].item()),)
                            kIndLoc += (cp.uint32(self.Ny[k].item()),)
                            kIndLoc += (cp.uint32(self.Nz[k].item()),)
                            kIndLoc += (cp.float32(self.dx[k].item()),)
                            kIndLoc += (cp.float32(self.dy[k].item()),)
                            kIndLoc += (cp.float32(self.dz[k].item()),)
                            kIndLoc += (cp.float32(self.bx[k].item()),)
                            kIndLoc += (cp.float32(self.by[k].item()),)
                            kIndLoc += (cp.float32(self.bz[k].item()),)
                            kIndLoc += (cp.float32(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item()),)
                            kIndLoc += (cp.float32(self.by[k].item() + self.Ny[k].item() * self.dy[k].item()),)
                            kIndLoc += (cp.float32(self.bz[k].item() + self.Nz[k].item() * self.dz[k].item()),)
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                kIndLoc += (self.d_xyindex[subset],)
                                kIndLoc += (self.d_zindex[subset],)
                            if self.useIndexBasedReconstruction and self.listmode > 0:
                                if not self.loadTOF:
                                    kIndLoc += (self.d_trIndex[0],)
                                    kIndLoc += (self.d_axIndex[0],)
                                else:
                                    kIndLoc += (self.d_trIndex[subset],)
                                    kIndLoc += (self.d_axIndex[subset],)
                            # if self.useTorch:
                            #     kIndLoc += (fD,)
                            # else:
                            if isinstance(f,list):
                                if self.useImages:
                                    chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                    array = cp.cuda.texture.CUDAarray(chl, self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
                                    if self.useTorch:
                                        array.copy_from(fD.reshape((self.Nz[k].item(), self.Ny[k].item(), self.Nx[k].item())))
                                    else:
                                        array.copy_from(f[k].reshape((self.Nz[k].item(), self.Ny[k].item(), self.Nx[k].item())))
                                    res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                    tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                            filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                                    ff = cp.cuda.texture.TextureObject(res, tdes)
                                    kIndLoc += (ff,)
                                else:
                                    if self.useTorch:
                                        kIndLoc += (fD,)
                                    else:
                                        kIndLoc += (f[k],)
                            else:
                                if self.useImages:
                                    chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                    array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item())
                                    if self.useTorch:
                                        array.copy_from(fD.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                                    else:
                                        array.copy_from(f.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                                    res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                    tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                            filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                                    ff = cp.cuda.texture.TextureObject(res, tdes)
                                    kIndLoc += (ff,)
                                else:
                                    if self.useTorch:
                                        kIndLoc += (fD,)
                                    else:
                                        kIndLoc += (f,)
                            if self.useTorch:
                                kIndLoc += (yD,)
                            else:
                                kIndLoc += (y,)
                            kIndLoc += (cp.uint8(self.no_norm),)
                            kIndLoc += (cp.uint64(self.nMeasSubset[subset].item()),)
                            kIndLoc += (cp.uint32(subset),)
                            kIndLoc += (cp.int32(k),)
                        self.knlF((self.globalSizeFP[subset][0] // self.localSizeFP[0], self.globalSizeFP[subset][1] // self.localSizeFP[1], self.globalSizeFP[subset][2]), (self.localSizeFP[0], self.localSizeFP[1], 1),kIndLoc)
                else:
                    import pycuda as cuda
                    if self.useTorch:
                        import torch
                        class Holder(cuda.driver.PointerHolderBase):
                            def __init__(self, t):
                                super(Holder, self).__init__()
                                self.t = t
                                self.gpudata = t.data_ptr()
                            def get_pointer(self):
                                return self.t.data_ptr()
                        if self.subsetType > 7 or self.subsets == 1:
                            y = torch.zeros(self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=torch.float32, device='cuda')
                        else:
                            y = torch.zeros(self.nMeasSubset[subset].item(), dtype=torch.float32, device='cuda')
                        yD = Holder(y)
                    else:
                        if self.subsetType > 7 or self.subsets == 1:
                            y = cuda.gpuarray.zeros(self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=np.float32)
                        else:
                            y = cuda.gpuarray.zeros(self.nMeasSubset[subset].item(), dtype=np.float32)
                    for k in range(self.nMultiVolumes + 1):
                        if isinstance(f,list):
                            if self.use_psf:
                                f[k] = self.computeConvolution(f[k])
                        else:
                            if self.use_psf:
                                f = self.computeConvolution(f)
                        if self.useTorch:
                            fD = Holder(f)
                        kIndLoc = self.kIndF
                        if self.FPType == 1 or self.FPType == 2 or self.FPType == 3 or self.FPType == 4:
                            if (self.attenuation_correction and not self.CTAttenuation):
                                kIndLoc += (self.d_atten[subset].gpudata,)
                        if self.FPType in [1, 2, 3]:
                            if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                                kIndLoc += (np.int64(self.nProjSubset[subset].item()),)
                            if (self.listmode == 0 and not self.CT):
                                kIndLoc += (self.d_x[0].gpudata,)
                            else:
                                kIndLoc += (self.d_x[subset].gpudata, )
                            if (self.CT or self.PET or self.listmode > 0):
                                kIndLoc += (self.d_z[subset].gpudata,)
                            else:
                                kIndLoc += (self.d_z[0].gpudata,)
                            if (self.normalization_correction):
                                kIndLoc += (self.d_norm[subset].gpudata,)
                            elif (self.additionalCorrection):
                                kIndLoc += (self.d_corr[subset].gpudata,)
                            kIndLoc += (self.d_Sens.gpudata,)
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['x'].item()),)
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['y'].item()),)
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_d[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_d[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_d[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['z'].item()),)
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                kIndLoc += (self.d_xyindex[subset].gpudata,)
                                kIndLoc += (self.d_zindex[subset].gpudata,)
                            if self.useTorch:
                                kIndLoc += (fD,)
                            else:
                                if isinstance(f,list):
                                    kIndLoc += (f[k].gpudata,)
                                else:
                                    kIndLoc += (f.gpudata,)
                            if self.useTorch:
                                kIndLoc += (yD,)
                            else:
                                kIndLoc += (y.gpudata,)
                            kIndLoc += (np.uint8(self.no_norm),)
                            kIndLoc += (np.uint64(self.nMeasSubset[subset].item()),)
                            kIndLoc += (np.uint32(subset),)
                            kIndLoc += (np.int32(k),)
                            self.knlF(*kIndLoc, grid=(self.globalSizeFP[subset][0] // self.localSizeFP[0], self.globalSizeFP[subset][1] // self.localSizeFP[1], self.globalSizeFP[subset][2]), block=(self.localSizeFP[0], self.localSizeFP[1], 1))
                if self.useTorch:
                    torch.cuda.synchronize()
                #     if self.useAF:
                #         if isinstance(f,list):
                #             af.device.unlock_array(f[k])
                #         else:
                #             af.device.unlock_array(f)
                # if self.useAF:
                #     af.device.unlock_array(y)
            else:
                import pyopencl as cl
                if not self.loadTOF:
                    if self.useIndexBasedReconstruction and self.listmode > 0:
                        self.d_trIndex[0] = cl.array.to_device(self.queue, self.trIndex[self.nMeas[subset] * 2 : self.nMeas[subset + 1] * 2])
                        self.d_axIndex[0] = cl.array.to_device(self.queue, self.axIndex[self.nMeas[subset] * 2 : self.nMeas[subset + 1] * 2])
                    elif self.listmode > 0:
                        apu = self.x.ravel()
                        self.d_x[0] = cl.array.to_device(self.queue, apu[self.nMeas[subset] * 6 : self.nMeas[subset + 1] * 6])
                if self.useAF:
                    import arrayfire as af
                    if self.subsetType > 7 or self.subsets == 1:
                        y = af.data.constant(0., self.nRowsD * self.nColsD * self.nProjSubset[subset].item())
                    else:
                        y = af.data.constant(0., self.nMeasSubset[subset].item())
                    yPtr = y.raw_ptr()
                    yD = cl.MemoryObject.from_int_ptr(yPtr)
                else:
                    if self.subsetType > 7 or self.subsets == 1:
                        y = cl.array.zeros(self.queue, self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=cl.cltypes.float)
                    else:
                        y = cl.array.zeros(self.queue, self.nMeasSubset[subset].item(), dtype=cl.cltypes.float)
                imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
                mf = cl.mem_flags
                for k in range(self.nMultiVolumes + 1):
                    if self.FPType < 5:
                        d_im = cl.Image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()))
                    else:
                        d_imInt = cl.Image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()))
                        d_im = cl.Image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()))
                    if isinstance(f,list):
                        if self.use_psf:
                            f[k] = self.computeConvolution(f[k])
                        if self.useAF:
                            if self.FPType < 5:
                                fPtr = f[k].raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                                af.device.unlock_array(f[k])
                            else:
                                intIm = af.data.constant(0., self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item())
                                if self.meanFP:
                                    im = af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0)
                                    d_meanFP = af.data.constant(0., self.Nx[k].item() + self.Ny[k].item())
                                    d_meanFP[0:self.Nx[k].item()] = af.flat(af.mean(af.mean(im, dim=0), dim=1))
                                    im -= af.tile(d_meanFP[0:self.Nx[k].item()], d0=im.shape[0], d1=im.shape[1], d3=1)
                                    intIm[1:,1:,:] = af.sat(im)
                                    im.eval()
                                else:
                                    intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0))
                                intIm.eval()
                                intIm = af.flat(intIm)
                                fPtr = intIm.raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_imInt, fD, offset=(0), origin=(0,0,0), region=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()));
                                af.device.unlock_array(intIm)
                                intIm = af.data.constant(0., self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item())
                                intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1))
                                intIm.eval()
                                intIm = af.flat(intIm)
                                fPtr = intIm.raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()));
                                af.device.unlock_array(intIm)
                        else:
                            cl.enqueue_copy(self.queue, d_im, f[k].data, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                    else:
                        if self.use_psf:
                            f = self.computeConvolution(f)
                        if self.useAF:
                            if self.FPType < 5:
                                fPtr = f.raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                                af.device.unlock_array(f)
                            else:
                                intIm = af.data.constant(0., self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item())
                                if self.meanFP:
                                    im = af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0)
                                    d_meanFP = af.data.constant(0., self.Nx[k].item() + self.Ny[k].item())
                                    d_meanFP[0:self.Nx[k].item()] = af.flat(af.mean(af.mean(im, dim=0), dim=1))
                                    im -= af.tile(d_meanFP[0:self.Nx[k].item()], d0=im.shape[0], d1=im.shape[1], d3=1)
                                    intIm[1:,1:,:] = af.sat(im)
                                    af.eval(im)
                                else:
                                    intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0))
                                af.eval(intIm)
                                intIm = af.flat(intIm)
                                fPtr = intIm.raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_imInt, fD, offset=(0), origin=(0,0,0), region=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()));
                                af.device.unlock_array(intIm)
                                intIm = af.data.constant(0., self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item())
                                intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1))
                                af.eval(intIm)
                                intIm = af.flat(intIm)
                                af.sync()
                                fPtr = intIm.raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()));
                                af.device.unlock_array(intIm)
                        else:
                            cl.enqueue_copy(self.queue, d_im, f.data, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                    kIndLoc = self.kIndF
                    if self.FPType == 1 or self.FPType == 2 or self.FPType == 3 or self.FPType == 4:
                        if (self.attenuation_correction and not self.CTAttenuation):
                            self.knlF.set_arg(kIndLoc, self.d_atten[subset].data)
                            kIndLoc += 1
                        # elif self.attenuation_correction and self.CTAttenuation:
                        #     self.knlF.set_arg(kIndLoc, self.d_atten.data)
                        #     kIndLoc += 1
                    if self.FPType == 5 or self.FPType == 4:
                        self.knlF.set_arg(kIndLoc, self.d_Nxyz[k])
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_b[k])
                        kIndLoc += 1
                        if self.FPType == 5:
                            self.knlF.set_arg(kIndLoc, self.dSize[k])
                            kIndLoc += 1
                            self.knlF.set_arg(kIndLoc, self.d_d[k])
                            kIndLoc += 1
                            self.knlF.set_arg(kIndLoc, self.d_Scale[k])
                            kIndLoc += 1
                        else:
                            self.knlF.set_arg(kIndLoc, self.d_bmax[k])
                            kIndLoc += 1
                            self.knlF.set_arg(kIndLoc, self.d_Scale4[k])
                            kIndLoc += 1
                    if self.FPType == 4:
                        self.knlF.set_arg(kIndLoc, d_im)
                        kIndLoc += 1
                        if self.useAF:
                            self.knlF.set_arg(kIndLoc, yD)
                        else:
                            self.knlF.set_arg(kIndLoc, y.data)
                        kIndLoc += 1
                        if (self.listmode == 0 and not self.CT):
                            self.knlF.set_arg(kIndLoc, self.d_x[0].data)
                        else:
                            self.knlF.set_arg(kIndLoc, self.d_x[subset].data)
                        kIndLoc += 1
                        if (self.CT or self.PET or self.listmode > 0):
                            self.knlF.set_arg(kIndLoc, self.d_z[subset].data)
                        else:
                            self.knlF.set_arg(kIndLoc, self.d_z[0].data)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[subset].item()))
                        kIndLoc += 1
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            self.knlF.set_arg(kIndLoc, self.d_xyindex[subset].data)
                            kIndLoc += 1
                            self.knlF.set_arg(kIndLoc, self.d_zindex[subset].data)
                            kIndLoc += 1
                        if (self.normalization_correction):
                            self.knlF.set_arg(kIndLoc, self.d_norm[subset].data)
                            kIndLoc += 1
                        elif (self.additionalCorrection):
                            self.knlF.set_arg(kIndLoc, self.d_corr[subset].data)
                            kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[subset].item()))
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.int)(k))
                    elif self.FPType == 5:
                        self.knlF.set_arg(kIndLoc, self.d_x[subset].data)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_z[subset].data)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, d_im)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, d_imInt)
                        kIndLoc += 1
                        if self.useAF:
                            self.knlF.set_arg(kIndLoc, yD)
                        else:
                            self.knlF.set_arg(kIndLoc, y.data)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[subset].item()))
                        # if self.meanFP:
                    elif self.FPType in [1, 2, 3]:
                        if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                            self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[subset].item()))
                            kIndLoc += 1
                        if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not self.CT) or (not self.loadTOF and self.listmode > 0):
                            self.knlF.set_arg(kIndLoc, self.d_x[0].data)
                        else:
                            self.knlF.set_arg(kIndLoc, self.d_x[subset].data)
                        kIndLoc += 1
                        if (self.CT or self.PET or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                            self.knlF.set_arg(kIndLoc, self.d_z[subset].data)
                        else:
                            self.knlF.set_arg(kIndLoc, self.d_z[0].data)
                        kIndLoc += 1
                        if (self.normalization_correction):
                            self.knlF.set_arg(kIndLoc, self.d_norm[subset].data)
                            kIndLoc += 1
                        elif (self.additionalCorrection):
                            self.knlF.set_arg(kIndLoc, self.d_corr[subset].data)
                            kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_Sens.data)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_Nxyz[k])
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_d[k])
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_b[k])
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_bmax[k])
                        kIndLoc += 1
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            self.knlF.set_arg(kIndLoc, self.d_xyindex[subset].data)
                            kIndLoc += 1
                            self.knlF.set_arg(kIndLoc, self.d_zindex[subset].data)
                            kIndLoc += 1
                        if self.useIndexBasedReconstruction and self.listmode > 0:
                            if not self.loadTOF:
                                self.knlF.set_arg(kIndLoc, self.d_trIndex[0].data)
                                kIndLoc += 1
                                self.knlF.set_arg(kIndLoc, self.d_axIndex[0].data)
                                kIndLoc += 1
                            else:
                                self.knlF.set_arg(kIndLoc, self.d_trIndex[subset].data)
                                kIndLoc += 1
                                self.knlF.set_arg(kIndLoc, self.d_axIndex[subset].data)
                                kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, d_im)
                        kIndLoc += 1
                        if self.useAF:
                            self.knlF.set_arg(kIndLoc, yD)
                        else:
                            self.knlF.set_arg(kIndLoc, y.data)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[subset].item()))
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.int)(k))
                    cl.enqueue_nd_range_kernel(self.queue, self.knlF, self.globalSizeFP[subset], self.localSizeFP)
                    self.queue.finish()
            if self.useAF:
                af.device.unlock_array(y)
        return y

    def backwardProject(self, y, subset = -1):
        if subset == -1:
            subset = self.subset
        if self.nMultiVolumes > 0:
            f = [None] * (self.nMultiVolumes + 1)
        if self.projector_type == 6:
            if not self.useCUDA:
                import arrayfire as af
                fProj = af.data.moddims(y, self.nRowsD, d1=self.nColsD, d2=self.nProjSubset[subset].item())
                for ii in range(self.nMultiVolumes + 1):
                    u1 = np.sum(self.nProjSubset[:subset])
                    if self.nMultiVolumes > 0:
                        f[ii] = af.data.constant(0, self.Nx[ii].item() * self.Ny[ii].item() * self.Nz[ii].item(), d1 = self.nProjSubset[subset].item())
                    else:
                        f = af.data.constant(0, self.Nx[ii].item() * self.Ny[ii].item() * self.Nz[ii].item(), d1 = self.nProjSubset[subset].item())
                    for kk in range(self.nProjSubset[subset].item()):
                        fApu = af.data.constant(0, self.Nx[ii].item(), d1=self.Ny[ii].item(), d2=self.Nz[ii].item()) # [128, 128, 96]
                        kuvaRot = fProj[:,:,kk] # [128, 96]
                        kuvaRot = af.data.reorder(kuvaRot, 1, 0, 2) # [96, 128]
                        kuvaRot = af.signal.convolve2(kuvaRot, self.d_gFilter[:, :, :, u1]) # [96, 128, 128]
                        kuvaRot = kuvaRot[:,:, self.blurPlanes[u1].item():] # [96, 128, 90]
                        kuvaRot = af.data.reorder(kuvaRot, 2, 1, 0) # [90, 128, 96]
                        fApu[self.blurPlanes[u1].item():,:,:] = kuvaRot
                        fApu = af.image.rotate(fApu, self.angles[u1], method=af.INTERP.BILINEAR)
                        if isinstance(f, list):
                            f[ii][:, kk] = af.flat(fApu)
                        else:
                            f[:, kk] = af.flat(fApu)
                        u1 += 1
                    if isinstance(f, list):
                        f[ii] = af.sum(f[ii], 1)
                    else:
                        f = af.sum(f, 1)
            else:
                if not self.useCuPy:
                    raise ValueError('SPECT is only supported on CUDA with CuPy!')
                else:
                    import cupy as cp
                    if self.useTorch:
                        import torch
                        from torchvision.transforms.functional import rotate
                        from torchvision.transforms import InterpolationMode
                        import torch.nn.functional as F
                        fProj = torch.reshape(y, (self.nProjSubset[subset].item(), self.nColsD, self.nRowsD))
                    for ii in range(self.nMultiVolumes + 1):
                        u1 = np.sum(self.nProjSubset[:subset])
                        if self.nMultiVolumes > 0:
                            f[ii] = torch.zeros((self.nProjSubset[subset].item(), self.Nx[ii].item() * self.Ny[ii].item() * self.Nz[ii].item()), dtype=torch.float32).cuda()
                        else:
                            f = torch.zeros((self.nProjSubset[subset].item(), self.Nx[ii].item() * self.Ny[ii].item() * self.Nz[ii].item()), dtype=torch.float32).cuda()
                        for kk in range(self.nProjSubset[subset].item()):
                            fApu = torch.zeros((self.Nz[ii].item(), self.Ny[ii].item(), self.Nx[ii].item()), dtype=torch.float32).cuda()
                            kuvaRot = fProj[kk,:,:]
                            kuvaRot = torch.permute(kuvaRot, (1, 0)).unsqueeze(0).unsqueeze(0)
                            kuvaRot = F.conv2d(kuvaRot, self.d_gFilter[:, :, :, :, u1], padding=(self.d_gFilter.shape[2] // 2, self.d_gFilter.shape[3] // 2))
                            kuvaRot = kuvaRot.squeeze(0)
                            kuvaRot = kuvaRot[self.blurPlanes[u1].item():, :,:]
                            kuvaRot = torch.permute(kuvaRot, (2, 1, 0))
                            # kuvaRot = cp.transpose(kuvaRot, (2, 1, 0))
                            fApu[:,:,self.blurPlanes[u1].item():] = kuvaRot
                            fApu = rotate(fApu, self.angles[u1].item(), InterpolationMode.BILINEAR)
                            if isinstance(f, list):
                                f[ii][kk,:] = torch.ravel(fApu)
                            else:
                                f[kk,:] = torch.ravel(fApu)
                            u1 += 1
                        if isinstance(f, list):
                            f[ii] = torch.sum(f[ii], 0),
                        else:
                            f = torch.sum(f, 0)
        else:
            if self.useCUDA:
                if self.useCuPy:
                    import cupy as cp
                    for k in range(self.nMultiVolumes + 1):
                        if self.useTorch:
                            import torch
                            if self.nMultiVolumes > 0:
                                f[k] = torch.zeros(self.N[k].item(), dtype=torch.float32, device='cuda')
                                fD = cp.asarray(f[k])
                            else:
                                f = torch.zeros(self.N[k].item(), dtype=torch.float32, device='cuda')
                                fD = cp.asarray(f)
                            yD = cp.asarray(y)
                        else:
                            if self.nMultiVolumes > 0:
                                f[k] = cp.zeros(self.N[k].item(), dtype=cp.float32)
                            else:
                                f = cp.zeros(self.N[k].item(), dtype=cp.float32)
                        if self.BPType == 5:
                            # y1 = cp.asarray(np.load('testi.npy'))
                            yy = cp.zeros((self.nRowsD+1,self.nColsD+1,self.nProjSubset[subset].item()), dtype=cp.float32, order='F')
                            if self.useTorch:
                                yy[1:,1:,:] = yD.reshape((self.nRowsD,self.nColsD,self.nProjSubset[subset].item()), order='F')
                            else:
                                yy[1:,1:,:] = y.reshape((self.nRowsD,self.nColsD,self.nProjSubset[subset].item()), order='F')
                            yy = yy.cumsum(0)
                            yy = yy.cumsum(1)
                            yy = yy.ravel(order='F')
                        kIndLoc = self.kIndB
                        if self.BPType in [1, 2, 3]:
                            if (self.attenuation_correction and not self.CTAttenuation):
                                kIndLoc += (self.d_atten[subset],)
                            if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                                kIndLoc += ((self.nProjSubset[subset].item()),)
                            if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not self.CT) or (not self.loadTOF and self.listmode > 0):
                                kIndLoc += (self.d_x[0],)
                            else:
                                kIndLoc += (self.d_x[subset],)
                            if (self.CT or self.PET or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                                kIndLoc += (self.d_z[subset],)
                            else:
                                kIndLoc += (self.d_z[0],)
                            if (self.normalization_correction):
                                kIndLoc += (self.d_norm[subset],)
                            if (self.additionalCorrection):
                                kIndLoc += (self.d_corr[subset],)
                            kIndLoc += (self.d_Sens,)
                            kIndLoc += (cp.uint32(self.Nx[k].item()),)
                            kIndLoc += (cp.uint32(self.Ny[k].item()),)
                            kIndLoc += (cp.uint32(self.Nz[k].item()),)
                            kIndLoc += (cp.float32(self.dx[k].item()),)
                            kIndLoc += (cp.float32(self.dy[k].item()),)
                            kIndLoc += (cp.float32(self.dz[k].item()),)
                            kIndLoc += (cp.float32(self.bx[k].item()),)
                            kIndLoc += (cp.float32(self.by[k].item()),)
                            kIndLoc += (cp.float32(self.bz[k].item()),)
                            kIndLoc += (cp.float32(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item()),)
                            kIndLoc += (cp.float32(self.by[k].item() + self.Ny[k].item() * self.dy[k].item()),)
                            kIndLoc += (cp.float32(self.bz[k].item() + self.Nz[k].item() * self.dz[k].item()),)
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                kIndLoc += (self.d_xyindex[subset],)
                                kIndLoc += (self.d_zindex[subset],)
                            if self.useIndexBasedReconstruction and self.listmode > 0:
                                if not self.loadTOF:
                                    kIndLoc += (self.d_trIndex[0],)
                                    kIndLoc += (self.d_axIndex[0],)
                                else:
                                    kIndLoc += (self.d_trIndex[subset],)
                                    kIndLoc += (self.d_axIndex[subset],)
                            if self.useTorch:
                                kIndLoc += (yD,)
                            else:
                                kIndLoc += (y,)
                            if self.useTorch:
                                kIndLoc += (fD,)
                            else:
                                if self.nMultiVolumes > 0:
                                    kIndLoc += (f[k],)
                                else:
                                    kIndLoc += (f,)
                            kIndLoc += (cp.uint8(self.no_norm),)
                            kIndLoc += (cp.uint64(self.nMeasSubset[subset].item()),)
                            kIndLoc += (cp.uint32(subset),)
                            kIndLoc += (cp.int32(k),)
                        else:
                            if self.CT:
                                if self.OffsetLimit.size > 0:
                                    kIndLoc += (self.d_T[subset],)
                                if self.BPType == 5 or self.BPType == 4:
                                    kIndLoc += (cp.uint32(self.Nx[k].item()),)
                                    kIndLoc += (cp.uint32(self.Ny[k].item()),)
                                    kIndLoc += (cp.uint32(self.Nz[k].item()),)
                                    kIndLoc += (cp.float32(self.bx[k].item()),)
                                    kIndLoc += (cp.float32(self.by[k].item()),)
                                    kIndLoc += (cp.float32(self.bz[k].item()),)
                                    kIndLoc += (cp.float32(self.dx[k].item()),)
                                    kIndLoc += (cp.float32(self.dy[k].item()),)
                                    kIndLoc += (cp.float32(self.dz[k].item()),)
                                    if self.BPType == 5:
                                        kIndLoc += (cp.float32(self.dScaleX[k].item()),)
                                        kIndLoc += (cp.float32(self.dScaleY[k].item()),)
                                        kIndLoc += (cp.float32(self.dScaleZ[k].item()),)
                                        kIndLoc += (cp.float32(self.dSizeXBP),)
                                        kIndLoc += (cp.float32(self.dSizeZBP),)
                                    else:
                                        kIndLoc += (cp.float32(self.kerroin[k].item()),)
                                if self.BPType == 4:
                                    # if self.useTorch:
                                    #     kIndLoc += (yD,)
                                    #     kIndLoc += (fD,)
                                    # else:
                                    if self.useImages:
                                        chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                        array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD, self.nProjSubset[subset].item())
                                        if self.useTorch:
                                            array.copy_from(yD.reshape((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD)))
                                        else:
                                            array.copy_from(y.reshape((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD)))
                                        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                                filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=1)
                                        yy = cp.cuda.texture.TextureObject(res, tdes)
                                        kIndLoc += (yy,)
                                    else:
                                        if self.useTorch:
                                            kIndLoc += (yD,)
                                        else:
                                            kIndLoc += (y,)
                                    if self.useTorch:
                                        kIndLoc += (fD,)
                                    else:
                                        if isinstance(f, list):
                                            kIndLoc += (f[k],)
                                        else:
                                            kIndLoc += (f,)
                                    kIndLoc += (self.d_x[subset],)
                                    kIndLoc += (self.d_z[subset],)
                                    kIndLoc += (self.d_Sens,)
                                else:
                                    kIndLoc += (self.d_x[subset],)
                                    kIndLoc += (self.d_z[subset],)
                                    # if self.useTorch:
                                    #     kIndLoc += (yD,)
                                    #     kIndLoc += (fD,)
                                    # else:
                                    if self.useImages:
                                        chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                        array = cp.cuda.texture.CUDAarray(chl, self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[subset].item())
                                        array.copy_from(yy.reshape((self.nProjSubset[subset].item(), self.nColsD + 1, self.nRowsD + 1)))
                                        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                                filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                                        yy = cp.cuda.texture.TextureObject(res, tdes)
                                        kIndLoc += (yy,)
                                    else:
                                        if self.useTorch:
                                            kIndLoc += (yD,)
                                        else:
                                            kIndLoc += (y,)
                                    if self.useTorch:
                                        kIndLoc += (fD,)
                                    else:
                                        if isinstance(f, list):
                                            kIndLoc += (f[k],)
                                        else:
                                            kIndLoc += (f,)
                                    kIndLoc += (self.d_Sens,)
                                    # if self.meanBP:
                                    #     kIndLoc += (dMeanBP)
                            else:
                                kIndLoc += (cp.uint32(self.Nx[k].item()),)
                                kIndLoc += (cp.uint32(self.Ny[k].item()),)
                                kIndLoc += (cp.uint32(self.Nz[k].item()),)
                                kIndLoc += (cp.float32(self.bx[k].item()),)
                                kIndLoc += (cp.float32(self.by[k].item()),)
                                kIndLoc += (cp.float32(self.bz[k].item()),)
                                kIndLoc += (cp.float32(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item()),)
                                kIndLoc += (cp.float32(self.by[k].item() + self.Ny[k].item() * self.dy[k].item()),)
                                kIndLoc += (cp.float32(self.bz[k].item() + self.Nz[k].item() * self.dz[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleX4[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleY4[k].item()),)
                                kIndLoc += (cp.float32(self.dScaleZ4[k].item()),)
                                # if self.useTorch:
                                #     kIndLoc += (yD,)
                                #     kIndLoc += (fD,)
                                # else:
                                if self.useImages:
                                    chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                    array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD, self.nProjSubset[subset].item())
                                    if self.useTorch:
                                        array.copy_from(yD.reshape((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD)))
                                    else:
                                        array.copy_from(y.reshape((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD)))
                                    res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                    tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                            filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=1)
                                    yy = cp.cuda.texture.TextureObject(res, tdes)
                                    kIndLoc += (yy,)
                                else:
                                    if self.useTorch:
                                        kIndLoc += (yD,)
                                    else:
                                        kIndLoc += (y,)
                                if self.useTorch:
                                    kIndLoc += (fD,)
                                else:
                                    if isinstance(f, list):
                                        kIndLoc += (f[k],)
                                    else:
                                        kIndLoc += (f,)
                                if self.listmode == 0 and not self.CT:
                                    kIndLoc += (self.d_x[0],)
                                else:
                                    kIndLoc += (self.d_x[subset],)
                                if (self.CT or self.PET or self.listmode > 0):
                                    kIndLoc += (self.d_z[subset],)
                                else:
                                    kIndLoc += (self.d_z[0],)
                                kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                                if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                    kIndLoc += (self.d_xyindex[subset],)
                                    kIndLoc += (self.d_zindex[subset],)
                                if (self.normalization_correction):
                                    kIndLoc += (self.d_norm[subset],)
                                elif (self.additionalCorrection):
                                    kIndLoc += (self.d_corr[subset],)
                                kIndLoc += (self.d_Sens,)
                            kIndLoc += (cp.uint8(self.no_norm),)
                            if self.CT:
                                kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                            else:
                                kIndLoc += (cp.uint64(self.nMeasSubset[subset].item()),)
                                kIndLoc += (cp.uint32(subset),)
                            kIndLoc += (cp.int32(k),)
                        self.knlB((self.globalSizeBP[subset][k][0] // self.localSizeBP[0], self.globalSizeBP[subset][k][1] // self.localSizeBP[1], self.globalSizeBP[subset][k][2]), (self.localSizeBP[0], self.localSizeBP[1], 1), kIndLoc)
                else:
                    import pycuda as cuda
                    if self.useTorch:
                        import torch
                        class Holder(cuda.driver.PointerHolderBase):
                            def __init__(self, t):
                                super(Holder, self).__init__()
                                self.t = t
                                self.gpudata = t.data_ptr()
                            def get_pointer(self):
                                return self.t.data_ptr()
                        yD = Holder(y)
                    for k in range(self.nMultiVolumes + 1):
                        if self.useTorch:
                            if self.nMultiVolumes > 0:
                                f[k] = torch.zeros(self.N[k].item(), dtype=torch.float32, device='cuda')
                                fD = Holder(f[k])
                            else:
                                f = torch.zeros(self.N[k].item(), dtype=torch.float32, device='cuda')
                                fD = Holder(f)
                        else:
                            if self.nMultiVolumes > 0:
                                f[k] = cuda.gpuarray.zeros(self.N[k].item(), dtype=np.float32)
                            else:
                                f = cuda.gpuarray.zeros(self.N[k].item(), dtype=np.float32)
                        kIndLoc = self.kIndB
                        if self.BPType in [1, 2, 3]:
                            if (self.attenuation_correction and not self.CTAttenuation):
                                kIndLoc += (self.d_atten[subset].gpudata,)
                            if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                                kIndLoc += ((self.nProjSubset[subset].item()),)
                            if (self.listmode == 0 and not self.CT):
                                kIndLoc += (self.d_x[0].gpudata,)
                            else:
                                kIndLoc += (self.d_x[subset].gpudata,)
                            if (self.CT or self.PET or self.listmode > 0):
                                kIndLoc += (self.d_z[subset].gpudata,)
                            else:
                                kIndLoc += (self.d_z[0].gpudata,)
                            if (self.normalization_correction):
                                kIndLoc += (self.d_norm[subset].gpudata,)
                            if (self.additionalCorrection):
                                kIndLoc += (self.d_corr[subset].gpudata,)
                            kIndLoc += (self.d_Sens.gpudata,)
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['x'].item()),)
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['y'].item()),)
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_d[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_d[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_d[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['z'].item()),)
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                kIndLoc += (self.d_xyindex[subset].gpudata,)
                                kIndLoc += (self.d_zindex[subset].gpudata,)
                            if self.useTorch:
                                kIndLoc += (yD,)
                                kIndLoc += (fD,)
                            else:
                                kIndLoc += (y.gpudata,)
                                if self.nMultiVolumes > 0:
                                    kIndLoc += (f[k].gpudata,)
                                else:
                                    kIndLoc += (f.gpudata,)
                            kIndLoc += (np.uint8(self.no_norm),)
                            kIndLoc += (np.uint64(self.nMeasSubset[subset].item()),)
                            kIndLoc += (np.uint32(subset),)
                            kIndLoc += (np.int32(k),)
                        else:
                            if self.CT:
                                if self.OffsetLimit.size > 0:
                                    kIndLoc += (self.d_T[subset].gpudata,)
                                if self.BPType == 5 or self.BPType == 4:
                                    kIndLoc += (np.uint32(self.d_Nxyz[k]['x'].item()),)
                                    kIndLoc += (np.uint32(self.d_Nxyz[k]['y'].item()),)
                                    kIndLoc += (np.uint32(self.d_Nxyz[k]['z'].item()),)
                                    kIndLoc += (np.float32(self.d_b[k]['x'].item()),)
                                    kIndLoc += (np.float32(self.d_b[k]['y'].item()),)
                                    kIndLoc += (np.float32(self.d_b[k]['z'].item()),)
                                    kIndLoc += (np.float32(self.d_d[k]['x'].item()),)
                                    kIndLoc += (np.float32(self.d_d[k]['y'].item()),)
                                    kIndLoc += (np.float32(self.d_d[k]['z'].item()),)
                                    if self.BPType == 5:
                                        kIndLoc += (np.float32(self.dScaleX[k].item()),)
                                        kIndLoc += (np.float32(self.dScaleY[k].item()),)
                                        kIndLoc += (np.float32(self.dScaleZ[k].item()),)
                                        kIndLoc += (np.float32(self.dSizeXBP),)
                                        kIndLoc += (np.float32(self.dSizeZBP),)
                                    else:
                                        kIndLoc += (np.float32(self.kerroin[k].item()),)
                                if self.BPType == 4:
                                    if self.useTorch:
                                        kIndLoc += (yD,)
                                        kIndLoc += (fD,)
                                    else:
                                        kIndLoc += (y.gpudata,)
                                        if isinstance(f, list):
                                            kIndLoc += (f[k].gpudata,)
                                        else:
                                            kIndLoc += (f.gpudata,)
                                    kIndLoc += (self.d_x[subset].gpudata,)
                                    kIndLoc += (self.d_z[subset].gpudata,)
                                    kIndLoc += (self.d_Sens.gpudata,)
                                else:
                                    kIndLoc += (self.d_x[subset].gpudata,)
                                    kIndLoc += (self.d_z[subset].gpudata,)
                                    if self.useTorch:
                                        kIndLoc += (yD,)
                                        kIndLoc += (fD,)
                                    else:
                                        kIndLoc += (y.gpudata,)
                                        if isinstance(f, list):
                                            kIndLoc += (f[k].gpudata,)
                                        else:
                                            kIndLoc += (f.gpudata,)
                                    kIndLoc += (self.d_Sens.gpudata,)
                                    # if self.meanBP:
                                    #     kIndLoc += (dMeanBP)
                            else:
                                kIndLoc += (np.uint32(self.d_Nxyz[k]['x'].item()),)
                                kIndLoc += (np.uint32(self.d_Nxyz[k]['y'].item()),)
                                kIndLoc += (np.uint32(self.d_Nxyz[k]['z'].item()),)
                                kIndLoc += (np.float32(self.d_b[k]['x'].item()),)
                                kIndLoc += (np.float32(self.d_b[k]['y'].item()),)
                                kIndLoc += (np.float32(self.d_b[k]['z'].item()),)
                                kIndLoc += (np.float32(self.d_bmax[k]['x'].item()),)
                                kIndLoc += (np.float32(self.d_bmax[k]['y'].item()),)
                                kIndLoc += (np.float32(self.d_bmax[k]['z'].item()),)
                                kIndLoc += (np.float32(self.d_Scale4[k]['x'].item()),)
                                kIndLoc += (np.float32(self.d_Scale4[k]['y'].item()),)
                                kIndLoc += (np.float32(self.d_Scale4[k]['z'].item()),)
                                if self.useTorch:
                                    kIndLoc += (yD,)
                                    kIndLoc += (fD,)
                                else:
                                    kIndLoc += (y.gpudata,)
                                    if isinstance(f, list):
                                        kIndLoc += (f[k].gpudata,)
                                    else:
                                        kIndLoc += (f.gpudata,)
                                if self.listmode == 0 and not self.CT:
                                    kIndLoc += (self.d_x[0].gpudata,)
                                else:
                                    kIndLoc += (self.d_x[subset].gpudata,)
                                if (self.CT or self.PET or self.listmode > 0):
                                    kIndLoc += (self.d_z[subset].gpudata,)
                                else:
                                    kIndLoc += (self.d_z[0].gpudata,)
                                kIndLoc += (np.int64(self.nProjSubset[subset].item()),)
                                if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                    kIndLoc += (self.d_xyindex[subset].gpudata,)
                                    kIndLoc += (self.d_zindex[subset].gpudata,)
                                if (self.normalization_correction):
                                    kIndLoc += (self.d_norm[subset].gpudata,)
                                elif (self.additionalCorrection):
                                    kIndLoc += (self.d_corr[subset].gpudata,)
                                kIndLoc += (self.d_Sens.gpudata,)
                            kIndLoc += (np.uint8(self.no_norm),)
                            if self.CT:
                                kIndLoc += (np.int64(self.nProjSubset[subset].item()),)
                            else:
                                kIndLoc += (np.uint64(self.nMeasSubset[subset].item()),)
                                kIndLoc += (np.uint32(subset),)
                            kIndLoc += (np.int32(k),)
                        self.knlB(*kIndLoc, grid=(self.globalSizeBP[subset][k][0] // self.localSizeBP[0], self.globalSizeBP[subset][k][1] // self.localSizeBP[1], self.globalSizeBP[subset][k][2]), block=(self.localSizeBP[0], self.localSizeBP[1], 1))
                    # if self.useAF:
                    #     if self.nMultiVolumes > 0:
                    #         af.device.unlock_array(f[k])
                    #     else:
                    #         af.device.unlock_array(f)
                    #     af.device.unlock_array(y)
                    
                    if self.use_psf:
                        if self.nMultiVolumes > 0:
                            f[k] = self.computeConvolution(f[k])
                        else:
                            f = self.computeConvolution(f)
                if self.useTorch:
                    torch.cuda.synchronize()
            else:
                import pyopencl as cl
                if self.useAF:
                    import arrayfire as af
                    cltype = af.Dtype.f32
                    if self.use_64bit_atomics:
                        cltype = af.Dtype.u64
                    elif self.use_32bit_atomics:
                        cltype = af.Dtype.u32
                    yPtr = y.raw_ptr()
                    yD = cl.MemoryObject.from_int_ptr(yPtr)
                else:
                    cltype = cl.cltypes.float
                    if self.use_64bit_atomics:
                        cltype = cl.cltypes.ulong
                    elif self.use_32bit_atomics:
                        cltype = cl.cltypes.uint
                if self.CT and self.BPType in [4,5]:
                    imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
                    if self.BPType < 5:
                        d_im = cl.Image(self.clctx, cl.mem_flags.READ_ONLY, imformat, shape=(self.nRowsD, self.nColsD, self.nProjSubset[subset].item()))
                        if self.useAF:
                            cl.enqueue_copy(self.queue, d_im, yD, offset=(0), origin=(0,0,0), region=(self.nRowsD, self.nColsD, self.nProjSubset[subset].item()));
                        else:
                            cl.enqueue_copy(self.queue, d_im, y.data, offset=(0), origin=(0,0,0), region=(self.nRowsD, self.nColsD, self.nProjSubset[subset].item()));
                    else:
                        d_im = cl.Image(self.clctx, cl.mem_flags.READ_ONLY, imformat, shape=(self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[subset].item()))
                        y = af.moddims(y, self.nRowsD, d1=self.nColsD, d2=self.nProjSubset[subset].item())
                        if self.meanBP:
                            d_meanBP = af.mean(af.mean(y, dim=0), dim=1)
                            y -= af.tile(d_meanBP, self.nRowsD, d1=self.nColsD)
                            mPtr = d_meanBP.raw_ptr()
                            dMeanBP = cl.MemoryObject.from_int_ptr(mPtr)
                        y = af.sat(y)
                        y = af.join(0, af.data.constant(0, 1, d1=y.shape[1], d2=y.shape[2]), y)
                        y = af.flat(af.join(1, af.data.constant(0, y.shape[0], 1, y.shape[2]), y))
                        yPtr = y.raw_ptr()
                        yD = cl.MemoryObject.from_int_ptr(yPtr)
                        cl.enqueue_copy(self.queue, d_im, yD, offset=(0), origin=(0,0,0), region=(self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[subset].item()));
                
                for k in range(self.nMultiVolumes + 1):
                    if self.useAF:
                        if self.nMultiVolumes > 0:
                            f[k] = af.data.constant(0, self.N[k].item(), dtype=cltype)
                            fPtr = f[k].raw_ptr()
                        else:
                            f = af.data.constant(0, self.N[k].item(), dtype=cltype)
                            fPtr = f.raw_ptr()
                        fD = cl.MemoryObject.from_int_ptr(fPtr)
                    else:
                        if self.nMultiVolumes > 0:
                            f[k] = cl.array.zeros(self.queue, self.N[k].item(), dtype=cltype)
                        else:
                            f = cl.array.zeros(self.queue, self.N[k].item(), dtype=cltype)
                    kIndLoc = self.kIndB
                    if self.BPType in [1, 2, 3]:
                        if (self.attenuation_correction and not self.CTAttenuation):
                            self.knlB.set_arg(kIndLoc, self.d_atten[subset].data)
                            kIndLoc += 1
                        if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                            self.knlB.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[subset].item()))
                            kIndLoc += 1
                        if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not self.CT) or (not self.loadTOF and self.listmode > 0):
                            self.knlB.set_arg(kIndLoc, self.d_x[0].data)
                        else:
                            self.knlB.set_arg(kIndLoc, self.d_x[subset].data)
                        kIndLoc += 1
                        if (self.CT or self.PET or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                            self.knlB.set_arg(kIndLoc, self.d_z[subset].data)
                        else:
                            self.knlB.set_arg(kIndLoc, self.d_z[0].data)
                        kIndLoc += 1
                        if (self.normalization_correction):
                            self.knlB.set_arg(kIndLoc, self.d_norm[subset].data)
                            kIndLoc += 1
                        if (self.additionalCorrection):
                            self.knlB.set_arg(kIndLoc, self.d_corr[subset].data)
                            kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_Nxyz[k])
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_d[k])
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_b[k])
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_bmax[k])
                        kIndLoc += 1
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            self.knlB.set_arg(kIndLoc, self.d_xyindex[subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_zindex[subset].data)
                            kIndLoc += 1
                        if self.useIndexBasedReconstruction and self.listmode > 0:
                            if not self.loadTOF:
                                self.knlB.set_arg(kIndLoc, self.d_trIndex[0].data)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_axIndex[0].data)
                                kIndLoc += 1
                            else:
                                self.knlB.set_arg(kIndLoc, self.d_trIndex[subset].data)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_axIndex[subset].data)
                                kIndLoc += 1
                        if self.useAF:
                            self.knlB.set_arg(kIndLoc, yD)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, fD)
                        else:
                            self.knlB.set_arg(kIndLoc, y.data)
                            kIndLoc += 1
                            if self.nMultiVolumes > 0:
                                self.knlB.set_arg(kIndLoc, f[k].data)
                            else:
                                self.knlB.set_arg(kIndLoc, f.data)
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[subset].item()))
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.int)(k))
                    else:
                        if self.CT:
                            if self.OffsetLimit.size > 0:
                                self.knlB.set_arg(kIndLoc, self.d_T[subset].data)
                                kIndLoc += 1
                            if self.BPType == 5 or self.BPType == 4:
                                self.knlB.set_arg(kIndLoc, self.d_Nxyz[k])
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_b[k])
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_d[k])
                                kIndLoc += 1
                                if self.BPType == 5:
                                    self.knlB.set_arg(kIndLoc, self.d_Scale[k])
                                    kIndLoc += 1
                                    self.knlB.set_arg(kIndLoc, self.dSizeBP)
                                    kIndLoc += 1
                                else:
                                    self.knlB.set_arg(kIndLoc, (cl.cltypes.float)(self.kerroin[k].item()))
                                    kIndLoc += 1
                            if self.BPType == 4:
                                self.knlB.set_arg(kIndLoc, d_im)
                                kIndLoc += 1
                                if self.useAF:
                                    self.knlB.set_arg(kIndLoc, fD)
                                else:
                                    if isinstance(f, list):
                                        self.knlB.set_arg(kIndLoc, f[k].data)
                                    else:
                                        self.knlB.set_arg(kIndLoc, f.data)
                                kIndLoc += 1
                                if not self.loadTOF and self.listmode > 0:
                                    self.knlB.set_arg(kIndLoc, self.d_x[0].data)
                                else:
                                    self.knlB.set_arg(kIndLoc, self.d_x[subset].data)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_z[subset].data)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                                kIndLoc += 1
                            else:
                                if not self.loadTOF and self.listmode > 0:
                                    self.knlB.set_arg(kIndLoc, self.d_x[0].data)
                                else:
                                    self.knlB.set_arg(kIndLoc, self.d_x[subset].data)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_z[subset].data)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, d_im)
                                kIndLoc += 1
                                if self.useAF:
                                    self.knlB.set_arg(kIndLoc, fD)
                                else:
                                    if isinstance(f, list):
                                        self.knlB.set_arg(kIndLoc, f[k].data)
                                    else:
                                        self.knlB.set_arg(kIndLoc, f.data)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                                kIndLoc += 1
                                if self.meanBP:
                                    self.knlB.set_arg(kIndLoc, dMeanBP)
                                    kIndLoc += 1
                        else:
                            self.knlB.set_arg(kIndLoc, self.d_Nxyz[k])
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_b[k])
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_bmax[k])
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_Scale4[k])
                            kIndLoc += 1
                            if self.useAF:
                                self.knlB.set_arg(kIndLoc, yD)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, fD)
                            else:
                                self.knlB.set_arg(kIndLoc, y.data)
                                kIndLoc += 1
                                if isinstance(f, list):
                                    self.knlB.set_arg(kIndLoc, f[k].data)
                                else:
                                    self.knlB.set_arg(kIndLoc, f.data)
                            kIndLoc += 1
                            if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not self.CT) or (not self.loadTOF and self.listmode > 0):
                                self.knlB.set_arg(kIndLoc, self.d_x[0].data)
                            else:
                                self.knlB.set_arg(kIndLoc, self.d_x[subset].data)
                            kIndLoc += 1
                            if (self.CT or self.PET or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                                self.knlB.set_arg(kIndLoc, self.d_z[subset].data)
                            else:
                                self.knlB.set_arg(kIndLoc, self.d_z[0].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nProjSubset[subset].item()))
                            kIndLoc += 1
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                self.knlB.set_arg(kIndLoc, self.d_xyindex[subset].data)
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.d_zindex[subset].data)
                                kIndLoc += 1
                            if (self.normalization_correction):
                                self.knlB.set_arg(kIndLoc, self.d_norm[subset].data)
                                kIndLoc += 1
                            elif (self.additionalCorrection):
                                self.knlB.set_arg(kIndLoc, self.d_corr[subset].data)
                                kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                            kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                        kIndLoc += 1
                        if self.CT:
                            self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nProjSubset[subset].item()))
                            kIndLoc += 1
                        else:
                            self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[subset].item()))
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                            kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.int)(k))
                                
                    cl.enqueue_nd_range_kernel(self.queue, self.knlB, self.globalSizeBP[subset][k], self.localSizeBP)
                    self.queue.finish()
                    if self.useAF:
                        if self.nMultiVolumes > 0:
                            af.device.unlock_array(f[k])
                        else:
                            af.device.unlock_array(f)
                        af.device.unlock_array(y)
                        if self.use_64bit_atomics:
                            if self.nMultiVolumes > 0:
                                f[k] = f[k].as_type(af.Dtype.f32) / self.TH
                            else:
                                f = f.as_type(af.Dtype.f32) / self.TH
                        elif self.use_32bit_atomics:
                            if self.nMultiVolumes > 0:
                                f[k] = f[k].as_type(af.Dtype.f32) / self.TH32
                            else:
                                f = f.as_type(af.Dtype.f32) / self.TH32
                    else:
                        if self.use_64bit_atomics:
                            if self.nMultiVolumes > 0:
                                f[k] = f[k].astype(cl.cltypes.float) / self.TH
                            else:
                                f = f.astype(cl.cltypes.float) / self.TH
                        elif self.use_32bit_atomics:
                            if self.nMultiVolumes > 0:
                                f[k] = f[k].astype(cl.cltypes.float) / self.TH32
                            else:
                                f = f.astype(cl.cltypes.float) / self.TH32
                    if self.use_psf:
                        if self.nMultiVolumes > 0:
                            f[k] = self.computeConvolution(f[k])
                        else:
                            f = self.computeConvolution(f)
        return f
    
    def T(self):
        self.trans = True
        return self
        
    
    def __imul__(self, B):
            return self.__mul__(self, B)
    
    def __mul__(self, B):
        if self.trans:
            self.trans = False
            return self.backwardProject(B, self.subset)
        else:
            return self.forwardProject(B, self.subset)
        
    def __rmul__(self, B):
        return self
    
    def __truediv__(self, B):
        return self
        
    class parameters(ctypes.Structure):
        _pack_  = 1
        _fields_ = [
            ('use_raw_data', ctypes.c_uint8),
            ('listmode', ctypes.c_uint8),
            ('verbose', ctypes.c_int8),
            ('n_rays_transaxial', ctypes.c_uint16),
            ('n_rays_axial', ctypes.c_uint16),
            ('projector_type', ctypes.c_uint32),
            ('attenuation_correction', ctypes.c_uint32),
            ('additionalCorrection', ctypes.c_uint32),
            ('normalization_correction', ctypes.c_uint32),
            ('randoms_correction', ctypes.c_uint32),
            ('nColsD', ctypes.c_uint32),
            ('nRowsD', ctypes.c_uint32),
            ('Nang', ctypes.c_uint32),
            ('Ndist', ctypes.c_uint32),
            ('subsets', ctypes.c_uint32),
            ('det_per_ring', ctypes.c_uint32),
            ('rings', ctypes.c_uint32),
            ('NxOrig', ctypes.c_uint32),
            ('NyOrig', ctypes.c_uint32),
            ('NzOrig', ctypes.c_uint32),
            ('NxPrior', ctypes.c_uint32),
            ('NyPrior', ctypes.c_uint32),
            ('NzPrior', ctypes.c_uint32),
            ('Niter', ctypes.c_uint32),
            ('Nt', ctypes.c_uint32),
            ('subsetType', ctypes.c_uint32),
            ('nMultiVolumes', ctypes.c_uint32),
            ('nLayers', ctypes.c_uint32),
            ('PDAdaptiveType', ctypes.c_uint32),
            ('powerIterations', ctypes.c_uint32),
            ('deblur_iterations', ctypes.c_uint32),
            ('gradInitIter', ctypes.c_uint32),
            ('gradLastIter', ctypes.c_uint32),
            ('filteringIterations', ctypes.c_uint32),
            ('mean_type', ctypes.c_uint32),
            ('Ndx', ctypes.c_uint32),
            ('Ndy', ctypes.c_uint32),
            ('Ndz', ctypes.c_uint32),
            ('Nlx', ctypes.c_uint32),
            ('Nly', ctypes.c_uint32),
            ('Nlz', ctypes.c_uint32),
            ('g_dim_x', ctypes.c_uint32),
            ('g_dim_y', ctypes.c_uint32),
            ('g_dim_z', ctypes.c_uint32),
            ('NiterAD', ctypes.c_uint32),
            ('inffi', ctypes.c_uint32),
            ('Nf', ctypes.c_uint32),
            ('deviceNum', ctypes.c_uint32),
            ('platform', ctypes.c_uint32),
            ('derivativeType', ctypes.c_uint32),
            ('TVtype', ctypes.c_uint32),
            ('FluxType', ctypes.c_uint32),
            ('DiffusionType', ctypes.c_uint32),
            ('POCS_NgradIter', ctypes.c_uint32),
            ('maskFPZ', ctypes.c_uint32),
            ('maskBPZ', ctypes.c_uint32),
            ('FISTAType', ctypes.c_uint32),
            ('nProjections', ctypes.c_int64),
            ('TOF_bins', ctypes.c_int64),
            ('tau', ctypes.c_float),
            ('tube_radius', ctypes.c_float),
            ('epps', ctypes.c_float),
            ('sigma_x', ctypes.c_float),
            ('tube_width_z', ctypes.c_float),
            ('tube_width_xy', ctypes.c_float),
            ('bmin', ctypes.c_float),
            ('bmax', ctypes.c_float),
            ('Vmax', ctypes.c_float),
            ('global_factor', ctypes.c_float),
            ('dL', ctypes.c_float),
            ('flat', ctypes.c_float),
            ('U', ctypes.c_float),
            ('h_ACOSEM', ctypes.c_float),
            ('dPitchX', ctypes.c_float),
            ('dPitchY', ctypes.c_float),
            ('cr_p', ctypes.c_float),
            ('cr_pz', ctypes.c_float),
            ('NLMsigma', ctypes.c_float),
            ('NLAdaptiveConstant', ctypes.c_float),
            ('w_sum', ctypes.c_float),
            ('KAD', ctypes.c_float),
            ('TimeStepAD', ctypes.c_float),
            ('RDP_gamma', ctypes.c_float),
            ('huber_delta', ctypes.c_float),
            ('gradV1', ctypes.c_float),
            ('gradV2', ctypes.c_float),
            ('alpha0TGV', ctypes.c_float),
            ('alpha1TGV', ctypes.c_float),
            ('GGMRF_p', ctypes.c_float),
            ('GGMRF_q', ctypes.c_float),
            ('GGMRF_c', ctypes.c_float),
            ('beta', ctypes.c_float),
            ('T', ctypes.c_float),
            ('dSizeXBP', ctypes.c_float),
            ('dSizeZBP', ctypes.c_float),
            ('TVsmoothing', ctypes.c_float),
            ('C', ctypes.c_float),
            ('SATVPhi', ctypes.c_float),
            ('eta', ctypes.c_float),
            ('APLSsmoothing', ctypes.c_float),
            ('hyperbolicDelta', ctypes.c_float),
            ('sourceToCRot', ctypes.c_float),
            ('POCS_alpha', ctypes.c_float),
            ('POCS_rMax', ctypes.c_float),
            ('POCS_alphaRed', ctypes.c_float),
            ('POCSepps', ctypes.c_float),
            ('use_psf', ctypes.c_bool),
            ('TOF', ctypes.c_bool),
            ('pitch', ctypes.c_bool),
            ('SPECT', ctypes.c_bool),
            ('PET', ctypes.c_bool),
            ('CT', ctypes.c_bool),
            ('largeDim', ctypes.c_bool),
            ('loadTOF', ctypes.c_bool),
            ('storeResidual', ctypes.c_bool),
            ('FISTA_acceleration', ctypes.c_bool),
            ('meanFP', ctypes.c_bool),
            ('meanBP', ctypes.c_bool),
            ('useMaskFP', ctypes.c_bool),
            ('useMaskBP', ctypes.c_bool),
            ('orthTransaxial', ctypes.c_bool),
            ('orthAxial', ctypes.c_bool),
            ('enforcePositivity', ctypes.c_bool),
            ('useMultiResolutionVolumes', ctypes.c_bool),
            ('save_iter', ctypes.c_bool),
            ('deblurring', ctypes.c_bool),
            ('useMAD', ctypes.c_bool),
            ('useImages', ctypes.c_bool),
            ('useEFOV', ctypes.c_bool),
            ('CTAttenuation', ctypes.c_bool),
            ('offsetCorrection', ctypes.c_bool),
            ('relaxationScaling', ctypes.c_bool),
            ('computeRelaxationParameters', ctypes.c_bool),
            ('storeFP', ctypes.c_bool),
            ('use2DTGV', ctypes.c_bool),
            ('med_no_norm', ctypes.c_bool),
            ('NLM_MRP', ctypes.c_bool),
            ('NLTV', ctypes.c_bool),
            ('NLRD', ctypes.c_bool),
            ('NLLange', ctypes.c_bool),
            ('NLGGMRF', ctypes.c_bool),
            ('NLM_use_anatomical', ctypes.c_bool),
            ('NLAdaptive', ctypes.c_bool),
            ('TV_use_anatomical', ctypes.c_bool),
            ('RDPIncludeCorners', ctypes.c_bool),
            ('RDP_use_anatomical', ctypes.c_bool),
            ('useL2Ball', ctypes.c_bool),
            ('saveSens', ctypes.c_bool),
            ('use_64bit_atomics', ctypes.c_bool),
            ('use_32bit_atomics', ctypes.c_bool),
            ('compute_sensitivity_image', ctypes.c_bool),
            ('useFDKWeights', ctypes.c_bool),
            ('useIndexBasedReconstruction', ctypes.c_bool),
            ('stochasticSubsetSelection', ctypes.c_bool),
            ('useTotLength', ctypes.c_bool),
            ('OSEM', ctypes.c_bool),
            ('LSQR', ctypes.c_bool),
            ('CGLS', ctypes.c_bool),
            ('SART', ctypes.c_bool),
            ('FISTA', ctypes.c_bool),
            ('FISTAL1', ctypes.c_bool),
            ('MRAMLA', ctypes.c_bool),
            ('RAMLA', ctypes.c_bool),
            ('ROSEM', ctypes.c_bool),
            ('RBI', ctypes.c_bool),
            ('DRAMA', ctypes.c_bool),
            ('COSEM', ctypes.c_bool),
            ('ECOSEM', ctypes.c_bool),
            ('ACOSEM', ctypes.c_bool),
            ('OSL_OSEM', ctypes.c_bool),
            ('MBSREM', ctypes.c_bool),
            ('BSREM', ctypes.c_bool),
            ('ROSEM_MAP', ctypes.c_bool),
            ('OSL_RBI', ctypes.c_bool),
            ('OSL_COSEM', ctypes.c_bool),
            ('PKMA', ctypes.c_bool),
            ('SPS', ctypes.c_bool),
            ('PDHG', ctypes.c_bool),
            ('PDHGKL', ctypes.c_bool),
            ('PDHGL1', ctypes.c_bool),
            ('CV', ctypes.c_bool),
            ('PDDY', ctypes.c_bool),
            ('POCS', ctypes.c_bool),
            ('FDK', ctypes.c_bool),
            ('SAGA', ctypes.c_bool),
            ('MRP', ctypes.c_bool),
            ('quad', ctypes.c_bool),
            ('Huber', ctypes.c_bool),
            ('L', ctypes.c_bool),
            ('FMH', ctypes.c_bool),
            ('weighted_mean', ctypes.c_bool),
            ('TV', ctypes.c_bool),
            ('hyperbolic', ctypes.c_bool),
            ('AD', ctypes.c_bool),
            ('APLS', ctypes.c_bool),
            ('TGV', ctypes.c_bool),
            ('NLM', ctypes.c_bool),
            ('RDP', ctypes.c_bool),
            ('GGMRF', ctypes.c_bool),
            ('ProxTV', ctypes.c_bool),
            ('ProxRDP', ctypes.c_bool),
            ('ProxNLM', ctypes.c_bool),
            ('MAP', ctypes.c_bool),
            ('custom', ctypes.c_bool),
            ('mDim', ctypes.c_uint64),
            ('nIterSaved', ctypes.c_uint64),
            ('sizeScat', ctypes.c_uint64),
            ('eFOV', ctypes.c_uint64),
            ('sizeX', ctypes.c_uint64),
            ('sizeZ', ctypes.c_uint64),
            ('sizeAtten', ctypes.c_uint64),
            ('sizeNorm', ctypes.c_uint64),
            ('sizePSF', ctypes.c_uint64),
            ('sizeXYind', ctypes.c_uint64),
            ('sizeZind', ctypes.c_uint64),
            ('xCenterSize', ctypes.c_uint64),
            ('yCenterSize', ctypes.c_uint64),
            ('zCenterSize', ctypes.c_uint64),
            ('sizeV', ctypes.c_uint64),
            ('measElem', ctypes.c_uint64),
            ('x', ctypes.POINTER(ctypes.c_float)),
            ('z', ctypes.POINTER(ctypes.c_float)),
            ('uV', ctypes.POINTER(ctypes.c_float)),
            ('dx', ctypes.POINTER(ctypes.c_float)),
            ('dy', ctypes.POINTER(ctypes.c_float)),
            ('dz', ctypes.POINTER(ctypes.c_float)),
            ('bx', ctypes.POINTER(ctypes.c_float)),
            ('by', ctypes.POINTER(ctypes.c_float)),
            ('bz', ctypes.POINTER(ctypes.c_float)),
            ('atten', ctypes.POINTER(ctypes.c_float)),
            ('norm', ctypes.POINTER(ctypes.c_float)),
            ('pituus', ctypes.POINTER(ctypes.c_int64)),
            ('xy_index', ctypes.POINTER(ctypes.c_uint32)),
            ('z_index', ctypes.POINTER(ctypes.c_uint16)),
            ('x_center', ctypes.POINTER(ctypes.c_float)),
            ('y_center', ctypes.POINTER(ctypes.c_float)),
            ('z_center', ctypes.POINTER(ctypes.c_float)),
            ('V', ctypes.POINTER(ctypes.c_float)),
            ('gaussPSF', ctypes.POINTER(ctypes.c_float)),
            ('saveNiter', ctypes.POINTER(ctypes.c_uint32)),
            ('Nx', ctypes.POINTER(ctypes.c_uint32)),
            ('Ny', ctypes.POINTER(ctypes.c_uint32)),
            ('Nz', ctypes.POINTER(ctypes.c_uint32)),
            ('randoms', ctypes.POINTER(ctypes.c_float)),
            ('corrVector', ctypes.POINTER(ctypes.c_float)),
            ('x0', ctypes.POINTER(ctypes.c_float)),
            ('offsetVal', ctypes.POINTER(ctypes.c_float)),
            ('dScaleX4', ctypes.POINTER(ctypes.c_float)),
            ('dScaleY4', ctypes.POINTER(ctypes.c_float)),
            ('dScaleZ4', ctypes.POINTER(ctypes.c_float)),
            ('dSizeX', ctypes.POINTER(ctypes.c_float)),
            ('dSizeY', ctypes.POINTER(ctypes.c_float)),
            ('dScaleX', ctypes.POINTER(ctypes.c_float)),
            ('dScaleY', ctypes.POINTER(ctypes.c_float)),
            ('dScaleZ', ctypes.POINTER(ctypes.c_float)),
            ('kerroin4', ctypes.POINTER(ctypes.c_float)),
            ('lam_drama', ctypes.POINTER(ctypes.c_float)),
            ('maskFP', ctypes.POINTER(ctypes.c_uint8)),
            ('maskBP', ctypes.POINTER(ctypes.c_uint8)),
            ('eFOVIndices', ctypes.POINTER(ctypes.c_uint8)),
            ('maskPrior', ctypes.POINTER(ctypes.c_uint8)),
            ('TOFIndices', ctypes.POINTER(ctypes.c_uint8)),
            ('angles', ctypes.POINTER(ctypes.c_float)),
            ('blurPlanes', ctypes.POINTER(ctypes.c_uint32)),
            ('gFilter', ctypes.POINTER(ctypes.c_float)),
            ('gFSize', ctypes.POINTER(ctypes.c_uint64)),
            ('precondTypeImage', ctypes.POINTER(ctypes.c_bool)),
            ('precondTypeMeas', ctypes.POINTER(ctypes.c_bool)),
            ('referenceImage', ctypes.POINTER(ctypes.c_float)),
            ('filterIm', ctypes.POINTER(ctypes.c_float)),
            ('filter', ctypes.POINTER(ctypes.c_float)),
            ('filter2', ctypes.POINTER(ctypes.c_float)),
            ('Ffilter', ctypes.POINTER(ctypes.c_float)),
            ('s', ctypes.POINTER(ctypes.c_float)),
            ('weights_quad', ctypes.POINTER(ctypes.c_float)),
            ('weights_huber', ctypes.POINTER(ctypes.c_float)),
            ('weighted_weights', ctypes.POINTER(ctypes.c_float)),
            ('APLS_ref_image', ctypes.POINTER(ctypes.c_float)),
            ('lambdaN', ctypes.POINTER(ctypes.c_float)),
            ('lambdaFiltered', ctypes.POINTER(ctypes.c_float)),
            ('alpha_PKMA', ctypes.POINTER(ctypes.c_float)),
            ('alphaPrecond', ctypes.POINTER(ctypes.c_float)),
            ('NLM_ref', ctypes.POINTER(ctypes.c_float)),
            ('RDP_ref', ctypes.POINTER(ctypes.c_float)),
            ('tauCP', ctypes.POINTER(ctypes.c_float)),
            ('tauCPFilt', ctypes.POINTER(ctypes.c_float)),
            ('sigmaCP', ctypes.POINTER(ctypes.c_float)),
            ('sigma2CP', ctypes.POINTER(ctypes.c_float)),
            ('thetaCP', ctypes.POINTER(ctypes.c_float)),
            ('TOFCenter', ctypes.POINTER(ctypes.c_float)),
            ('TV_ref', ctypes.POINTER(ctypes.c_float)),
            ('trIndices', ctypes.POINTER(ctypes.c_uint16)),
            ('axIndices', ctypes.POINTER(ctypes.c_uint16)),
            ('crXY',ctypes.c_float),
            ('rayShiftsDetector',ctypes.POINTER(ctypes.c_float)),
            ('rayShiftsSource',ctypes.POINTER(ctypes.c_float)),
        ]
