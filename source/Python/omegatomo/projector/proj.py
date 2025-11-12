# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 17:40:13 2024

@author: Ville-Veikko Wettenhovi



#########################################################################
# Copyright (C) 2024-2025 Ville-Veikko Wettenhovi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#########################################################################
"""

# from SimpleITK import GetSpacing
import numpy as np
import numpy.typing as npt
import ctypes
import math
import os
from .coordinates import computePixelCenters
from .coordinates import computePixelSize
from .coordinates import computeProjectorScalingValues
from .coordinates import computeVoxelVolumes
from .indices import indexMaker
from .indices import formSubsetIndices

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
    angles: npt.NDArray[np.float32] = np.empty(0, dtype = np.float32)
    blurPlanes: npt.NDArray[np.int32] = np.empty(0, dtype = np.int32)
    blurPlanes2: npt.NDArray[np.int32] = np.empty(0, dtype = np.int32)
    radiusPerProj: npt.NDArray[np.float32] = np.empty(0, dtype = np.float32)
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
    gaussianNLM = np.empty(0, dtype = np.float32)
    xy_index = np.empty(0, dtype=np.uint32)
    z_index = np.empty(0, dtype=np.uint16)
    tauCPFilt = 0.
    sigmaCP = 1.
    tauCP = 0.
    sigma2CP = 1.
    thetaCP = 1.
    machine_name = ""
    fpath = ""
    filterWindow = "hamming"
    sampling_interpolation_method = "linear"
    arc_interpolation = "linear"
    scatterPath = ""
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
    ring_difference = 0
    obtain_trues = False
    reconstruct_trues = False
    store_scatter = False
    reconstruct_scatter = False
    store_randoms = False
    source = False
    weighted_center_weight = 2.
    offangle = 0.
    binning = 1
    sourceToCRot = 1.
    sourceToDetector = 1.
    nBed = 1
    TOF_FWHM = 0.
    TOF_width = 0.
    TOF_bins_used = 0
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
    corrections_during_reconstruction = True
    ordinaryPoisson = None
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
    colL: float = 0. # Collimator hole length
    colR: float = 0. # Collimator hole radius
    colD: float = 0. # Distance from collimator hole centre to detector surface
    colFxy: float = np.inf # Collimator focal distance, XY plane
    colFz: float = np.inf # Collimator focal distance, Z direction
    iR: float = 0. # Detector intrinsic resolution
    implementation = 2
    rotateAttImage = 0.
    flipAttImageXY = False
    flipAttImageZ = False
    NSinos = 0
    TotSinos = NSinos
    V = np.empty(0, dtype = np.float32)
    uu = 0
    ub = 0
    det_w_pseudo = 0
    startAngle = 0.
    nHeads = 1
    angleIncrement = 0.
    flatFieldScaling = 1.
    CT_attenuation = True
    scatter_correction = 0
    subtract_scatter = True
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
    empty_weight = False
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
    usePriorMask = False
    offsetCorrection = False
    tube_width_z = 0.
    tube_width_xy = 0.
    use_psf = False
    save_iter = False
    deblurring = False
    use_64bit_atomics = True
    use_32bit_atomics = False
    n_rays_transaxial = 1
    n_rays_axial = 1
    RDPIncludeCorners = False
    meanFP = False
    meanBP = False
    Niter = 4
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
    Ndist = 400
    Nang = 1
    NSinos = 1
    dL = 0
    epps = 1e-5
    cr_p = 0.
    cr_pz = 0.
    verbose = 1
    partitions = 1
    Nt = 1
    randoms_correction = 0
    additionalCorrection = 0
    attenuation_correction = 0
    normalization_correction = 0
    global_correction_factor = 1.
    rings = 0
    linear_multip = 0
    detectors = 0
    NxOrig = 1
    NyOrig = 1
    NzOrig = 1
    NxPrior = 1
    NyPrior = 1
    NzPrior = 1
    det_per_ring = 0
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
    Nlx = 1
    Nly = 1
    Nlz = 1
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
    useMultiResolutionVolumes = False
    nMultiVolumes = 0
    storeMultiResolution = False
    extrapLength = None
    axialExtrapolation = False
    transaxialExtrapolation = False
    useExtrapolationWeighting = False
    transaxialEFOV = False
    axialEFOV = False
    eFOVLength = None
    NLM_gauss = 2.
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
    rayShiftsDetector: npt.NDArray[np.float32] = np.empty(0, dtype=np.float32)
    rayShiftsSource: npt.NDArray[np.float32] = np.empty(0, dtype=np.float32)
    CORtoDetectorSurface: float = 0 # Detector swivel radius
    swivelAngles: npt.NDArray[np.float32] = np.empty(0, dtype = np.float32)
    coneOfResponseStdCoeffA = 0
    coneOfResponseStdCoeffB = 0
    coneOfResponseStdCoeffC = 0
    FISTAType = 0
    maskFPZ = 1
    maskBPZ = 1
    stochasticSubsetSelection = False
    useTotLength = True
    TOFIndices = np.empty(0, dtype = np.uint8)
    useParallelBeam = False
    builtin = True
    precorrect = False
    seed = -1
    useHelical = False
    helicalRadius = 1.

    def __init__(self):
        # C-struct
        self.param = self.parameters()
    def addProjector(self):
        if self.OSL_OSEM or self.MBSREM or self.ROSEM_MAP or self.OSL_RBI or self.OSL_COSEM > 0 or self.PKMA or self.SPS or self.PDHG or self.PDHGKL or self.PDHGL1 or self.PDDY or self.CV or self.SAGA or self.SART or self.ASD_POCS:
            self.MAP = True
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
            if self.projector_type == 6 and len(self.SinM > 0):
                endSinogramRows = self.FOVa_x / self.dPitchX; # Desired amount of sinogram rows
                endSinogramCols = self.axial_fov / self.dPitchY; # Desired amount of sinogram columns
                padRows = int((endSinogramRows-self.nRowsD)/2) # Pad this amount on both sides
                padCols = int((endSinogramCols-self.nColsD)/2) # Pad this amount on both sides
                if padRows < 0:
                    self.SinM = self.SinM[-padRows:padRows, :, :]
                if padRows > 0:
                    self.SinM = np.pad(self.SinM, ((padRows, padRows), (0, 0), (0, 0)))
                if padCols < 0:
                    self.SinM = self.SinM[:, -padCols:padCols, :]
                if padCols > 0:
                    self.SinM = np.pad(self.SinM, ((0, 0), (padCols, padCols), (0, 0)))
                self.nRowsD = self.Nx; # Set new sinogram size
                self.nColsD = self.Nz; # Set new sinogram size
                
                # Now the sinogram and FOV XZ-plane match in physical dimensions but not in resolution.
                from skimage.transform import resize
                self.SinM = resize(self.SinM, (self.Nx, self.Nz), preserve_range=True)
                
            if self.swivelAngles.size == 0:
                self.swivelAngles = self.angles + 180
            if self.offangle != 0:
                self.angles += self.offangle
                self.swivelAngles += self.offangle
            
            if self.vaimennus.size > 0:
                from skimage.transform import rotate
                self.vaimennus = rotate(self.vaimennus, self.offangle)
                if self.flipImageX:
                    self.vaimennus = np.flip(self.vaimennus, 2)
                if self.flipImageY:
                    self.vaimennus = np.flip(self.vaimennus, 1)
                if self.flipImageZ:
                    self.vaimennus = np.flip(self.vaimennus, 3)
                    
            if self.flipImageZ:
                self.SinM = np.flip(self.SinM, axis=2)
            
            self.n_rays_transaxial = self.nRays
            self.NSinos = self.nProjections
            self.TotSinos = self.nProjections
            self.det_per_ring = self.nRowsD * self.nProjections
            self.arc_correction = False
            self.det_w_pseudo = self.det_per_ring
            self.cryst_per_block = self.nColsD * self.nRowsD
            self.span = 3
            self.Ndist = self.nRowsD
            self.Nang = self.nColsD
            self.use_raw_data = 0
        
        if (type(self.cryst_per_block_axial) == np.ndarray and self.cryst_per_block_axial.size == 0) and ((type(self.cryst_per_block) == np.ndarray and self.cryst_per_block.size > 0) 
                                                                                                          or (type(self.cryst_per_block) != np.ndarray and self.cryst_per_block > 0)):
            self.cryst_per_block_axial = self.cryst_per_block
        if self.rings == 0 and self.cryst_per_block_axial >= 1 and self.linear_multip >= 1:
            self.rings = self.cryst_per_block_axial * self.linear_multip
        if self.dPitchX > 0. and self.dPitchY == 0.:
            self.dPitchY = self.dPitchX
        if self.cr_p > 0. and self.cr_pz == 0.:
            self.cr_pz = self.cr_p
        if self.cr_p == 0. and self.dPitchX > 0.:
            self.cr_p = self.dPitchX
        if self.cr_pz == 0. and self.dPitchY > 0.:
            self.cr_pz = self.dPitchY
        if self.axial_fov == 0. and (self.cr_pz > 0. or self.dPitchY > 0.):
            if self.dPitchY > 0.:
                self.axial_fov = self.dPitchY
            else:
                self.axial_fov = self.cr_pz
        rings = self.rings
        if self.det_per_ring == 0 and self.blocks_per_ring >= 1 and self.cryst_per_block >= 1:
            self.det_per_ring = self.blocks_per_ring * self.cryst_per_block
        if self.det_w_pseudo == 0:
            self.det_w_pseudo = self.det_per_ring
        if self.use_raw_data and self.x.size > 1:
            det_per_ring = self.x.size
        else:
            det_per_ring = self.det_per_ring
        if self.detectors == 0:
            self.detectors = self.det_w_pseudo * self.rings
        if self.ring_difference == 0 and self.rings > 0:
            self.ring_difference = self.rings - 1;
        if self.segment_table.size == 0 and self.span > 0 and self.rings > 0:
            self.segment_table = np.concatenate((np.array(self.rings*2-1,ndmin=1), np.arange(self.rings*2-1 - (self.span + 1), max(self.Nz - self.ring_difference*2, self.rings - self.ring_difference), -self.span*2)))
            self.segment_table = np.insert(np.repeat(self.segment_table[1:], 2), 0, self.segment_table[0])
        if self.Nang == 1 and self.det_per_ring > 1:
            self.Nang = self.det_per_ring // 2
        if self.span == 1:
            self.TotSinos = self.rings**2
            self.NSinos = self.TotSinos
        elif self.NSinos == 0:
            self.NSinos = np.sum(self.segment_table)
        if self.TotSinos == 0:
            self.TotSinos = self.NSinos
        self.OMEGAErrorCheck()
        self.TOF = self.TOF_bins > 1 and (self.projector_type == 1 or self.projector_type == 11 or self.projector_type == 3 or self.projector_type == 33 
                                                      or self.projector_type == 13 or self.projector_type == 31 or self.projector_type == 4 
                                                      or self.projector_type == 14 or self.projector_type == 41 or self.projector_type == 44 
                                                      or self.projector_type == 34 or self.projector_type == 43)
        if self.TOF:
            if self.TOF_bins_used == 0 and self.TOF_bins > 0:
                self.TOF_bins_used = self.TOF_bins
            if self.TOF_bins_used != self.TOF_bins:
                # self.SinM = np.sum(self.SinM,3)
                self.sigma_x = 0.
                # self.TOF = False
            elif self.TOFCenter.size == 0:
                c = 2.99792458e11 # speed of light in mm/s
                self.sigma_x = (c*self.TOF_FWHM/2.) / (2. * math.sqrt(2. * math.log(2.)))
                edges_user = np.linspace(-self.TOF_width * self.TOF_bins/2, self.TOF_width * self.TOF_bins / 2, self.TOF_bins + 1, dtype=np.float32)
                edges_user = edges_user[0:-1] + self.TOF_width/2. # the most probable value where annihilation occured
                self.TOFCenter = np.zeros(np.shape(edges_user),dtype = np.float32, order='F')
                self.TOFCenter[0] = edges_user[math.floor(np.size(edges_user)/2)]
                self.TOFCenter[1::2] = edges_user[math.floor(np.size(edges_user)/2) + 1:]
                self.TOFCenter[2::2] = edges_user[math.floor(np.size(edges_user)/2) - 1:  : -1]
                if self.TOF_offset > 0:
                    self.TOFCenter = self.TOFCenter + self.TOF_offset
                self.TOFCenter = -self.TOFCenter * c / 2.
            else:
                self.TOFCenter = np.float32(self.TOFCenter)
        else:
            self.sigma_x = 0.
        if self.ordinaryPoisson == None:
            self.ordinaryPoisson = self.corrections_during_reconstruction
        if self.maskFP.size > 1 and ((not(self.maskFP.size == (self.nRowsD * self.nColsD)) and not(self.maskFP.size == (self.nRowsD * self.nColsD * self.nProjections)) 
                                      and (self.CT == True or self.SPECT == True)) or (not(self.maskFP.size == (self.Nang * self.Ndist)) and not(self.maskFP.size == (self.Nang * self.Ndist * self.NSinos)) and self.CT == False)):
            if self.CT == True or self.SPECT == True:
                raise ValueError('Incorrect size for the forward projection mask! Must be the size of a single projection image [' + str(self.nRowsD) + ' ' + str(self.nColsD) + ']  or full stack of [' + str(self.nRowsD) + ' ' + str(self.nColsD) + ' ' + str(self.nProjections) + ']')
            else:
                raise ValueError('Incorrect size for the forward projection mask! Must be the size of a single sinogram image [' + str(self.Nang) + ' ' + str(self.Ndist) + '] or 3D stack [' + str(self.Nang) + ' ' + str(self.Ndist) + ' ' + str(self.NSinos) + ']')
        elif self.maskFP.size > 1 and (self.maskFP.size == (self.nRowsD * self.nColsD) or self.maskFP.size == (self.nRowsD * self.nColsD * self.nProjections)):
            self.useMaskFP = True
            if (self.maskFP.ndim == 3):
                self.maskFPZ = self.maskFP.shape[2]
            self.maskFP = np.asfortranarray(self.maskFP)
        else:
            self.useMaskFP = False
        if self.maskFP.dtype != np.uint8:
            self.maskFP = self.maskFP.astype(np.uint8)
        
        if self.maskBP.size > 1 and not(self.maskBP.size == self.Nx * self.Ny) and not(self.maskBP.size == self.Nx * self.Ny * self.Nz):
            raise ValueError('Incorrect size for the backward projection mask! Must be the size of a single image [' + str(self.Nx) + ' ' + str(self.Ny) + '] or 3D stack [' + str(self.Nx) + ' ' + str(self.Ny) + ' ' + str(self.Nz) + ']')
        elif self.maskBP.size == self.Nx * self.Ny or self.maskBP.size == self.Nx * self.Ny * self.Nz:
            self.useMaskBP = True
            if (self.maskBP.ndim == 3):
                self.maskBPZ = self.maskBP.shape[2]
            self.maskBP = np.asfortranarray(self.maskBP)
        else:
            self.useMaskBP = False
        if self.maskBP.dtype != np.uint8:
            self.maskBP = self.maskBP.astype(np.uint8)
        list_mode_format = False
        
        if self.use_raw_data:
            rings = rings - np.sum(self.pseudot)
            self.rings = rings
            # self.detectors = det_per_ring * rings;
        
        if isinstance(self.partitions, np.ndarray):
            if self.partitions.size >= 1:
                self.Nt = self.partitions.size
            # else:
            #     self.Nt = self.partitions[0].item()
        else:
            self.Nt = self.partitions
        if self.Nt > 1 and self.subsetType == 3:
            raise ValueError('Subset type 3 is not supported with dynamic data!')
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
        if (self.CT == False and self.SPECT == False) and (((self.subsetType > 7 or (self.subsetType == 0 and self.listmode == 0)) and self.subsets > 1) or self.subsets == 1):
            self.nProjections = self.NSinos
            self.PET = True
        else:
            self.PET = False
            self.nProjections = self.NSinos
        if self.listmode and self.subsets > 1 and not(self.subsetType == 0)  and not(self.subsetType == 1) and not(self.subsetType == 3):
            print('Only subset types 0, 1, and 3 are supported with list-mode data! Switching to subset type 0.')
            self.subsetType = 0
        if self.listmode and self.Nt > 1:
            self.loadTOF = False
        indexMaker(self)
        self.setUpCorrections()
        self.x0 = self.x0.ravel('F')
        if self.CT and self.projector_type == 11 and not self.useCPU:
            self.projector_type = 4
        
        # Coordinates of the detectors
        if self.projector_type != 6:
            if self.listmode == False:
                if self.SPECT:
                    from .detcoord import getCoordinatesSPECT
                    x_det, z_det = getCoordinatesSPECT(self)
                    y = 0
                else:
                    from .detcoord import getCoordinates
                    x_det, y, z_det = getCoordinates(self)
            else:
                if self.TOF:
                    if self.TOFIndices.size != self.SinM.size:
                        raise ValueError('The number of TOF indices does not correspond to the number of events!')
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
                self.dPitch = self.dPitchX
                self.dPitchY = self.dPitchY
                self.dPitchX = self.dPitchX
                self.cr_p = self.dPitchX
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
        if ((self.CT or self.PET or self.SPECT) and self.projector_type != 6) and self.listmode == 0:
            if self.subsetType >= 8 and self.subsets > 1 and not self.FDK:
                if self.CT:
                    x_det = np.reshape(x_det, (self.nProjections, 6))
                    x_det = x_det[self.index,:]
                    x_det = x_det.ravel('C')
                    if self.pitch:
                        z_det = np.reshape(z_det, (self.nProjections, 6))
                        z_det = z_det[self.index,:]
                        z_det = z_det.ravel('C')
                    elif self.useHelical:
                        z_det = z_det[self.index]
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
                if self.CT and not self.useHelical:
                    self.uV = self.uV[self.index,:]
        if self.listmode == 0 and self.projector_type != 6:
            if self.SPECT:
                self.x = x_det.ravel('F')
                self.z = z_det.ravel('F')
            else:
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

        if self.SPECT:
            from .detcoord import SPECTParameters
            SPECTParameters(self)
        if self.projector_type == 6:
            if self.subsets > 1 and (self.subsetType == 8 or self.subsetType == 9 or self.subsetType == 10 or self.subsetType == 11):
                self.angles = self.angles[self.index]
                self.swivelAngles = self.swivelAngles[self.index]
                self.radiusPerProj = self.radiusPerProj[self.index]
                self.blurPlanes = self.blurPlanes[self.index]
                self.blurPlanes2 = self.blurPlanes2[self.index]
            #self.gFilter = self.gFilter.ravel('F').astype(dtype=np.float32)
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
            else:
                self.x = np.asfortranarray(self.x)
        # Compute PSF kernel
        self.PSFKernel()
        self.N = self.Nx.astype(np.uint64) * self.Ny.astype(np.uint64) * self.Nz.astype(np.uint64)
        if not isinstance(self.FOVa_x, np.ndarray):
            self.FOVa_x = np.array(self.FOVa_x, dtype=np.float32, ndmin=1)
        if not isinstance(self.FOVa_y, np.ndarray):
            self.FOVa_y = np.array(self.FOVa_y, dtype=np.float32, ndmin=1)
        if not isinstance(self.axial_fov, np.ndarray):
            self.axial_fov = np.array(self.axial_fov, dtype=np.float32, ndmin=1)
        if self.projector_type in [2, 3, 22, 33, 13, 12, 31, 32, 21, 23, 42, 43, 34, 24]:
            if self.projector_type in [3, 33, 13, 31, 32, 23, 43, 34]:
                self.orthTransaxial = True
            elif (self.projector_type in [2, 22, 12, 21, 24, 42]) and (self.tube_width_xy > 0 or self.SPECT):
                self.orthTransaxial = True
            else:
                self.orthTransaxial = False
        if self.projector_type in [2, 3, 22, 33, 13, 12, 31, 32, 21, 23, 42, 43, 34, 24]:
            if self.projector_type in [3, 33, 13, 31, 32, 23, 43, 34]:
                self.orthAxial = True
            elif (self.projector_type in [2, 22, 12, 21, 24, 42]) and (self.tube_width_z > 0 or self.SPECT):
                self.orthAxial = True
            else:
                self.orthAxial = False
        if self.use_32bit_atomics and self.use_64bit_atomics:
            self.use_64bit_atomics = False
        self.x0 = self.x0.astype(dtype=np.float32)
        if isinstance(self.x, int):
            self.x = np.zeros(1, dtype=np.float32)
            self.z = np.zeros(1, dtype=np.float32)
        else:
            self.x = self.x.astype(dtype=np.float32)
            self.z = self.z.astype(dtype=np.float32)
        if self.listmode and self.randoms_correction:
            self.randoms_correction = False
        
        
    def OMEGAErrorCheck(self):
        if self.FOVa_x > 0 and self.FOVa_y == 0:
            self.FOVa_y = self.FOVa_x
        if not self.CT and not self.SPECT and (self.FOVa_x >= self.diameter or self.FOVa_y >= self.diameter) and self.diameter > 0:
            print(f"Transaxial FOV is larger than the scanner diameter ({self.diameter})!")
        if not self.CT and not self.SPECT and self.axial_fov < (self.rings * self.cr_pz - self.cr_pz) and self.rings > 0:
            print("Axial FOV is too small, crystal ring(s) on the boundary have no slices!")
        
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
        
        if not self.CT and not self.SPECT and self.ring_difference >= self.rings and not self.use_raw_data and self.rings > 0:
            print(f"Ring difference can be at most {self.rings - 1}. Setting it to the maximum possible.")
            self.ring_difference = self.rings - 1
        
        if not self.CT and not self.SPECT and self.ring_difference < 0 and not self.use_raw_data:
            raise ValueError("Ring difference has to be at least 0!")
        
        if not self.CT and not self.SPECT and self.Nang > self.det_w_pseudo / 2 and not self.use_raw_data and self.Nang > 1:
            raise ValueError(f"Number of sinogram angles can be at most the number of detectors per ring divided by two ({self.det_w_pseudo / 2})!")
        
        if not self.CT and not self.SPECT and self.TotSinos < self.NSinos and not self.use_raw_data and self.TotSinos > 0:
            print(f"The number of sinograms used ({self.NSinos}) is larger than the total number of sinograms ({self.TotSinos})! Setting the total number of sinograms to that of sinograms used!")
            self.TotSinos = self.NSinos
        
        if not self.CT and not self.SPECT and (self.ndist_side > 1 and self.Ndist % 2 == 0 or self.ndist_side < -1 and self.Ndist % 2 == 0) and not self.use_raw_data:
            raise ValueError("ndist_side can be either 1 or -1!")
        
        if not self.CT and not self.SPECT and self.ndist_side == 0 and self.Ndist % 2 == 0 and not self.use_raw_data:
            raise ValueError("ndist_side cannot be 0 when Ndist is even!")
        if self.useIndexBasedReconstruction and self.projector_type not in [1, 2, 3, 4, 11, 14, 12, 13, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44]:
            raise ValueError('Index-based reconstruction only supports projector types 1-4!')
        
        if self.Nt < 1:
            if isinstance(self.partitions, np.ndarray) and self.partitions.size > 0:
                self.Nt = self.partitions.size
            else:
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
        
        if self.useMultiResolutionVolumes and not self.useEFOV:
            print("Multi-resolution reconstruction selected, but extended FOV is not selected! Disabling multi-resolution volumes.")
            self.useMultiResolutionVolumes = False
        
        if (self.LSQR or self.CGLS) and self.useMultiResolutionVolumes:
            raise ValueError("Multi-resolution reconstruction is not supported with LSQR or CGLS!")
        if self.FDK and self.storeFP:
            print('Forward projections cannot be stored with FDK/FBP!')
            
        # if not self.CT and not self.SPECT and self.det_per_ring == self.det_w_pseudo and self.fill_sinogram_gaps:
        #     raise ValueError('Gap filling is only supported with pseudo detectors!')
        if not self.largeDim and not self.x0.any():
            if not self.CT:
                print('Initial value is an empty array, using the default values (1)')
                self.x0 = np.ones((self.Nx, self.Ny, self.Nz), dtype=np.float32, order='F')
            else:
                print('Initial value is an empty array, using the default values (1e-4)')
                self.x0 = np.ones((self.Nx, self.Ny, self.Nz), dtype=np.float32, order='F') * 1e-4
        if not self.largeDim and self.x0.size < self.Nx * self.Ny * self.Nz and not self.useEFOV:
            raise ValueError(f"Initial value has a matrix size smaller ({np.prod(self.x0.shape)}) than the actual image size ({self.Nx*self.Ny*self.Nz})!")
        if not self.largeDim and self.x0.size > self.Nx * self.Ny * self.Nz and not self.useEFOV:
            print(f"Initial value has a matrix size larger ({np.prod(self.x0.shape)}) than the actual image size ({self.Nx*self.Ny*self.Nz})! Attempting automatic resize.")
            try:
                from skimage.transform import resize #scikit-image
                self.x0 = resize(self.x0, (self.Nx, self.Ny, self.Nz))
            except ModuleNotFoundError:
                print('skimage package not found! Unable to perform automatic resize! Install scikit-image package with "pip install scikit-image".')
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
            print('Using TV type = 2, but no anatomical reference set. Using TV type = 1 instead!')
            self.TVtype == 1
        if self.projector_type not in [1, 2, 3, 4, 5, 6, 11, 14, 12, 13, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44, 45, 51, 15, 54, 55]:
            raise ValueError('The selected projector type is not supported!')
        if self.APLS and not os.path.exists(self.APLS_ref_image) and self.MAP and not type(self.APLS_ref_image) == np.ndarray:
            raise FileNotFoundError('APLS selected, but the anatomical reference image was not found on path!')
        if self.epps <= 0:
            print('Epsilon value is zero or less than zero; must be a positive value. Using the default value (1e-6).')
            self.epps = 1e-6
        # if self.projector_type in [5, 4, 34, 24, 44, 14, 41, 42, 43, 15, 45, 54, 51, 55] and self.useCPU:
        #     raise ValueError('The selected projector type is not supported with CPU implementation!')
        if np.sum(self.precondTypeImage) == 0 and (self.PKMA or self.MRAMLA or self.MBSREM):
            print('No image-based preconditioner selected with PKMA/MRAMLA/MBSREM. EM preconditioner is highly recommended!')
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
        
        if self.useCPU and self.projector_type in [2, 3, 4, 5, 15, 25, 35, 45, 21, 22, 23, 24, 31, 32, 33, 34, 35, 41, 42, 43, 44, 45, 51, 52, 53, 54, 55, 14, 12, 13]:
            raise ValueError('Selected projector type is not supported with CPU implementation!')
        if self.projector_type == 2 and self.CT:
            raise ValueError('Orthogonal distance-based projector is NOT supported when using CT data!')
        
        if self.projector_type in [5, 15, 51, 45, 54] and not self.CT:
            raise ValueError('Projector type 5 is only supported with CT data!')
        
        if self.projector_type == 6 and not self.SPECT:
            raise ValueError('Projector type 6 is only supported with SPECT data!')
            
        if (not(self.projector_type == 6) and not(self.projector_type == 1) and not(self.projector_type == 11) and not(self.projector_type == 2) and not(self.projector_type == 22)) and self.SPECT:
            raise ValueError('SPECT only supports projector types 1, 2 and 6!')
        
        if self.projector_type == 6:
            if self.subsets > 1 and self.subsetType < 8:
                raise ValueError('Subset types 0-7 are not supported with projector type 6!')
        
        if self.FDK and (self.Niter > 1 or self.subsets > 1):
            if self.largeDim:
                self.Niter = 1
            else:
                print('When using FDK/FBP, the number of iterations and subsets must be set as 1. Setting both to 1.')
                self.subsets = 1
                self.Niter = 1
        if self.largeDim:
            if not self.PDHG and not self.FDK and not self.PKMA and not self.PDDY and not self.PDHGL1 and not self.PDHGKL:
                raise ValueError('Large dimension support is only available for PDHG, PKMA, and FDK!')
            if self.MRP or self.quad or self.Huber or self.weighted_mean or self.FMH or self.ProxTV or self.TGV or self.L or self.AD:
                raise ValueError('Large dimension support is only available for non-local methods, RDP, GGMRF, hyperbolic prior and TV!')
        if self.useCUDA and self.useCPU:
            raise ValueError('Both CUDA and CPU selected! Select only one!')
        
        if self.TOF_bins_used > 1 and (self.projector_type not in [1, 11, 3, 33, 31, 13, 4, 41, 14, 43, 34]) and not self.CT and not self.SPECT:
            raise ValueError('TOF is currently only supported with improved Siddon (projector_type = 1), interpolation-based projector (projector_type = 4) and volume of intersection (projector_type = 3)')
        
        if self.TOF_bins_used > 1 and self.TOF_width <= 0 and not self.CT and not self.SPECT:
            raise ValueError('TOF width (self.TOF_width) must be greater than zero.')
        
        if self.TOF_bins_used > 1 and self.TOF_FWHM == 0 and not self.CT and not self.SPECT:
            raise ValueError('TOF enabled, but the TOF FWHM (self.TOF_FWHM) is zero. FWHM must be nonzero.')
        
        if self.corrections_during_reconstruction and (self.scatter_correction or self.randoms_correction) and (self.PDHG or self.PDHGL1 or self.FISTA or self.LSQR or self.CGLS or self.FISTAL1 or self.PDDY):
            print('Randoms/scatter correction cannot be applied during the reconstruction with the selected algorithm! Attempting precorrection!')
            self.ordinaryPoisson = False
            # self.scatter_correction = False
            # self.randoms_correction = False     
        if self.useIndexBasedReconstruction and (self.randoms_correction or self.scatter_correction) and self.TOF_bins > 1:
            raise ValueError('Randoms and/or scatter correction cannot be used with index-based reconstruction with TOF data!')
        if self.subsetType == 3 and self.maskFP.size > 1:
            raise ValueError('Forward projection mask is not supported with subset type 3!')
        
        if self.precondTypeMeas[1] and (self.subsetType < 8 and not(self.subsetType == 4) and not(self.subsetType == 0)):
            raise ValueError('Filtering-based preconditioner only works with subset types 0, 4 and 8-11!')
        if self.useCPU and self.ECOSEM:
            raise ValueError('ECOSEM is not supported with CPU!')
        if self.NLM_use_anatomical and self.useCPU and self.NLM:
            raise ValueError('Reference image weighting for NLM is not supported with CPU!')
        if self.useIndexBasedReconstruction and self.useCPU:
            raise ValueError('Index-based reconstruction is not supported on CPU!')
        if self.useCPU:
            print('CPU functionality is limited and might not work correctly in all cases! Use at your own risk!')
        if self.useHelical and not self.projector_type == 4:
            raise ValueError('Only projector type 4 is supported with curved helical data!')
        if self.offangle > 0 and self.useHelical:
            print('Rotation is not yet supported with helical CT data')
            
        varNeg = ['LSQR','CGLS','FDK','SART']
        neg = [name for name in varNeg if getattr(self, name, False)]
        if len(neg) > 0:
            self.enforcePositivity = False
            
        if self.verbose > 0:
            if self.use_ASCII and self.use_machine == 0:
                dispi = 'Using ASCII data'
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
                dispi = 'Using user-input data'
            
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
                    import arrayfire as af
                    AFinstalled = True
                except ModuleNotFoundError:
                    print('ArrayFire package not found! ArrayFire features are not supported. You can install ArrayFire package with "pip install arrayfire".')
                    AFinstalled = False
                if AFinstalled and not self.useCPU:
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
                elif self.useCPU:
                    print('Using CPU-based reconstruction')
                else:
                    print('Selected device number is ' + str(self.deviceNum))
                    
                algorithms = [
                    "OSEM", "MRAMLA", "RAMLA", "RBI", "ROSEM", "DRAMA", "COSEM", "ECOSEM", "ACOSEM", "LSQR", "CGLS", "FDK", "FISTA", "FISTAL1",
                    "OSL_OSEM", "MBSREM", "BSREM", "OSL_RBI", "OSL_COSEM", "ROSEM_MAP", "PKMA", "SART", "ASD_POCS", "SAGA",
                    "PDHG", "PDHGL1", "PDDY", "PDHGKL", "CV" ]
                
                
                varPrior = ['MRP','quad','Huber','L','FMH','weighted_mean','TV','hyperbolic','AD','APLS','TGV','NLM','RDP','GGMRF','ProxTV','ProxRDP','ProxNLM','custom']
                priors = [name for name in varPrior if getattr(self, name, False)]
                
                enabled_algorithms = [name for name in algorithms if getattr(self, name, False)]
                
                # Check how many are True
                if len(enabled_algorithms) == 0:
                    print("No reconstruction method selected! If you are using custom reconstruction, ignore this warning")
                    self.builtin = False
                elif len(enabled_algorithms) > 1:
                    raise ValueError(f"Multiple reconstruction algorithms selected: {enabled_algorithms}")
                else:
                    print(f"{enabled_algorithms[0]} reconstruction method selected.")
                    
                if len(priors) > 1:
                    raise ValueError(f"Multiple priors selected: {priors}")
                
                if len(priors) > 0:
                    if not self.NLM:
                        print(f"{priors[0]} prior selected.")
                    else:
                        if self.NLTV:
                            print(f"{priors[0]} prior selected with NLTV.")
                        elif self.NLRD:
                            print(f"{priors[0]} prior selected with NLRD.")
                        elif self.NLLange:
                            print(f"{priors[0]} prior selected with NL Lange.")
                        elif self.NLGGMRF:
                            print(f"{priors[0]} prior selected with NLGGMRF.")
                        elif self.NLM_MRP:
                            print(f"{priors[0]} prior selected with filtering mode.")
                        else:
                            print(f"{priors[0]} prior selected.")
                        if self.NLAdaptive:
                            print('Using adaptive weighting.')
                    
                
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
                elif self.projector_type == 3 or self.projector_type == 33:
                    print('Volume of intersection based ray tracer selected.');
                elif self.projector_type == 4:
                    print('Interpolation-based projector selected.')
                elif self.projector_type == 5:
                    print('Branchless distance-driven based projector selected.')
                elif self.projector_type == 6:
                    print('Rotation-based projector selected (SPECT).')
                elif self.projector_type > 10:
                    fpType = self.projector_type // 10
                    bpType = self.projector_type - fpType * 10
                    dispi = ""
                    if fpType == 1:
                        dispi += "Improved Siddon's algorithm selected for forward projection, "
                    elif fpType == 2:
                        dispi += "Orthogonal distance-based ray-tracer selected for forward projection, "
                    elif fpType == 3:
                        dispi += "Volume of intersection based ray tracer selected for forward projection, "
                    elif fpType == 4:
                        dispi += "Interpolation-based projector selected for forward projection, "
                    elif fpType == 5:
                        dispi += "Branchless distance-driven projector selected for forward projection, "
                    if bpType == 1:
                        dispi += "improved Siddon's algorithm for backprojection."
                    elif bpType == 2:
                        dispi += "orthogonal distance-based ray-tracer for backprojection."
                    elif bpType == 3:
                        dispi += "volume of intersection based ray tracer for backprojection."
                    elif bpType == 4:
                        dispi += "interpolation-based projector for backprojection."
                    elif bpType == 5:
                        dispi += "branchless distance-driven projector for backprojection."
                    print(dispi)
        
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
                    if self.subsets % 100 in [11, 12, 13]:
                        abbr = 'th'
                    elif self.subsets % 10 == 1:
                        abbr = 'st'
                    elif self.subsets % 10 == 2:
                        abbr = 'nd'
                    elif self.subsets % 10 == 3:
                        abbr = 'rd'
                    else:
                        abbr = 'th'
                    if self.subsetType == 1:
                        print(f'Every {self.subsets}{abbr} column measurement is taken per subset.')
                    elif self.subsetType == 2:
                        print(f'Every {self.subsets}{abbr} row measurement is taken per subset.')
                    elif self.subsetType == 3:
                        print('Using random subset sampling.')
                    elif self.subsetType == 4:
                        print(f'Every {self.subsets}{abbr} sinogram column is taken per subset.')
                    elif self.subsetType == 5:
                        print(f'Every {self.subsets}{abbr} sinogram row is taken per subset.')
                    elif self.subsetType == 6:
                        print(f'Using angle-based subset sampling with {self.n_angles} angles combined per subset.')
                    elif self.subsetType == 7:
                        print('Using golden angle-based subset sampling.')
                    elif self.subsetType == 8:
                        print(f'Using every {self.subsets}{abbr} sinogram/projection image.')
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
                if self.axialEFOV and self.transaxialEFOV:
                    from scipy.ndimage import zoom
                    self.nMultiVolumes = 6
                
                    NxM = round(self.NxOrig * self.multiResolutionScale)
                    dxM = self.FOVxOrig / NxM
                    NyM = round(self.NyOrig * self.multiResolutionScale)
                    dyM = self.FOVyOrig / NyM
                    NzM = round(self.NzOrig * self.multiResolutionScale)
                    dzM = self.axialFOVOrig / NzM
                
                    NxM2 = round((self.FOVa_x - self.FOVxOrig) / 2 / dxM) * 2
                    FOVxM = NxM2 * dxM
                    NyM2 = round((self.FOVa_y - self.FOVyOrig) / 2 / dyM) * 2
                    FOVyM = NyM2 * dyM
                    NzM2 = round((self.axial_fov - self.axialFOVOrig) / 2 / dzM) * 2
                    FOVzM = NzM2 * dzM
                
                    self.FOVa_x = np.array([
                        self.FOVxOrig, self.FOVxOrig, self.FOVxOrig,
                        FOVxM / 2, FOVxM / 2,
                        self.FOVxOrig, self.FOVxOrig
                    ], dtype=np.float32)
                
                    self.FOVa_y = np.array([
                        self.FOVyOrig, self.FOVyOrig, self.FOVyOrig,
                        FOVyM + self.FOVyOrig, FOVyM + self.FOVyOrig,
                        FOVyM / 2, FOVyM / 2
                    ], dtype=np.float32)
                
                    self.axial_fov = np.array([
                        self.axialFOVOrig, FOVzM / 2, FOVzM / 2,
                        FOVzM + self.axialFOVOrig, FOVzM + self.axialFOVOrig,
                        FOVzM + self.axialFOVOrig, FOVzM + self.axialFOVOrig
                    ], dtype=np.float32)
                
                    self.Nx = np.array([
                        self.NxOrig, NxM, NxM,
                        NxM2 // 2, NxM2 // 2,
                        NxM, NxM
                    ], dtype=np.uint32)
                
                    self.Ny = np.array([
                        self.NyOrig, NyM, NyM,
                        NyM + NyM2, NyM + NyM2,
                        NyM2 // 2, NyM2 // 2
                    ], dtype=np.uint32)
                
                    self.Nz = np.array([
                        self.NzOrig, NzM2 // 2, NzM2 // 2,
                        NzM + NzM2, NzM + NzM2,
                        NzM + NzM2, NzM + NzM2
                    ], dtype=np.uint32)
                
                    print(f"Extended FOV is {(FOVxM + self.FOVxOrig) / self.FOVxOrig * 100:.2f} % of the original")
                
                    if self.x0.shape[0] == nx and np.min(self.x0) != np.max(self.x0):
                        apu = zoom(self.x0, self.multiResolutionScale, order=1).astype(np.float32)
                
                        x1 = apu[
                            self.Nx[3]:self.Nx[3] + self.Nx[1],
                            self.Ny[5]:self.Ny[5] + self.Ny[1],
                            :self.Nz[1]
                        ]
                        x2 = apu[self.Nx[4]:self.Nx[4] + self.Nx[2],
                            self.Ny[6]:self.Ny[6] + self.Ny[2],
                            -self.Nz[1]:]
                        x3 = apu[:self.Nx[3], :, :]
                        if apu.shape[0] % 2 == 0:
                            x4 = apu[self.Nx[4] + self.Nx[2]:, :, :]
                        else:
                            x4 = apu[1 + self.Nx[4] + self.Nx[2]:, :, :]
                        x5 = apu[
                            self.Nx[3]:self.Nx[3] + self.Nx[1],
                            :self.Ny[5],
                            :self.Nz[3]
                        ]
                        if apu.shape[1] % 2 == 0:
                            x6 = apu[
                                self.Nx[4]:self.Nx[4] + self.Nx[2],
                                self.Ny[6] + self.Ny[2]:,
                                :self.Nz[4]
                            ]
                        else:
                            x6 = apu[
                                self.Nx[4]:self.Nx[4] + self.Nx[2],
                                1 + self.Ny[6] + self.Ny[2]:,
                                :self.Nz[4]
                            ]
                
                        sx0 = self.x0.shape
                        self.x0 = self.x0[
                            int((sx0[0] - self.NxOrig) // 2):int((sx0[0] - self.NxOrig) // 2 + self.NxOrig),
                            int((sx0[1] - self.NyOrig) // 2):int((sx0[1] - self.NyOrig) // 2 + self.NyOrig),
                            int((sx0[2] - self.NzOrig) // 2):int((sx0[2] - self.NzOrig) // 2 + self.NzOrig)
                        ].astype(np.float32)
                
                        self.x0 = np.concatenate([
                            self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F'),
                            x3.ravel('F'), x4.ravel('F'), x5.ravel('F'),
                            x6.ravel('F')
                        ])
                
                    elif self.x0.shape[0] == self.NxOrig or np.min(self.x0) == np.max(self.x0):
                        val = np.min(self.x0)
                        self.x0 = self.x0[:self.Nx[0].item(), :self.Ny[0].item(),:self.Nz[0].item()]
                        x1 = np.ones((self.Nx[1], self.Ny[1], self.Nz[1]), dtype=np.float32, order='F') * val
                        x2 = np.ones((self.Nx[2], self.Ny[2], self.Nz[2]), dtype=np.float32, order='F') * val
                        x3 = np.ones((self.Nx[3], self.Ny[3], self.Nz[3]), dtype=np.float32, order='F') * val
                        x4 = np.ones((self.Nx[4], self.Ny[4], self.Nz[4]), dtype=np.float32, order='F') * val
                        x5 = np.ones((self.Nx[5], self.Ny[5], self.Nz[5]), dtype=np.float32, order='F') * val
                        x6 = np.ones((self.Nx[6], self.Ny[6], self.Nz[6]), dtype=np.float32, order='F') * val
                        self.x0 = np.concatenate([
                            self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F'),
                            x3.ravel('F'), x4.ravel('F'), x5.ravel('F'),
                            x6.ravel('F')])
                
                elif self.transaxialEFOV and not self.axialEFOV:
                    self.nMultiVolumes = 4
                
                    NxM = round(self.NxOrig * self.multiResolutionScale)
                    dxM = self.FOVxOrig / NxM
                    NyM = round(self.NyOrig * self.multiResolutionScale)
                    dyM = self.FOVyOrig / NyM
                    NzM = round(self.NzOrig * self.multiResolutionScale)
                
                    NxM2 = round((self.FOVa_x - self.FOVxOrig) / 2 / dxM) * 2
                    FOVxM = NxM2 * dxM
                    NyM2 = round((self.FOVa_y - self.FOVyOrig) / 2 / dyM) * 2
                    FOVyM = NyM2 * dyM
                
                    self.FOVa_x = np.array([self.FOVxOrig, FOVxM / 2, FOVxM / 2, self.FOVxOrig, self.FOVxOrig], dtype=np.float32)
                    self.FOVa_y = np.array([self.FOVyOrig, FOVyM + self.FOVyOrig, FOVyM + self.FOVyOrig, FOVyM / 2, FOVyM / 2], dtype=np.float32)
                    self.axial_fov = np.array([self.axialFOVOrig] * 5, dtype=np.float32)
                
                    self.Nx = np.array([self.NxOrig, NxM2 // 2, NxM2 // 2, NxM, NxM], dtype=np.uint32)
                    self.Ny = np.array([self.NyOrig, NyM + NyM2, NyM + NyM2, NyM2 // 2, NyM2 // 2], dtype=np.uint32)
                    self.Nz = np.array([self.NzOrig, NzM, NzM, NzM, NzM], dtype=np.uint32)
                
                    print(f"Extended FOV is {(FOVxM + self.FOVxOrig) / self.FOVxOrig * 100:.2f} % of the original")
                
                    if self.x0.shape[0] == nx and np.min(self.x0) != np.max(self.x0):
                        apu = zoom(self.x0, self.multiResolutionScale, order=1).astype(np.float32)
                        self.x1 = apu[self.Nx[1]:self.Nx[1] + self.Nx[0], :, :]
                        if apu.shape[0] % 2 == 0:
                            self.x2 = apu[self.Nx[2] + self.Nx[0]:, :, :]
                        else:
                            self.x2 = apu[1 + self.Nx[2] + self.Nx[0]:, :, :]
                
                        sx0 = self.x0.shape
                        self.x0 = self.x0[
                            int((sx0[0] - self.NxOrig) // 2):int((sx0[0] - self.NxOrig) // 2 + self.NxOrig),
                            int((sx0[1] - self.NyOrig) // 2):int((sx0[1] - self.NyOrig) // 2 + self.NyOrig),
                            int((sx0[2] - self.NzOrig) // 2):int((sx0[2] - self.NzOrig) // 2 + self.NzOrig)
                        ].astype(np.float32)
                
                        self.x0 = np.concatenate([self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F')])
                
                    elif self.x0.shape[0] == self.NxOrig or np.min(self.x0) == np.max(self.x0):
                        val = np.min(self.x0)
                        self.x0 = self.x0[:self.Nx[0].item(), :self.Ny[0].item(),:]
                        x1 = np.ones((self.Nx[1], self.Ny[1], self.Nz[1]), dtype=np.float32, order='F') * val
                        x2 = np.ones((self.Nx[2], self.Ny[2], self.Nz[2]), dtype=np.float32, order='F') * val
                        self.x0 = np.concatenate([self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F')])
                elif not self.transaxialEFOV and self.axialEFOV:
                    self.nMultiVolumes = 2
                
                    NxM = round(self.NxOrig * self.multiResolutionScale)
                    NyM = round(self.NyOrig * self.multiResolutionScale)
                    NzM = round(self.NzOrig * self.multiResolutionScale)
                    dzM = self.axialFOVOrig / NzM
                    NzM2 = round((self.axial_fov - self.axialFOVOrig) / 2 / dzM) * 2
                    FOVzM = NzM2 * dzM
                
                    self.FOVa_x = np.array([self.FOVxOrig] * 3, dtype=np.float32)
                    self.FOVa_y = np.array([self.FOVyOrig] * 3, dtype=np.float32)
                    self.axial_fov = np.array([
                        self.axialFOVOrig, FOVzM / 2, FOVzM / 2
                    ], dtype=np.float32)
                
                    self.Nx = np.array([self.NxOrig, NxM, NxM], dtype=np.uint32)
                    self.Ny = np.array([self.NyOrig, NyM, NyM], dtype=np.uint32)
                    self.Nz = np.array([self.NzOrig, NzM2 // 2, NzM2 // 2], dtype=np.uint32)
                
                    print(f"Extended FOV is {(FOVzM + self.axialFOVOrig) / self.axialFOVOrig * 100:.2f} % of the original")
                
                    if self.x0.shape[0] == nx and np.min(self.x0) != np.max(self.x0):
                        apu = zoom(self.x0, self.multiResolutionScale, order=1).astype(np.float32)
                        x1 = apu[:, :, :self.Nz[1]]
                        x2 = apu[:, :, -self.Nz[2]:]
                
                        sx0 = self.x0.shape
                        self.x0 = self.x0[
                            int((sx0[0] - self.NxOrig) // 2):int((sx0[0] - self.NxOrig) // 2 + self.NxOrig),
                            int((sx0[1] - self.NyOrig) // 2):int((sx0[1] - self.NyOrig) // 2 + self.NyOrig),
                            int((sx0[2] - self.NzOrig) // 2):int((sx0[2] - self.NzOrig) // 2 + self.NzOrig)
                        ].astype(np.float32)
                
                        self.x0 = np.concatenate([self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F')])
                
                    elif self.x0.shape[0] == self.NxOrig or np.min(self.x0) == np.max(self.x0):
                        val = np.min(self.x0)
                        self.x0 = self.x0[:, :,:self.Nz[0].item()]
                        x1 = np.ones((self.Nx[1], self.Ny[1], self.Nz[1]), dtype=np.float32, order='F') * val
                        x2 = np.ones((self.Nx[2], self.Ny[2], self.Nz[2]), dtype=np.float32, order='F') * val
                        self.x0 = np.concatenate([self.x0.ravel('F'), x1.ravel('F'), x2.ravel('F')])
                
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
        else:
            self.NxPrior = self.Nx
            self.NyPrior = self.Ny
            self.NzPrior = self.Nz
        if self.offsetCorrection:
            self.OffsetLimit = np.zeros((self.nProjections, 1), dtype=np.float32);
            for kk in range(self.nProjections):
                sx = self.x[kk, 0];
                sy = self.x[kk, 1];
                dx = self.x[kk, 3];
                dy = self.x[kk, 4];
                ii = np.arange(0., self.nRowsD, 0.25, dtype=np.float32) - self.nRowsD / 2
                if self.pitch:
                    dx = dx + self.z[kk, 0] * ii + self.z[kk, 3] * ii
                    dy = dy + self.z[kk, 1] * ii + self.z[kk, 4] * ii
                else:
                    dx = dx + self.z[kk, 0] * ii
                    dy = dy + self.z[kk, 1] * ii
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
        from omegatomo.projector.init import initProjector
        initProjector(self)
                
    def computeConvolution(self, f, ii = 0):
        from omegatomo.projector.projfunctions import conv3D
        return conv3D(self, f, ii)
        
    def forwardProject(self, f, subset = -1):
        from omegatomo.projector.projfunctions import forwardProjection
        return forwardProjection(self, f, subset)
        
    def backwardProject(self, y, subset = -1):
        from omegatomo.projector.projfunctions import backwardProjection
        return backwardProjection(self, y, subset)
    
    
    def T(self):
        self.trans = True
        return self
        
    
    def __imul__(self, B):
            return self.__mul__(self, B)
    
    def __mul__(self, B):
        if self.trans:
            from omegatomo.projector.projfunctions import backwardProjection
            self.trans = False
            return backwardProjection(self, B, self.subset)
        else:
            from omegatomo.projector.projfunctions import forwardProjection
            return forwardProjection(self, B, self.subset)
        
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
            ('helicalRadius', ctypes.c_float),
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
            ('useParallelBeam', ctypes.c_bool),
            ('useHelical', ctypes.c_bool),
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
            ('PDDY', ctypes.c_bool),
            ('CV', ctypes.c_bool),
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
            ('seed', ctypes.c_int64),
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
            ('gaussianNLM', ctypes.POINTER(ctypes.c_float)),
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
            ('swivelAngles', ctypes.POINTER(ctypes.c_float)),
            ('blurPlanes', ctypes.POINTER(ctypes.c_int32)),
            ('blurPlanes2', ctypes.POINTER(ctypes.c_int32)),
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
            ('tauCP', ctypes.POINTER(ctypes.c_float)),
            ('tauCPFilt', ctypes.POINTER(ctypes.c_float)),
            ('sigmaCP', ctypes.POINTER(ctypes.c_float)),
            ('sigma2CP', ctypes.POINTER(ctypes.c_float)),
            ('thetaCP', ctypes.POINTER(ctypes.c_float)),
            ('TOFCenter', ctypes.POINTER(ctypes.c_float)),
            ('TV_ref', ctypes.POINTER(ctypes.c_float)),
            ('trIndices', ctypes.POINTER(ctypes.c_uint16)),
            ('axIndices', ctypes.POINTER(ctypes.c_uint16)),
            ('rayShiftsDetector',ctypes.POINTER(ctypes.c_float)),
            ('rayShiftsSource',ctypes.POINTER(ctypes.c_float)),
            ('coneOfResponseStdCoeffA',ctypes.c_float),
            ('coneOfResponseStdCoeffB',ctypes.c_float),
            ('coneOfResponseStdCoeffC',ctypes.c_float),
            ('NLM_ref', ctypes.POINTER(ctypes.c_float)),
            ('RDP_ref', ctypes.POINTER(ctypes.c_float)),
        ]
