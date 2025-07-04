# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:50:49 2024

@author: Ville-Veikko Wettenhovi
"""
import time
import numpy as np
import ctypes
import os
from .prepass import prepassPhase
from .prepass import parseInputs
from .prepass import loadCorrections
from .prepass import sinogramToX

def transferData(options):
    options.param.use_raw_data = ctypes.c_uint8(options.use_raw_data)
    options.param.listmode = ctypes.c_uint8(options.listmode)
    options.param.verbose = ctypes.c_int8(options.verbose)
    options.param.n_rays_transaxial = ctypes.c_uint16(options.n_rays_transaxial)
    options.param.n_rays_axial = ctypes.c_uint16(options.n_rays_axial)
    options.param.projector_type = ctypes.c_uint32(options.projector_type)
    options.param.attenuation_correction = ctypes.c_uint32(options.attenuation_correction)
    options.param.additionalCorrection = ctypes.c_uint32(options.additionalCorrection)
    options.param.normalization_correction = ctypes.c_uint32(options.normalization_correction)
    options.param.randoms_correction = ctypes.c_uint32(options.randoms_correction)
    options.param.nColsD = ctypes.c_uint32(options.nColsD)
    options.param.nRowsD = ctypes.c_uint32(options.nRowsD)
    options.param.Nang = ctypes.c_uint32(options.Nang)
    options.param.Ndist = ctypes.c_uint32(options.Ndist)
    options.param.subsets = ctypes.c_uint32(options.subsets)
    options.param.det_per_ring = ctypes.c_uint32(options.det_per_ring)
    options.param.rings = ctypes.c_uint32(options.rings)
    options.param.NxOrig = ctypes.c_uint32(options.NxOrig)
    options.param.NyOrig = ctypes.c_uint32(options.NyOrig)
    options.param.NzOrig = ctypes.c_uint32(options.NzOrig)
    options.param.NxPrior = ctypes.c_uint32(options.NxPrior)
    options.param.NyPrior = ctypes.c_uint32(options.NyPrior)
    options.param.NzPrior = ctypes.c_uint32(options.NzPrior)
    options.param.Niter = ctypes.c_uint32(options.Niter)
    options.param.Nt = ctypes.c_uint32(options.Nt)
    options.param.subsetType = ctypes.c_uint32(options.subsetType)
    options.param.nMultiVolumes = ctypes.c_uint32(options.nMultiVolumes)
    options.param.nLayers = ctypes.c_uint32(options.nLayers)
    options.param.PDAdaptiveType = ctypes.c_uint32(options.PDAdaptiveType)
    options.param.powerIterations = ctypes.c_uint32(options.powerIterations)
    options.param.deblur_iterations = ctypes.c_uint32(options.deblur_iterations)
    options.param.gradInitIter = ctypes.c_uint32(options.gradInitIter)
    options.param.gradLastIter = ctypes.c_uint32(options.gradLastIter)
    options.param.filteringIterations = ctypes.c_uint32(options.filteringIterations)
    options.param.mean_type = ctypes.c_uint32(options.mean_type)
    options.param.Ndx = ctypes.c_uint32(options.Ndx)
    options.param.Ndy = ctypes.c_uint32(options.Ndy)
    options.param.Ndz = ctypes.c_uint32(options.Ndz)
    options.param.Nlx = ctypes.c_uint32(options.Nlx)
    options.param.Nly = ctypes.c_uint32(options.Nly)
    options.param.Nlz = ctypes.c_uint32(options.Nlz)
    options.param.g_dim_x = ctypes.c_uint32(options.g_dim_x)
    options.param.g_dim_y = ctypes.c_uint32(options.g_dim_y)
    options.param.g_dim_z = ctypes.c_uint32(options.g_dim_z)
    options.param.NiterAD = ctypes.c_uint32(options.NiterAD)
    if isinstance(options.inffi, np.ndarray):
        options.param.inffi = ctypes.c_uint32(options.inffi.item())
    else:
        options.param.inffi = ctypes.c_uint32(options.inffi)
    options.param.Nf = ctypes.c_uint32(options.Nf)
    options.param.deviceNum = ctypes.c_uint32(options.deviceNum)
    options.param.platform = ctypes.c_uint32(options.platform)
    options.param.derivativeType = ctypes.c_uint32(options.derivativeType)
    options.param.TVtype = ctypes.c_uint32(options.TVtype)
    options.param.FluxType = ctypes.c_uint32(options.FluxType)
    options.param.DiffusionType = ctypes.c_uint32(options.DiffusionType)
    options.param.POCS_NgradIter = ctypes.c_uint32(options.POCS_NgradIter)
    options.param.maskFPZ = ctypes.c_uint32(options.maskFPZ)
    options.param.maskBPZ = ctypes.c_uint32(options.maskBPZ)
    options.param.FISTAType = ctypes.c_uint32(options.FISTAType)
    options.param.nProjections = ctypes.c_int64(options.nProjections)
    options.param.TOF_bins = ctypes.c_int64(options.TOF_bins)
    options.param.tau = ctypes.c_float(options.tau)
    options.param.tube_radius = ctypes.c_float(options.tube_radius)
    options.param.epps = ctypes.c_float(options.epps)
    options.param.sigma_x = ctypes.c_float(options.sigma_x)
    options.param.tube_width_z = ctypes.c_float(options.tube_width_z)
    options.param.tube_width_xy = ctypes.c_float(options.tube_width_xy)
    options.param.bmin = ctypes.c_float(options.bmin)
    options.param.bmax = ctypes.c_float(options.bmax)
    options.param.Vmax = ctypes.c_float(options.Vmax)
    options.param.global_factor = ctypes.c_float(options.global_factor)
    options.param.dL = ctypes.c_float(options.dL)
    options.param.flat = ctypes.c_float(options.flat)
    options.param.U = ctypes.c_float(options.U)
    options.param.h_ACOSEM = ctypes.c_float(options.h_ACOSEM)
    options.param.dPitchX = ctypes.c_float(options.dPitchX)
    options.param.dPitchY = ctypes.c_float(options.dPitchY)
    options.param.cr_p = ctypes.c_float(options.cr_p)
    options.param.cr_pz = ctypes.c_float(options.cr_pz)
    options.param.NLMsigma = ctypes.c_float(options.NLMsigma)
    options.param.NLAdaptiveConstant = ctypes.c_float(options.NLAdaptiveConstant)
    options.param.w_sum = ctypes.c_float(options.w_sum)
    options.param.KAD = ctypes.c_float(options.KAD)
    options.param.TimeStepAD = ctypes.c_float(options.TimeStepAD)
    options.param.RDP_gamma = ctypes.c_float(options.RDP_gamma)
    options.param.huber_delta = ctypes.c_float(options.huber_delta)
    options.param.gradV1 = ctypes.c_float(options.gradV1)
    options.param.gradV2 = ctypes.c_float(options.gradV2)
    options.param.alpha0TGV = ctypes.c_float(options.alpha0TGV)
    options.param.alpha1TGV = ctypes.c_float(options.alpha1TGV)
    options.param.GGMRF_p = ctypes.c_float(options.GGMRF_p)
    options.param.GGMRF_q = ctypes.c_float(options.GGMRF_q)
    options.param.GGMRF_c = ctypes.c_float(options.GGMRF_c)
    options.param.beta = ctypes.c_float(options.beta)
    options.param.T = ctypes.c_float(options.B)
    options.param.dSizeXBP = ctypes.c_float(options.dSizeXBP)
    options.param.dSizeZBP = ctypes.c_float(options.dSizeZBP)
    options.param.TVsmoothing = ctypes.c_float(options.TVsmoothing)
    options.param.C = ctypes.c_float(options.C)
    options.param.SATVPhi = ctypes.c_float(options.SATVPhi)
    options.param.eta = ctypes.c_float(options.eta)
    options.param.APLSsmoothing = ctypes.c_float(options.APLSsmoothing)
    options.param.hyperbolicDelta = ctypes.c_float(options.hyperbolicDelta)
    options.param.sourceToCRot = ctypes.c_float(options.sourceToCRot)
    options.param.POCS_alpha = ctypes.c_float(options.POCS_alpha)
    options.param.POCS_rMax = ctypes.c_float(options.POCS_rMax)
    options.param.POCS_alphaRed = ctypes.c_float(options.POCS_alphaRed)
    options.param.POCSepps = ctypes.c_float(options.POCSepps)
    options.param.use_psf = ctypes.c_bool(options.use_psf)
    options.param.TOF = ctypes.c_bool(options.TOF)
    options.param.pitch = ctypes.c_bool(options.pitch)
    options.param.SPECT = ctypes.c_bool(options.SPECT)
    options.param.PET = ctypes.c_bool(options.PET)
    options.param.CT = ctypes.c_bool(options.CT)
    options.param.largeDim = ctypes.c_bool(options.largeDim)
    options.param.loadTOF = ctypes.c_bool(options.loadTOF)
    options.param.storeResidual = ctypes.c_bool(options.storeResidual)
    options.param.FISTA_acceleration = ctypes.c_bool(options.FISTA_acceleration)
    options.param.meanFP = ctypes.c_bool(options.meanFP)
    options.param.meanBP = ctypes.c_bool(options.meanBP)
    options.param.useMaskFP = ctypes.c_bool(options.useMaskFP)
    options.param.useMaskBP = ctypes.c_bool(options.useMaskBP)
    options.param.orthTransaxial = ctypes.c_bool(options.orthTransaxial)
    options.param.orthAxial = ctypes.c_bool(options.orthAxial)
    options.param.enforcePositivity = ctypes.c_bool(options.enforcePositivity)
    options.param.useMultiResolutionVolumes = ctypes.c_bool(options.useMultiResolutionVolumes)
    options.param.save_iter = ctypes.c_bool(options.save_iter)
    options.param.deblurring = ctypes.c_bool(options.deblurring)
    options.param.useMAD = ctypes.c_bool(options.useMAD)
    options.param.useImages = ctypes.c_bool(options.useImages)
    options.param.useEFOV = ctypes.c_bool(options.useEFOV)
    options.param.CTAttenuation = ctypes.c_bool(options.CTAttenuation)
    options.param.offsetCorrection = ctypes.c_bool(options.offsetCorrection)
    options.param.relaxationScaling = ctypes.c_bool(options.relaxationScaling)
    options.param.computeRelaxationParameters = ctypes.c_bool(options.computeRelaxationParameters)
    options.param.storeFP = ctypes.c_bool(options.storeFP)
    options.param.use2DTGV = ctypes.c_bool(options.use2DTGV)
    options.param.med_no_norm = ctypes.c_bool(options.med_no_norm)
    options.param.NLM_MRP = ctypes.c_bool(options.NLM_MRP)
    options.param.NLTV = ctypes.c_bool(options.NLTV)
    options.param.NLRD = ctypes.c_bool(options.NLRD)
    options.param.NLLange = ctypes.c_bool(options.NLLange)
    options.param.NLGGMRF = ctypes.c_bool(options.NLGGMRF)
    options.param.NLM_use_anatomical = ctypes.c_bool(options.NLM_use_anatomical)
    options.param.NLAdaptive = ctypes.c_bool(options.NLAdaptive)
    options.param.TV_use_anatomical = ctypes.c_bool(options.TV_use_anatomical)
    options.param.RDPIncludeCorners = ctypes.c_bool(options.RDPIncludeCorners)
    options.param.RDP_use_anatomical = ctypes.c_bool(options.RDP_use_anatomical)
    options.param.useL2Ball = ctypes.c_bool(options.useL2Ball)
    options.param.saveSens = ctypes.c_bool(options.saveSens)
    options.param.use_64bit_atomics = ctypes.c_bool(options.use_64bit_atomics)
    options.param.use_32bit_atomics = ctypes.c_bool(options.use_32bit_atomics)
    options.param.compute_sensitivity_image = ctypes.c_bool(options.compute_sensitivity_image)
    options.param.useFDKWeights = ctypes.c_bool(options.useFDKWeights)
    options.param.useIndexBasedReconstruction = ctypes.c_bool(options.useIndexBasedReconstruction)
    options.param.stochasticSubsetSelection = ctypes.c_bool(options.stochasticSubsetSelection)
    options.param.useTotLength = ctypes.c_bool(options.useTotLength)
    options.param.OSEM = ctypes.c_bool(options.OSEM)
    options.param.LSQR = ctypes.c_bool(options.LSQR)
    options.param.CGLS = ctypes.c_bool(options.CGLS)
    options.param.SART = ctypes.c_bool(options.SART)
    options.param.FISTA = ctypes.c_bool(options.FISTA)
    options.param.FISTAL1 = ctypes.c_bool(options.FISTAL1)
    options.param.MRAMLA = ctypes.c_bool(options.MRAMLA)
    options.param.RAMLA = ctypes.c_bool(options.RAMLA)
    options.param.ROSEM = ctypes.c_bool(options.ROSEM)
    options.param.RBI = ctypes.c_bool(options.RBI)
    options.param.DRAMA = ctypes.c_bool(options.DRAMA)
    options.param.COSEM = ctypes.c_bool(options.COSEM)
    options.param.ECOSEM = ctypes.c_bool(options.ECOSEM)
    options.param.ACOSEM = ctypes.c_bool(options.ACOSEM)
    options.param.OSL_OSEM = ctypes.c_bool(options.OSL_OSEM)
    options.param.MBSREM = ctypes.c_bool(options.MBSREM)
    options.param.BSREM = ctypes.c_bool(options.BSREM)
    options.param.ROSEM_MAP = ctypes.c_bool(options.ROSEM_MAP)
    options.param.OSL_RBI = ctypes.c_bool(options.OSL_RBI)
    options.param.OSL_COSEM = ctypes.c_bool(options.OSL_COSEM)
    options.param.PKMA = ctypes.c_bool(options.PKMA)
    options.param.SPS = ctypes.c_bool(options.SPS)
    options.param.PDHG = ctypes.c_bool(options.PDHG)
    options.param.PDHGKL = ctypes.c_bool(options.PDHGKL)
    options.param.PDHGL1 = ctypes.c_bool(options.PDHGL1)
    options.param.PDDY = ctypes.c_bool(options.PDDY)
    options.param.CV = ctypes.c_bool(options.CV)
    options.param.POCS = ctypes.c_bool(options.ASD_POCS)
    options.param.FDK = ctypes.c_bool(options.FDK)
    options.param.SAGA = ctypes.c_bool(options.SAGA)
    options.param.MRP = ctypes.c_bool(options.MRP)
    options.param.quad = ctypes.c_bool(options.quad)
    options.param.Huber = ctypes.c_bool(options.Huber)
    options.param.L = ctypes.c_bool(options.L)
    options.param.FMH = ctypes.c_bool(options.FMH)
    options.param.weighted_mean = ctypes.c_bool(options.weighted_mean)
    options.param.TV = ctypes.c_bool(options.TV)
    options.param.hyperbolic = ctypes.c_bool(options.hyperbolic)
    options.param.AD = ctypes.c_bool(options.AD)
    options.param.APLS = ctypes.c_bool(options.APLS)
    options.param.TGV = ctypes.c_bool(options.TGV)
    options.param.NLM = ctypes.c_bool(options.NLM)
    options.param.RDP = ctypes.c_bool(options.RDP)
    options.param.GGMRF = ctypes.c_bool(options.GGMRF)
    options.param.ProxTV = ctypes.c_bool(options.ProxTV)
    options.param.ProxRDP = ctypes.c_bool(options.ProxRDP)
    options.param.ProxNLM = ctypes.c_bool(options.ProxNLM)
    options.param.MAP = ctypes.c_bool(options.MAP)
    options.param.custom = ctypes.c_bool(options.custom)
    options.param.mDim = ctypes.c_uint64(options.SinM.size // options.Nt)
    options.param.nIterSaved = ctypes.c_uint64(options.saveNIter.size)
    options.param.sizeScat = ctypes.c_uint64(options.corrVector.size)
    options.param.eFOV = ctypes.c_uint64(options.eFOVIndices.size)
    if options.listmode and options.compute_sensitivity_image:
        options.param.sizeX = ctypes.c_uint64(options.uV.size)
    else:
        options.param.sizeX = ctypes.c_uint64(options.x.size)
    options.param.sizeZ = ctypes.c_uint64(options.z.size)
    options.param.sizeAtten = ctypes.c_uint64(options.vaimennus.size)
    options.param.sizeNorm = ctypes.c_uint64(options.normalization.size)
    options.param.sizePSF = ctypes.c_uint64(options.gaussK.size)
    options.param.sizeXYind = ctypes.c_uint64(options.xy_index.size)
    options.param.sizeZind = ctypes.c_uint64(options.z_index.size)
    options.param.xCenterSize = ctypes.c_uint64(options.x_center.size)
    options.param.yCenterSize = ctypes.c_uint64(options.y_center.size)
    options.param.zCenterSize = ctypes.c_uint64(options.z_center.size)
    options.param.sizeV = ctypes.c_uint64(options.V.size)
    options.param.measElem = ctypes.c_uint64(options.SinM.size)
    options.param.x = options.x.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.z = options.z.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.uV = options.uV.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dx = options.dx.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dy = options.dy.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dz = options.dz.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.bx = options.bx.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.by = options.by.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.bz = options.bz.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.atten = options.vaimennus.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.norm = options.normalization.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.pituus = options.nMeas.ctypes.data_as(ctypes.POINTER(ctypes.c_int64))
    options.param.xy_index = options.xy_index.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    options.param.z_index = options.z_index.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    options.param.x_center = options.x_center.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.y_center = options.y_center.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.z_center = options.z_center.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.V = options.V.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.gaussPSF = options.gaussK.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.gaussianNLM = options.gaussianNLM.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.saveNiter = options.saveNIter.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    options.param.Nx = options.Nx.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    options.param.Ny = options.Ny.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    options.param.Nz = options.Nz.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    if not options.SinDelayed.dtype == 'single':
        options.SinDelayed = options.SinDelayed.astype(np.float32)
    options.param.randoms = options.SinDelayed.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.corrVector = options.corrVector.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.x0 = options.x0.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.offsetVal = options.OffsetLimit.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dScaleX4 = options.dScaleX4.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dScaleY4 = options.dScaleY4.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dScaleZ4 = options.dScaleZ4.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dSizeX = options.dSizeX.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dSizeY = options.dSizeY.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dScaleX = options.dScaleX.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dScaleY = options.dScaleY.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.dScaleZ = options.dScaleZ.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.kerroin4 = options.kerroin.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.lam_drama = options.lam_drama.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.maskFP = options.maskFP.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    options.param.maskBP = options.maskBP.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    options.param.eFOVIndices = options.eFOVIndices.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    options.param.maskPrior = options.maskPrior.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    options.param.TOFIndices = options.TOFIndices.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    options.param.angles = options.angles.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.blurPlanes = options.blurPlanes.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    options.param.gFilter = options.gFilter.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.gFSize = np.array(options.gFilter.shape, dtype=np.uint64)
    options.param.gFSize = options.gFSize.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64))
    options.param.precondTypeImage = options.precondTypeImage.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    options.param.precondTypeMeas = options.precondTypeMeas.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    options.param.referenceImage = options.referenceImage.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.filterIm = options.filterIm.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.filter = options.filter0.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.filter2 = options.filter2.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.Ffilter = options.Ffilter.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.s = options.s.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.weights_quad = options.weights_quad.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.weights_huber = options.weights_huber.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.weighted_weights = options.weighted_weights.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.APLS_ref_image = options.APLS_ref_image.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.lambdaN = options.lambdaN.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.lambdaFiltered = options.lambdaFiltered.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.alpha_PKMA = options.alpha_PKMA.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.alphaPrecond = options.alphaPrecond.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.tauCP = options.tauCP.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.tauCPFilt = options.tauCPFilt.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.sigmaCP = options.sigmaCP.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.sigma2CP = options.sigma2CP.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.thetaCP = options.thetaCP.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.TOFCenter = options.TOFCenter.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.TV_ref = options.TV_referenceImage.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.trIndices = options.trIndex.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    options.param.axIndices = options.axIndex.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    #For SPECT...
    options.param.crXY = ctypes.c_float(options.crXY)
    options.param.rayShiftsDetector = options.rayShiftsDetector.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.rayShiftsSource = options.rayShiftsSource.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.coneOfResponseStdCoeffA = ctypes.c_float(options.coneOfResponseStdCoeffA)
    options.param.coneOfResponseStdCoeffB = ctypes.c_float(options.coneOfResponseStdCoeffB)
    options.param.coneOfResponseStdCoeffC = ctypes.c_float(options.coneOfResponseStdCoeffC)
    # ...until here
    options.param.NLM_ref = options.NLM_referenceImage.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    options.param.RDP_ref = options.RDP_referenceImage.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    
def reconstructions_mainCT(options):
    options.CT = True
    if options.storeResidual:
        pz, FPOutputP, residual = reconstructions_main(options)
        return pz, FPOutputP, residual
    else:
        pz, FPOutputP = reconstructions_main(options)
        return pz, FPOutputP

def reconstructions_mainSPECT(options):
    options.SPECT = True
    if options.storeResidual:
        pz, FPOutputP, residual = reconstructions_main(options)
        return pz, FPOutputP, residual
    else:
        pz, FPOutputP = reconstructions_main(options)
        return pz, FPOutputP

def reconstructions_main(options):
    import tkinter as tk
    from tkinter.filedialog import askopenfilename
    tic = time.perf_counter()
    options.addProjector()
    print('Preparing for reconstruction...')
    if np.size(options.weights) > 0:
        options.empty_weight = False
    fname, suffix = os.path.splitext(options.fpath)
    if options.SinM.size < 1 and (len(options.fpath) == 0 or len(suffix) == 0):
        root = tk.Tk()
        root.withdraw()
        options.fpath = askopenfilename(title='Select measurement datafile',filetypes=(('NPY, NPZ and MAT files','*.mat *.npy *.npz'),('All','*.*')))
        if len(options.fpath) == 0:
            raise ValueError('No file selected')
    if options.SinM.size < 1 and options.fpath[len(options.fpath)-3:len(options.fpath)+1:1] == 'mat':
        from pymatreader import read_mat
        try:
            var = read_mat(options.fpath)
        except OSError:
            print('File not found, please select the measurement data file')
            root = tk.Tk()
            root.withdraw()
            options.fpath = askopenfilename(title='Select measurement datafile',filetypes=(('MAT files','*.mat'),('All','*.*')))
            if len(options.fpath) == 0:
                raise ValueError('No file selected')
            var = read_mat(options.fpath)
        if options.reconstruct_trues == True:
            options.SinM = np.array(var["SinTrues"],order='F')
        elif options.reconstruct_scatter == True:
            options.SinM = np.array(var["SinScatter"],order='F')
        else:
            if ((options.randoms_correction or options.scatter_correction or options.normalization_correction) and options.corrections_during_reconstruction == False):
                options.SinM = np.array(var["SinM"],order='F')
            else:
                options.SinM = np.array(var["raw_SinM"],order='F')
        if options.randoms_correction and not options.reconstruct_scatter and not options.reconstruct_trues and options.SinDelayed.size < 1:
            options.SinDelayed = np.array(var["SinDelayed"],order='F')
    elif options.SinM.size < 1 and options.fpath[len(options.fpath)-3:len(options.fpath)+1:1] == 'npy':
        options.SinM = np.load(options.fpath)
    elif options.SinM.size < 1 and options.fpath[len(options.fpath)-3:len(options.fpath)+1:1] == 'npz':
        varList = np.load(options.fpath)
        if ((options.randoms_correction or options.scatter_correction or options.normalization_correction) and options.corrections_during_reconstruction == False):
            options.SinM = varList.files[varList.files.index['SinM']]
        else:
            options.SinM = varList.files[varList.files.index['raw_SinM']]
    if options.TOF and options.TOF_bins_used == 1:
        options.TOF_bins = options.TOF_bins_used
        options.SinM = np.sum(options.SinM, axis=3)
        options.TOF = False
    loadCorrections(options)
    # if options.normalization_correction and options.corrections_during_reconstruction == True:
    #     normdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', '..', 'mat-files')) + "/"
    #     if os.path.exists(normdir):
    #         var = read_mat(options.fpath)
    if options.CT and options.flat <= 0:
        print('No flat value input! Using the maximum value as the flat value. Alternatively, input the flat value into options.flat')
        options.flat = np.max(options.SinM).astype(dtype=np.float32)
    if not options.CT and not options.SPECT and not(options.SinM.size == options.Ndist * options.Nang * options.TotSinos) and options.listmode == 0:
        ValueError('The number of elements in the input data does not match the input number of angles, radial distances and total number of sinograms multiplied together!')
    if not options.usingLinearizedData and (options.LSQR or options.CGLS or options.FISTA or options.FISTAL1 or options.PDHG or options.PDHGL1 or options.PDDY or options.FDK or options.SART or options.ASD_POCS) and not options.largeDim and options.CT:
        from .prepass import linearizeData
        linearizeData(options)
        options.usingLinearizedData = True
    if options.listmode == False:
        options.SinM = np.reshape(options.SinM, (int(options.nRowsD), int(options.nColsD), options.nProjections, options.TOF_bins, options.Nt), order='F')
    elif options.listmode == True and options.compute_sensitivity_image:
        from omegatomo.projector.detcoord import getCoordinates
        options.use_raw_data = True
        x, y, z = getCoordinates(options)
        options.use_raw_data = False
        options.uV = np.float32(x)
        options.z = np.float32(z)
    if (options.quad or options.FMH or options.L or options.weighted_mean or options.Huber or options.GGMRF) and options.MAP:
        if hasattr(options, 'weights') and np.size(options.weights) > 0:
            weights_flat = np.array(options.weights).flatten()
            expected_length = ((options.Ndx * 2 + 1) * 
                (options.Ndy * 2 + 1) * 
                (options.Ndz * 2 + 1))
            if len(weights_flat) < expected_length:
                raise ValueError(
                    f'Weights vector is too short, needs to be {expected_length} in length'
                )
            elif len(weights_flat) > expected_length:
                raise ValueError(
                    f'Weights vector is too long, needs to be {expected_length} in length'
                )
            middle_index = int(np.ceil(expected_length / 2)) - 1
            if not np.isinf(weights_flat[middle_index]):
                weights_flat[middle_index] = np.inf
            options.weights = weights_flat
        else:
            options.empty_weight = True
    parseInputs(options, True)
    if not options.CT and (not options.LSQR and not options.CGLS):
        options.SinM[options.SinM < 0] = 0
    if options.FDK:
        options.precondTypeMeas[1] = True
    prepassPhase(options)
    options.tau = 2.5
    if options.use_32bit_atomics and options.use_64bit_atomics:
        options.use_64bit_atomics = False
    if options.use_64bit_atomics and (options.useCPU or options.useCUDA):
        options.use_64bit_atomics = False
    if options.use_32bit_atomics and (options.useCPU or options.useCUDA):
        options.use_32bit_atomics = False
    if options.storeMultiResolution:
        output = np.zeros(int(np.sum(options.N) * options.Nt), dtype=np.float32, order = 'F')
    elif options.useMultiResolutionVolumes:
        output = np.zeros(options.NxOrig * options.NyOrig * options.NzOrig * options.Nt, dtype=np.float32, order = 'F')
    else:
        output = np.zeros(options.Nx[0].item() * options.Ny[0].item() * options.Nz[0].item() * options.Nt, dtype=np.float32, order = 'F')
    if options.storeFP:
        FPOutput = np.zeros(options.SinM.size, dtype=np.float32, order = 'F')
    else:
        FPOutput = np.empty(0, dtype=np.float32)
    if options.storeResidual:
        residual = np.zeros(options.Niter * options.subsets, dtype=np.float32)
    else:
        residual = np.zeros(1, dtype=np.float32)
    options.headerDir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'opencl')) + "/"
    transferData(options)
    inStr = options.headerDir.encode('utf-8')
    libdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..'))
    # point_ptr = ctypes.pointer(options.param)
    if not options.SinM.dtype == 'float32' and not options.largeDim and options.loadTOF:
        options.SinM = options.SinM.astype(np.float32)
    elif not options.SinM.dtype == 'uint16':
        options.SinM = options.SinM.astype(np.float32)
    if options.SinM.ndim > 1:
        options.SinM = options.SinM.ravel('F')
    if options.useCUDA:
        if options.SinM.dtype == 'uint16':
            libN = 'CUDA_matrixfree_uint16_lib'
        else:
            libN = 'CUDA_matrixfree_lib'
        if os.name == 'posix':
            libname = str(os.path.join(libdir, libN + ".so"))
        elif os.name == 'nt':
            libname = str(os.path.join(libdir,libN + ".dll"))
        else:
            libname = str(os.path.join(libdir,libN + ".so"))
    elif options.useCPU:
        if os.name == 'posix':
            libname = str(os.path.join(libdir,"CPU_matrixfree_lib.so"))
        elif os.name == 'nt':
            libname = str(os.path.join(libdir,"CPU_matrixfree_lib.dll"))
        else:
            libname = str(os.path.join(libdir,"CPU_matrixfree_lib.so"))
    else:
        if options.SinM.dtype == 'uint16':
            libN = 'OpenCL_matrixfree_uint16_lib'
        else:
            libN = 'OpenCL_matrixfree_lib'
        if os.name == 'posix':
            libname = str(os.path.join(libdir, libN + ".so"))
        elif os.name == 'nt':
            libname = str(os.path.join(libdir,libN + ".dll"))
        else:
            libname = str(os.path.join(libdir,libN + ".so"))
    residualP = residual.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    if options.SinM.dtype == 'uint16':
        SinoP = options.SinM.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    else:
        SinoP = options.SinM.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    outputP = output.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    FPOutputP = FPOutput.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    c_lib = ctypes.CDLL(libname)
    c_lib.omegaMain(options.param, ctypes.c_char_p(inStr), SinoP, outputP, FPOutputP, residualP)
    if options.useMultiResolutionVolumes:
        output = output.reshape((options.NxOrig, options.NyOrig, options.NzOrig), order = 'F')
    else:
        output = output.reshape((options.Nx[0], options.Ny[0], options.Nz[0]), order = 'F')
    if options.subsets == 1 and options.storeFP == True:
        FPOutput = FPOutput.reshape((options.nRowsD, options.nColsD, options.nProjections, options.TOF_bins), order = 'F')
    toc = time.perf_counter()
    if (options.verbose > 0):
        print(f"Reconstruction took {toc - tic:0.4f} seconds")
    if options.storeResidual:
        return output, FPOutput, residual
    else:
        return output, FPOutput