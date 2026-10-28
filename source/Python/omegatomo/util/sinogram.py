# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 13:37:03 2025

@author: Ville-Veikko Wettenhovi
"""
import numpy as np

def saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, Nang, Ndist, ringDifference, span = 1, rings = 0, TOFbins = 1, segTable = np.empty(0, dtype=np.float32),
                 Nt = 1, detPerRing = 0, cryst_per_block = 0, nDistSide = 1, nPseudos = 0, nLayers = 1, layer1 = np.zeros(0, dtype=np.uint8), layer2 = np.zeros(0, dtype=np.uint8),
                 Sino = np.empty(0, dtype=np.uint16), time = np.zeros(0, dtype=np.uint16), bins = np.zeros(0, dtype=np.uint16), tIndex = np.empty(0, dtype=np.bool_),
                 trueIndex = None, scatterIndex = None, randomIndex = None, SinoT = None, SinoC = None, SinoR = None):
    """
    Bins the events into a sinogram. The events are appended to Sino, if it is non-empty.

    The sinogram is a flat uint16 array with the layout of the MATLAB array (Ndist, Nang, NSinos, nLayers^2, TOFbins, Nt)
    (column-major), i.e. its size is Ndist * Nang * NSinos * nLayers^2 * TOFbins * Nt.

    segTable: cumulative segment table, np.cumsum(options.segment_table). Default: computed from rings.
    time:     time frame of each event, ZERO-BASED (0, ..., Nt - 1). Only used when Nt > 1.
    bins:     TOF bin of each event, 0, ..., TOFbins - 1 (center-first order). Only used when TOFbins > 1.
    layer1, layer2: the layer (0, ..., nLayers - 1) of the first and second detector of each event. Only used when nLayers > 1.
    tIndex:   not used (kept for compatibility)

    Without trueIndex/scatterIndex/randomIndex the return value is the sinogram (Sino). If at least one of the boolean
    per-event arrays trueIndex, scatterIndex or randomIndex is given, the corresponding sinograms are also formed
    (optionally appended to SinoT, SinoC, SinoR, respectively) and the return value is the tuple (Sino, SinoT, SinoC, SinoR),
    where the sinograms of the classes that were not requested are empty arrays. All events are stored in Sino. An event
    that is in several classes is stored only in the first of them (trues, scatter, randoms).
    """

    import os, ctypes
    fPath = os.path.dirname( __file__ )
    if os.path.exists(os.path.join(fPath, 'usingPyPi.py')):
        libdir = os.path.join(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')), "libs")
    else:
        libdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..'))

    if os.name == 'nt':
        libname = str(os.path.join(libdir,"createSinogram.dll"))
    else:
        libname = str(os.path.join(libdir,"createSinogram.so"))
    # Plain Python integers: with NumPy 1.x the product of a NumPy unsigned integer and Python integers becomes uint64/float64,
    # which np.zeros and ctypes do not accept
    Ndist = int(Ndist)
    Nang = int(Nang)
    ringDifference = int(ringDifference)
    span = int(span)
    rings = int(rings)
    TOFbins = int(TOFbins)
    Nt = int(Nt)
    nLayers = int(nLayers)
    nDistSide = int(nDistSide)
    nPseudos = int(nPseudos)
    cryst_per_block = int(cryst_per_block)
    detPerRing = int(detPerRing)
    if segTable.size == 0 and rings == 0 and span > 1:
        raise ValueError('Input either the number of crystal rings (rings) or segment table (segTable)!')
    if span == 1 and rings == 0:
        raise ValueError('If span = 1, the number of crystal rings must be input!')
    if detPerRing == 0:
        detPerRing = Nang * 2
    if Nt < 1 or Nt > 65535:
        raise ValueError('The number of time frames (Nt) must be between 1 and 65535!')
    if nLayers < 1:
        raise ValueError('The number of layers (nLayers) must be at least 1!')

    # The sinogram creation code expects the CUMULATIVE segment table (the same as MATLAB and the other loaders pass, i.e.
    # np.cumsum(options.segment_table)): segTable[(k-1)*2] and segTable[(k-1)*2+1] are the index of the first plane of
    # segment -k and +k, respectively. The last element is thus the total number of planes.
    # The default table is built from the numbers of planes of the segments 0, -1, +1, -2, +2, ... (options.segment_table).
    if span > 1:
        if segTable.size == 0:
            segStop = max(rings*2-1 - 2*ringDifference, rings - ringDifference) - 1
            segSizes = np.concatenate((np.array(rings*2-1, ndmin=1, dtype=np.int64), np.arange(rings*2-1 - (span + 1), segStop, - span*2, dtype=np.int64)))
            segSizes = np.insert(np.repeat(segSizes[1:], 2), 0, segSizes[0])
            segTable = np.cumsum(segSizes)
        else:
            segTable = np.asarray(segTable).ravel()
            if np.any(np.diff(segTable.astype(np.int64)) < 0):
                raise ValueError('segTable must be the cumulative segment table, np.cumsum(options.segment_table), not the numbers of planes of each segment!')
        segTable = np.ascontiguousarray(segTable, dtype=np.uint32)
        NSinos = int(segTable[-1])
    else:
        segTable = np.ascontiguousarray(segTable, dtype=np.uint32)
        NSinos = rings ** 2

    # The memory layout is the same as in MATLAB: (Ndist, Nang, NSinos, nLayers^2, TOFbins, Nt)
    sinoSize = Ndist * Nang * NSinos
    TOFSize = sinoSize * TOFbins
    bufSize = TOFSize * nLayers**2 * Nt

    # The C++ code reads the events as these types
    ring_pos1 = np.ascontiguousarray(ring_pos1, dtype=np.uint16).ravel()
    ring_pos2 = np.ascontiguousarray(ring_pos2, dtype=np.uint16).ravel()
    ring_number1 = np.ascontiguousarray(ring_number1, dtype=np.uint16).ravel()
    ring_number2 = np.ascontiguousarray(ring_number2, dtype=np.uint16).ravel()
    koko = ring_pos1.size
    if ring_pos2.size != koko or ring_number1.size != koko or ring_number2.size != koko:
        raise ValueError('The detector and ring index arrays must have the same number of elements!')

    # Time frame of each event (zero-based), 16-bit unsigned integer in the C++ code. Only used if Nt > 1.
    if Nt > 1:
        timeIn = np.asarray(time).ravel()
        if timeIn.size != koko:
            raise ValueError('If Nt > 1, time must contain the time frame (0, ..., Nt - 1) of each event!')
        if koko > 0 and (timeIn.min() < 0 or timeIn.max() >= Nt):
            raise ValueError('The time frames (time) must be between 0 and Nt - 1!')
        time = np.ascontiguousarray(timeIn, dtype=np.uint16)
    else:
        time = np.zeros(0, dtype=np.uint16)
    # TOF bin of each event. Only used if TOFbins > 1.
    if TOFbins > 1:
        binsIn = np.asarray(bins).ravel()
        if binsIn.size != koko:
            raise ValueError('If TOFbins > 1, bins must contain the TOF bin of each event!')
        if koko > 0 and (binsIn.min() < 0 or binsIn.max() >= TOFbins):
            raise ValueError('The TOF bins (bins) must be between 0 and TOFbins - 1!')
        bins = np.ascontiguousarray(binsIn, dtype=np.uint16)
    else:
        bins = np.zeros(0, dtype=np.uint16)
    # Layers of the detectors of each event. Only used if nLayers > 1.
    if nLayers > 1:
        layer1 = np.ascontiguousarray(layer1, dtype=np.uint8).ravel()
        layer2 = np.ascontiguousarray(layer2, dtype=np.uint8).ravel()
        if layer1.size != koko or layer2.size != koko:
            raise ValueError('If nLayers > 1, layer1 and layer2 must contain the layer of each event!')
        if koko > 0 and (layer1.max() >= nLayers or layer2.max() >= nLayers):
            raise ValueError('The layers (layer1, layer2) must be between 0 and nLayers - 1!')
    else:
        layer1 = np.zeros(0, dtype=np.uint8)
        layer2 = np.zeros(0, dtype=np.uint8)

    def getSino(S, name):
        # An empty array is replaced with zeros, a non-empty one is appended to
        if S is None or S.size == 0:
            return np.zeros(bufSize, dtype=np.uint16)
        if S.dtype != np.uint16 or not S.flags.c_contiguous or S.size != bufSize:
            raise ValueError(name + ' must be a C-contiguous uint16 array with ' + str(bufSize) + ' elements (Ndist * Nang * NSinos * nLayers^2 * TOFbins * Nt)!')
        return S

    Sino = getSino(Sino, 'Sino')
    # Classes of the events (the class sinograms are formed only if requested)
    sIndex = np.empty(0, dtype=np.bool_)
    classes = []
    for idx, name in ((trueIndex, 'trueIndex'), (scatterIndex, 'scatterIndex'), (randomIndex, 'randomIndex')):
        if idx is not None:
            idx = np.ascontiguousarray(idx, dtype=np.bool_).ravel()
            if idx.size != koko:
                raise ValueError(name + ' must have one element per event!')
        classes.append(idx)
    classSinos = (SinoT, SinoC, SinoR)
    classNames = ('SinoT', 'SinoC', 'SinoR')
    sinos = []
    for k in range(3):
        if classes[k] is not None:
            sinos.append(getSino(classSinos[k], classNames[k]))
        else:
            sinos.append(np.empty(0, dtype=np.uint16))
    SinoT, SinoC, SinoR = sinos
    trueIdx = classes[0] if classes[0] is not None else sIndex
    scatterIdx = classes[1] if classes[1] is not None else sIndex
    randIdx = classes[2] if classes[2] is not None else sIndex

    SinoP = Sino.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    SinoTP = SinoT.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    SinoCP = SinoC.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    SinoRP = SinoR.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    ring1P = ring_number1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    ring2P = ring_number2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    pos1P = ring_pos1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    pos2P = ring_pos2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    truesP = trueIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    scatterP = scatterIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    randP = randIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    segP = segTable.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    timeP = time.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    binsP = bins.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    layer1P = layer1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    layer2P = layer2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    # The sizes of the class sinograms are only flags in the C++ code (> 1 means that the sinogram is formed)
    tSize = 2 if classes[0] is not None else 0
    sSize = 2 if classes[1] is not None else 0
    rSize = 2 if classes[2] is not None else 0
    from omegatomo.util.dllpath import addDLLDirectories
    addDLLDirectories()
    c_lib = ctypes.CDLL(libname)
    c_lib.sinoMain(pos1P, pos2P, ring1P, ring2P, truesP, scatterP, randP, ctypes.c_uint64(sinoSize), ctypes.c_uint32(Ndist),
                    ctypes.c_uint32(Nang), ctypes.c_uint32(ringDifference), ctypes.c_uint32(span), segP, ctypes.c_uint64(sinoSize),
                    ctypes.c_uint64(TOFSize), timeP, ctypes.c_uint64(Nt), ctypes.c_int32(detPerRing), ctypes.c_int32(rings),
                    binsP, ctypes.c_int32(nDistSide), ctypes.c_int32(detPerRing), ctypes.c_int32(nPseudos), ctypes.c_int32(cryst_per_block),
                    ctypes.c_int32(nLayers), layer1P, layer2P, ctypes.c_int64(koko), ctypes.c_int64(tSize), ctypes.c_int64(sSize),
                    ctypes.c_int64(rSize), SinoP, SinoTP, SinoCP, SinoRP)

    if tSize > 0 or sSize > 0 or rSize > 0:
        return Sino, SinoT, SinoC, SinoR
    return Sino