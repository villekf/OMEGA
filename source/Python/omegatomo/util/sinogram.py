# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 13:37:03 2025

@author: Ville-Veikko Wettenhovi
"""
import numpy as np

def saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, Nang, Ndist, ringDifference, span = 1, rings = 0, TOFbins = 1, segTable = np.empty(0, dtype=np.float32), 
                 Nt = 1, detPerRing = 0, cryst_per_block = 0, nDistSide = 1, nPseudos = 0, nLayers = 1, layer1 = np.zeros(0, dtype=np.uint8), layer2 = np.zeros(0, dtype=np.uint8), 
                 Sino = np.empty(0, dtype=np.uint16), time = np.zeros(0, dtype=np.uint16), bins = np.zeros(0, dtype=np.uint16), tIndex = np.empty(0, dtype=np.bool_)):
    
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
    if segTable.size == 0 and rings == 0 and span > 1:
        raise ValueError('Input either the number of crystal rings (rings) or segment table (segTable)!')
    if span == 1 and rings == 0:
        raise ValueError('If span = 1, the number of crystal rings must be input!')
    if detPerRing == 0:
        detPerRing = Nang * 2
        
    if segTable.size == 0 and span > 0 and rings > 0:
        segTable = np.concatenate((np.array(rings*2-1,ndmin=1, dtype=np.uint32), np.arange(rings*2-1 - (span + 1), rings - ringDifference, - span*2, dtype=np.uint32)))
        segTable = np.insert(np.repeat(segTable[1:], 2), 0, segTable[0])
        
    if span == 1:
        NSinos = rings ** 2
    else:
        NSinos = np.sum(segTable)
        
    sinoSize = Ndist * Nang * NSinos
    TOFSize = Ndist * Nang * NSinos * TOFbins
    koko = ring_pos1.size
    if Sino.size == 0:
        Sino = np.zeros(TOFSize, dtype=np.uint16)
    if tIndex.size > 0:
        SinoT = np.zeros(sinoSize, dtype=np.uint16)
    else:
        SinoT = np.empty(0, dtype=np.uint16)
    
    sIndex = np.empty(0, dtype=np.bool_)
    rIndex = np.empty(0, dtype=np.bool_)
    SinoC = np.empty(0, dtype=np.uint16)
    SinoR = np.empty(0, dtype=np.uint16)
    SinoP = Sino.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    SinoTP = SinoT.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    SinoCP = SinoC.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    SinoRP = SinoR.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    ring1P = ring_number1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    ring2P = ring_number2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    pos1P = ring_pos1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    pos2P = ring_pos2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    truesP = sIndex.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    scatterP = sIndex.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    randP = rIndex.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    segP = segTable.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    timeP = time.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    binsP = bins.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    layer1P = layer1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    layer2P = layer2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
    c_lib = ctypes.CDLL(libname)
    c_lib.sinoMain(pos1P, pos2P, ring1P, ring2P, truesP, scatterP, randP, ctypes.c_uint64(sinoSize), ctypes.c_uint32(Ndist), 
                    ctypes.c_uint32(Nang), ctypes.c_uint32(ringDifference), ctypes.c_uint32(span), segP, ctypes.c_uint64(sinoSize), 
                    ctypes.c_uint64(TOFSize), timeP, ctypes.c_uint64(Nt), ctypes.c_int32(detPerRing), ctypes.c_int32(rings), 
                    binsP, ctypes.c_int32(nDistSide), ctypes.c_int32(detPerRing), ctypes.c_int32(nPseudos), ctypes.c_int32(cryst_per_block),
                    ctypes.c_int32(nLayers), layer1P, layer2P, ctypes.c_int64(koko), ctypes.c_int64(0), ctypes.c_int64(0), 
                    ctypes.c_int64(0), SinoP, SinoTP, SinoCP, SinoRP)
    
    return Sino