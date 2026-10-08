# -*- coding: utf-8 -*-


def _toArray(value, dtype, nMin = 1):
    # Contiguous 1D array of the given type. Repeated up to nMin elements if needed (per-layer parameters)
    import numpy as np
    arr = np.ascontiguousarray(np.atleast_1d(np.asarray(value)).ravel().astype(dtype))
    if arr.size < nMin:
        arr = np.ascontiguousarray(np.repeat(arr[:1], nMin))
    return arr


def _scalar(value):
    # First element of a scalar or array-like parameter as a Python number
    import numpy as np
    return np.asarray(value).ravel()[0].item()


def _ptr(arr, ctype):
    import ctypes
    return arr.ctypes.data_as(ctypes.POINTER(ctype))


def _loadRootLibrary():
    import os
    import ctypes
    fPath = os.path.dirname( __file__ )
    if os.path.exists(os.path.join(fPath, '..', 'util', 'usingPyPi.py')):
        libdir = os.path.join(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')), "libs")
    else:
        libdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..'))
    if os.name == 'nt':
        libname = str(os.path.join(libdir,"libRoot.dll"))
    else:
        libname = str(os.path.join(libdir,"libRoot.so"))
    from omegatomo.util.dllpath import addDLLDirectories
    addDLLDirectories()
    c_lib = ctypes.CDLL(libname)
    u16 = ctypes.POINTER(ctypes.c_uint16)
    u32 = ctypes.POINTER(ctypes.c_uint32)
    c_lib.rootEntries.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)]
    c_lib.rootEntries.restype = ctypes.c_int
    # Must match rootMain in libRoot.h exactly
    c_lib.rootMain.argtypes = [
        ctypes.c_char_p, ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_bool, ctypes.c_uint32, # rootFile, tPoints, alku, loppu, source, linear_multip
        u32, ctypes.c_uint32, u32, # cryst_per_block, blocks_per_ring, det_per_ring
        u16, u16, u16, u16, u16, u16, u16, # S, SC, RA, trIndex, axIndex, DtrIndex, DaxIndex
        ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.POINTER(ctypes.c_uint8), # obtain_trues, store_scatter, store_randoms, scatter_components
        ctypes.c_bool, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_bool, # randoms_correction, coord, Dcoord, store_coordinates
        u32, ctypes.c_uint32, u32, ctypes.POINTER(ctypes.c_uint64), ctypes.c_uint32, u32, ctypes.c_uint32, ctypes.c_uint32, # cryst_per_block_z, transaxial_multip, rings, sinoSize, Ndist, Nang, ringDifference, span
        u32, ctypes.c_int64, ctypes.c_uint64, ctypes.c_int32, # seg, Nt, TOFSize, nDistSide
        u16, u16, u16, u16, u16, # Sino, SinoT, SinoC, SinoR, SinoD
        u32, ctypes.c_int32, ctypes.c_double, ctypes.c_double, ctypes.c_bool, ctypes.c_int32, # detWPseudo, nPseudos, binSize, FWHM, verbose, nLayers
        ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, # dx, dy, dz, bx, by, bz
        ctypes.c_int64, ctypes.c_int64, ctypes.c_int64, ctypes.c_bool, ctypes.c_bool, # Nx, Ny, Nz, dualLayerSubmodule, indexBased
        u16, ctypes.POINTER(ctypes.c_uint8)] # tIndex, TOFIndex
    c_lib.rootMain.restype = ctypes.c_int
    return c_lib


def _findRootFiles(fpath):
    import os
    import glob
    if os.path.isfile(fpath):
        return [fpath]
    if os.path.isdir(fpath):
        return sorted(glob.glob(os.path.join(fpath, '*.root')))
    return []


def loadROOT(options, store_coordinates = False):
    """
    Loads GATE ROOT data.

    Returns
    -------
    Sino, SinoT, SinoC, SinoR, SinoD, Fcoord, FDcoord, DtrIndex, DaxIndex, TOFIndices
    Sinogram mode: the sinograms (Fcoord and FDcoord are empty).
    Coordinate mode (store_coordinates = True): Fcoord is the (6, N) coordinate array (list of Nt arrays in dynamic mode).
    Index-based mode (options.useIndexBasedReconstruction): Fcoord/FDcoord are the (2, N) transaxial/axial detector
    indices (lists in dynamic mode), DtrIndex/DaxIndex the same for the delayed coincidences (static only).
    TOFIndices is the (N,) TOF bin index of each event (list in dynamic mode) when TOF is used with coordinates or indices.
    """
    import os
    import ctypes
    import numpy as np
    import math
    from omegatomo.projector import computePixelSize

    indexBased = bool(options.useIndexBasedReconstruction)
    if indexBased:
        store_coordinates = False
    store_coordinates = bool(store_coordinates)
    nLayers = int(_scalar(options.nLayers))
    ringsArr = _toArray(options.rings, np.uint32, nLayers)
    NangArr = _toArray(options.Nang, np.uint32)
    Ndist = int(_scalar(options.Ndist))
    totSinos = int(_scalar(options.TotSinos))
    if options.span == 1:
        totSinos = int(ringsArr[0]) ** 2
    sinoSize = Ndist * int(NangArr[0]) * totSinos
    TOF_bins = int(_scalar(options.TOF_bins))
    TOFSize = sinoSize * TOF_bins

    xx, yy, zz = computePixelSize(options)

    if isinstance(options.Nx, np.ndarray):
        Nx = options.Nx[0].item()
        Ny = options.Ny[0].item()
        Nz = options.Nz[0].item()
    else:
        Nx = options.Nx
        Ny = options.Ny
        Nz = options.Nz
    dx = options.dx[0].item()
    dy = options.dy[0].item()
    dz = options.dz[0].item()
    bx = options.bx[0].item()
    by = options.by[0].item()
    bz = options.bz[0].item()

    TOF = TOF_bins > 1


    if TOF:
        FWHM = options.TOF_FWHM / (2. * math.sqrt(2. * math.log(2.)))
    else:
        FWHM = 0.

    if isinstance(options.pseudot, np.ndarray):
        if options.pseudot.size == 0:
            nPseudos = 0
        elif options.pseudot.size == 1:
            nPseudos = options.pseudot[0].item()
        else:
            nPseudos = options.pseudot.size
    else:
        nPseudos = options.pseudot

    alku = options.start
    loppu = options.end
    if np.isinf(loppu):
        loppu = 1e9
    if isinstance(options.partitions, np.ndarray):
        if options.partitions.size > 1:
            Nt = options.partitions.size
        else:
            vali = (loppu - alku) / options.partitions[0].item()
            options.partitions = np.repeat(vali, options.partitions)
            Nt = options.partitions.size
    else:
        Nt = options.partitions
        vali = (loppu - alku) / options.partitions
        options.partitions = np.repeat(vali, options.partitions)
    Nt = int(Nt)
    dynamic = Nt > 1
    # Time-step durations
    tPoints = np.ascontiguousarray(np.asarray(options.partitions, dtype=np.float64).ravel())
    if options.source:
        C = np.zeros((Nx, Ny, Nz, Nt), dtype=np.uint16, order='F')
        if options.store_scatter:
            SC = np.zeros((Nx, Ny, Nz, Nt), dtype=np.uint16, order='F')
        else:
            SC = np.zeros((1, 1, 1), dtype=np.uint16, order='F')
        if options.store_randoms:
            RA = np.zeros((Nx, Ny, Nz, Nt), dtype=np.uint16, order='F')
        else:
            RA = np.zeros((1, 1, 1), dtype=np.uint16, order='F')
    else:
        RA = np.zeros((1, 1, 1), dtype=np.uint16, order='F')
        SC = np.zeros((1, 1, 1), dtype=np.uint16, order='F')
        C = np.zeros((1, 1, 1), dtype=np.uint16, order='F')
    if indexBased:
        # The sinograms are not formed when saving the detector indices (they are when saving coordinates)
        Sino = np.zeros(1, dtype=np.uint16, order='F')
        SinoT = np.zeros(1, dtype=np.uint16, order='F')
        SinoC = np.zeros(1, dtype=np.uint16, order='F')
        SinoR = np.zeros(1, dtype=np.uint16, order='F')
        SinoD = np.zeros(1, dtype=np.uint16, order='F')
    else:
        # The memory layout matches the index computed in saveSinogram.h
        Sino = np.zeros((Ndist, int(NangArr[0]), totSinos, nLayers**2, TOF_bins, Nt), dtype=np.uint16, order='F')
        if options.obtain_trues:
            SinoT = np.zeros((Ndist, int(NangArr[0]), totSinos, nLayers**2, TOF_bins, Nt), dtype=np.uint16, order='F')
        else:
            SinoT = np.zeros(1, dtype=np.uint16, order='F')
        if options.store_scatter:
            SinoC = np.zeros((Ndist, int(NangArr[0]), totSinos, nLayers**2, TOF_bins, Nt), dtype=np.uint16, order='F')
        else:
            SinoC = np.zeros(1, dtype=np.uint16, order='F')
        if options.store_randoms:
            SinoR = np.zeros((Ndist, int(NangArr[0]), totSinos, nLayers**2, TOF_bins, Nt), dtype=np.uint16, order='F')
        else:
            SinoR = np.zeros(1, dtype=np.uint16, order='F')
        if options.randoms_correction:
            SinoD = np.zeros((Ndist, int(NangArr[0]), totSinos, nLayers**2, Nt), dtype=np.uint16, order='F')
        else:
            SinoD = np.zeros(1, dtype=np.uint16, order='F')
    seg = np.ascontiguousarray(np.cumsum(options.segment_table), dtype=np.uint32)

    c_lib = _loadRootLibrary()

    files = _findRootFiles(options.fpath)
    if len(files) == 0:
        print('No files found! Please select a ROOT file')
        import tkinter as tk
        from tkinter.filedialog import askopenfilename
        import glob
        root = tk.Tk()
        root.withdraw()
        filename = askopenfilename(title='Select first ROOT file',filetypes=([('ROOT Files','*.root')]))
        if not filename:
            raise ValueError('No file was selected')
        files = sorted(glob.glob(os.path.join(os.path.split(filename)[0], '*.root')))
    nFiles = len(files)
    if nFiles == 0:
        raise ValueError('No ROOT files found!')

    # Parameters that are the same for every file
    cryst_per_block = _toArray(options.cryst_per_block, np.uint32, nLayers)
    det_per_ring = _toArray(options.det_per_ring, np.uint32, nLayers)
    cryst_per_block_axial = _toArray(options.cryst_per_block_axial, np.uint32, nLayers)
    rings = _toArray(options.rings, np.uint32, nLayers)
    Nang = _toArray(options.Nang, np.uint32, nLayers)
    det_w_pseudo = _toArray(options.det_w_pseudo, np.uint32, nLayers)
    sinoSizeArr = np.array([sinoSize], dtype=np.uint64)
    scatter_components = np.zeros(4, dtype=np.uint8)
    sc = np.atleast_1d(np.asarray(options.scatter_components)).ravel()
    scatter_components[:min(4, sc.size)] = sc[:4].astype(np.uint8)
    useTOFInd = TOF and (indexBased or store_coordinates)
    delays = bool(options.randoms_correction)
    if dynamic and delays and (indexBased or store_coordinates):
        if store_coordinates:
            print('Coordinates for delayed coincidences are not supported for list-mode data in dynamic mode!')
        else:
            print('Indices for delayed coincidences are not supported for list-mode data in dynamic mode!')

    # Per-file results
    fTr, fAx, fTOF, fCoord = [], [], [], []
    fDTr, fDAx, fDCoord = [], [], []
    fTime = []

    for lk in range(0, nFiles):

        rootFile = files[lk]

        nCoinc = ctypes.c_int64(0)
        nDelay = ctypes.c_int64(0)
        if c_lib.rootEntries(rootFile.encode('utf-8'), ctypes.byref(nCoinc), ctypes.byref(nDelay)) != 0:
            raise ValueError('Unable to open ROOT file ' + rootFile)
        Nentries = nCoinc.value
        Ndelays = nDelay.value if delays else 0

        if indexBased:
            trIndex = np.zeros(2 * Nentries, dtype=np.uint16)
            axIndex = np.zeros(2 * Nentries, dtype=np.uint16)
            if delays:
                DtrIndex = np.zeros(2 * Ndelays, dtype=np.uint16)
                DaxIndex = np.zeros(2 * Ndelays, dtype=np.uint16)
            else:
                DtrIndex = np.zeros(1, dtype=np.uint16)
                DaxIndex = np.zeros(1, dtype=np.uint16)
        else:
            trIndex = np.zeros(1, dtype=np.uint16)
            axIndex = np.zeros(1, dtype=np.uint16)
            DtrIndex = np.zeros(1, dtype=np.uint16)
            DaxIndex = np.zeros(1, dtype=np.uint16)
        if store_coordinates:
            coord = np.zeros(6 * Nentries, dtype=np.float32)
            if delays:
                Dcoord = np.zeros(6 * Ndelays, dtype=np.float32)
            else:
                Dcoord = np.zeros(1, dtype=np.float32)
        else:
            coord = np.zeros(1, dtype=np.float32)
            Dcoord = np.zeros(1, dtype=np.float32)
        if useTOFInd:
            TOFIndex = np.zeros(Nentries, dtype=np.uint8)
        else:
            TOFIndex = np.zeros(1, dtype=np.uint8)
        if dynamic and (indexBased or store_coordinates):
            tIndex = np.full(Nentries, 32768, dtype=np.uint16)
        else:
            tIndex = np.zeros(1, dtype=np.uint16)

        c_lib.rootMain(rootFile.encode('utf-8'), _ptr(tPoints, ctypes.c_double), alku, loppu, bool(options.source), int(options.linear_multip), _ptr(cryst_per_block, ctypes.c_uint32),
                       int(options.blocks_per_ring), _ptr(det_per_ring, ctypes.c_uint32), _ptr(C, ctypes.c_uint16), _ptr(SC, ctypes.c_uint16), _ptr(RA, ctypes.c_uint16),
                       _ptr(trIndex, ctypes.c_uint16), _ptr(axIndex, ctypes.c_uint16), _ptr(DtrIndex, ctypes.c_uint16), _ptr(DaxIndex, ctypes.c_uint16),
                       bool(options.obtain_trues), bool(options.store_scatter), bool(options.store_randoms), _ptr(scatter_components, ctypes.c_uint8),
                       delays, _ptr(coord, ctypes.c_float), _ptr(Dcoord, ctypes.c_float), store_coordinates, _ptr(cryst_per_block_axial, ctypes.c_uint32),
                       int(options.transaxial_multip), _ptr(rings, ctypes.c_uint32), _ptr(sinoSizeArr, ctypes.c_uint64), Ndist, _ptr(Nang, ctypes.c_uint32),
                       int(options.ring_difference), int(options.span), _ptr(seg, ctypes.c_uint32), Nt, TOFSize, int(options.ndist_side),
                       _ptr(Sino, ctypes.c_uint16), _ptr(SinoT, ctypes.c_uint16), _ptr(SinoC, ctypes.c_uint16), _ptr(SinoR, ctypes.c_uint16), _ptr(SinoD, ctypes.c_uint16),
                       _ptr(det_w_pseudo, ctypes.c_uint32), int(nPseudos), float(options.TOF_width), float(FWHM), bool(options.verbose), nLayers, dx, dy, dz, bx, by, bz,
                       int(Nx), int(Ny), int(Nz), bool(options.dualLayerSubmodule), indexBased, _ptr(tIndex, ctypes.c_uint16), _ptr(TOFIndex, ctypes.c_uint8))

        # Remove the events that were skipped (outside the time window or TOF range), they are left as zeros
        if indexBased:
            tr = trIndex.reshape((-1, 2))
            ax = axIndex.reshape((-1, 2))
            keep = ~((tr[:, 0] == 0) & (tr[:, 1] == 0) & (ax[:, 0] == 0) & (ax[:, 1] == 0))
            fTr.append(tr[keep, :])
            fAx.append(ax[keep, :])
            if delays and not dynamic:
                Dtr = DtrIndex.reshape((-1, 2))
                Dax = DaxIndex.reshape((-1, 2))
                keepD = ~((Dtr[:, 0] == 0) & (Dtr[:, 1] == 0) & (Dax[:, 0] == 0) & (Dax[:, 1] == 0))
                fDTr.append(Dtr[keepD, :])
                fDAx.append(Dax[keepD, :])
        elif store_coordinates:
            cc = coord.reshape((-1, 6))
            keep = np.any(cc != 0, axis=1)
            fCoord.append(cc[keep, :])
            if delays and not dynamic:
                Dcc = Dcoord.reshape((-1, 6))
                fDCoord.append(Dcc[np.any(Dcc != 0, axis=1), :])
        if useTOFInd:
            fTOF.append(TOFIndex[keep])
        if dynamic and (indexBased or store_coordinates):
            fTime.append(tIndex[keep])
        print('File ' + rootFile + ' loaded')

    def _cat(parts, cols, dtype):
        if len(parts) == 0:
            return np.zeros((0, cols), dtype=dtype)
        return np.concatenate(parts, axis=0)

    Fcoord = np.empty(0, dtype=np.float32)
    FDcoord = np.empty(0, dtype=np.float32)
    DtrIndex = np.empty(0, dtype=np.uint16)
    DaxIndex = np.empty(0, dtype=np.uint16)
    TOFIndices = np.empty(0, dtype=np.uint8)
    if indexBased or store_coordinates:
        if indexBased:
            first = _cat(fTr, 2, np.uint16)
            second = _cat(fAx, 2, np.uint16)
        else:
            first = _cat(fCoord, 6, np.float32)
            second = None
        TOFall = np.concatenate(fTOF) if useTOFInd and len(fTOF) > 0 else np.empty(0, dtype=np.uint8)
        if not dynamic:
            Fcoord = np.asfortranarray(first.T)
            if indexBased:
                FDcoord = np.asfortranarray(second.T)
                if delays:
                    DtrIndex = np.asfortranarray(_cat(fDTr, 2, np.uint16).T)
                    DaxIndex = np.asfortranarray(_cat(fDAx, 2, np.uint16).T)
            elif delays:
                FDcoord = np.asfortranarray(_cat(fDCoord, 6, np.float32).T)
            TOFIndices = TOFall
        else:
            timeAll = np.concatenate(fTime) if len(fTime) > 0 else np.empty(0, dtype=np.uint16)
            Fcoord = list()
            if indexBased:
                FDcoord = list()
            if useTOFInd:
                TOFIndices = list()
            for uu in range(0, Nt):
                ind = timeAll == uu
                Fcoord.append(np.asfortranarray(first[ind, :].T))
                if indexBased:
                    FDcoord.append(np.asfortranarray(second[ind, :].T))
                if useTOFInd:
                    TOFIndices.append(TOFall[ind])
    return Sino, SinoT, SinoC, SinoR, SinoD, Fcoord, FDcoord, DtrIndex, DaxIndex, TOFIndices
