# -*- coding: utf-8 -*-



def loadROOT(options, store_coordinates = False):
    import ROOT
    import os
    import ctypes
    import numpy as np
    import math
    import glob
    from omegatomo.projector import computePixelSize
    
    totSinos = options.TotSinos;
    if options.span == 1:
        totSinos = options.rings**2
    sinoSize = options.Ndist * options.Nang * totSinos
    TOFSize = sinoSize * options.TOF_bins
    
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
    
    TOF = options.TOF_bins > 1
    
    
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
    if isinstance(options.cryst_per_block, np.ndarray):
        cryst_per_block = options.cryst_per_block[0].item()
    else:
        cryst_per_block = options.cryst_per_block
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
    Sino = np.zeros((options.Ndist, options.Nang, totSinos, options.nLayers**2, options.TOF_bins, Nt), dtype=np.uint16, order='F')
    if options.obtain_trues:
        SinoT = np.zeros((options.Ndist, options.Nang, totSinos, options.nLayers**2, options.TOF_bins, Nt), dtype=np.uint16, order='F')
    else:
        SinoT = np.zeros(1, dtype=np.uint16, order='F')
    if options.store_scatter:
        SinoC = np.zeros((options.Ndist, options.Nang, totSinos, options.nLayers**2, options.TOF_bins, Nt), dtype=np.uint16, order='F')
    else:
        SinoC = np.zeros(1, dtype=np.uint16, order='F')
    if options.store_randoms:
        SinoR = np.zeros((options.Ndist, options.Nang, totSinos, options.nLayers**2, options.TOF_bins, Nt), dtype=np.uint16, order='F')
    else:
        SinoR = np.zeros(1, dtype=np.uint16, order='F')
    if options.randoms_correction:
        SinoD = np.zeros((options.Ndist, options.Nang, totSinos, options.nLayers**2, Nt), dtype=np.uint16, order='F')
    else:
        SinoD = np.zeros(1, dtype=np.uint16, order='F')
    seg = np.cumsum(options.segment_table)
    
    fPath = os.path.dirname( __file__ )
    if os.path.exists(os.path.join(fPath, '..', 'util', 'usingPyPi.py')):
        libdir = os.path.join(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')), "libs")
    else:
        libdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..'))
    # dir2 = os.path.join(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..','..')), 'cpp')
    if os.name == 'posix':
        libname = str(os.path.join(libdir,"libRoot.so"))
    elif os.name == 'nt':
        libname = str(os.path.join(libdir,"libRoot.dll"))
    else:
        libname = str(os.path.join(libdir,"libRoot.so"))
    c_lib = ctypes.CDLL(libname)
    DtrIndex = np.empty(0, dtype=np.uint16)
    DaxIndex = np.empty(0, dtype=np.uint16)
    if store_coordinates:
        Fcoord = np.empty(0, dtype=np.float32)
        FDcoord = np.empty(0, dtype=np.float32)
    else:
        Fcoord = np.empty(0, dtype=np.uint16)
        FDcoord = np.empty(0, dtype=np.uint16)
    if store_coordinates and Nt > 1:
        Fcoord = list()
        FDcoord = list()
        for uu in range(0, Nt):
            Fcoord.append(np.empty(0, dtype=np.float32))
            if options.randoms_correction:
                FDcoord.append(np.empty(0, dtype=np.float32))
    
    filename, file_extension = os.path.splitext(options.fpath)
    if file_extension == 'root':
        nFiles = 1
        rootFile = list([options.fpath])
    else:
        files = glob.glob(os.path.join(options.fpath, '*.root'))
        nFiles = len(files)
    if nFiles == 0:
        print('No files found! Please select a ROOT file')
        import tkinter as tk
        from tkinter.filedialog import askopenfilename
        root = tk.Tk()
        root.withdraw()
        filename = askopenfilename(title='Select first ROOT file',filetypes=([('ROOT Files','*.root')]))
        if not filename:
            raise ValueError('No file was selected')
        files = glob.glob(os.path.join(os.path.split(filename)[0], '*.root'))
        nFiles = len(files)
        
    
    
    for lk in range(0, nFiles):
        
        rootFile = files[lk]
    
        file = ROOT.TFile.Open(rootFile)
        
        joku = file.Coincidences
        
        Nentries = joku.GetEntries()
        
        if options.randoms_correction:
            joku = file.delay
            Dentries = joku.GetEntries()
            
        file.Close()
        
        inStr = rootFile.encode('utf-8')
        if store_coordinates:
            tPoints = np.zeros(Nentries, dtype=np.uint16)
        else:
            tPoints = np.zeros(1, dtype=np.uint16)
            if options.useIndexBasedReconstruction:
                trIndices = np.zeros(Nentries, dtype=np.uint16)
                axIndices = np.zeros(Nentries, dtype=np.uint16)
                if options.randoms_correction:
                    DtrIndices = np.zeros(Dentries, dtype=np.uint16)
                    DaxIndices = np.zeros(Dentries, dtype=np.uint16)
                else:
                    DtrIndices = np.zeros(1, dtype=np.uint16)
                    DaxIndices = np.zeros(1, dtype=np.uint16)
            else:
                trIndices = np.zeros(1, dtype=np.uint16)
                axIndices = np.zeros(1, dtype=np.uint16)
                DtrIndices = np.zeros(1, dtype=np.uint16)
                DaxIndices = np.zeros(1, dtype=np.uint16)
        tPointP = tPoints.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        segP = seg.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        CP = C.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        RAP = RA.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        SCP = SC.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        SinoP = Sino.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        SinoTP = SinoT.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        SinoCP = SinoC.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        SinoRP = SinoR.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        SinoDP = SinoD.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        trIndicesP = trIndices.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        axIndicesP = axIndices.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        DtrIndicesP = DtrIndices.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        DaxIndicesP = DaxIndices.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        scatterComP = options.scatter_components.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
        partitionsP = options.partitions.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        if not isinstance(cryst_per_block, np.ndarray):
            cryst_per_block = np.array(cryst_per_block,dtype=np.uint32)
        cryst_per_blockP = cryst_per_block.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        if not isinstance(options.det_per_ring, np.ndarray):
            det_per_ring = np.array(options.det_per_ring,dtype=np.uint32)
        else:
            det_per_ring = np.uint32(options.det_per_ring)
        det_per_ringP = det_per_ring.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        if not isinstance(options.cryst_per_block_axial, np.ndarray):
            cryst_per_block_axial = np.array(options.cryst_per_block_axial,dtype=np.uint32)
        else:
            cryst_per_block_axial = np.uint32(options.cryst_per_block_axial)
        cryst_per_block_axialP = cryst_per_block_axial.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        if not isinstance(options.rings, np.ndarray):
            rings = np.array(options.rings,dtype=np.uint32)
        else:
            rings = np.uint32(options.rings)
        ringsP = rings.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        if not isinstance(sinoSize, np.ndarray):
            sinoSize = np.array(sinoSize,dtype=np.uint64)
        sinoSizeP = sinoSize.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64))
        if not isinstance(options.Nang, np.ndarray):
            Nang = np.array(options.Nang,dtype=np.uint32)
        else:
            Nang = np.uint32(options.Nang)
        NangP = Nang.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        if not isinstance(options.det_w_pseudo, np.ndarray):
            det_w_pseudo = np.array(options.det_w_pseudo,dtype=np.uint32)
        else:
            det_w_pseudo = np.uint32(options.det_w_pseudo)
        det_w_pseudoP = det_w_pseudo.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        if store_coordinates:
            coord = np.zeros(Nentries, dtype=np.float32)
            if options.randoms_correction:
                Rcoord = np.zeros(Dentries, dtype=np.float32)
            else:
                Rcoord = np.zeros(1, dtype=np.float32)
        else:
            coord = np.zeros(1, dtype=np.float32)
            Rcoord = np.zeros(1, dtype=np.float32)
        coordP = coord.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        RcoordP = Rcoord.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        c_lib.rootMain(ctypes.c_char_p(inStr), partitionsP, ctypes.c_double(alku), ctypes.c_double(loppu), ctypes.c_bool(options.source), ctypes.c_uint32(options.linear_multip), cryst_per_blockP, 
                       ctypes.c_uint32(options.blocks_per_ring), det_per_ringP, CP, SCP, RAP, trIndicesP, axIndicesP, DtrIndicesP, DaxIndicesP, ctypes.c_bool(options.obtain_trues), ctypes.c_bool(options.store_scatter), ctypes.c_bool(options.store_randoms), 
                       scatterComP, ctypes.c_bool(options.randoms_correction), coordP, RcoordP, store_coordinates, cryst_per_block_axialP, ctypes.c_uint32(options.transaxial_multip), 
                       ringsP, sinoSizeP, ctypes.c_uint32(options.Ndist), NangP, ctypes.c_uint32(options.ring_difference), ctypes.c_uint32(options.span),
                       segP, ctypes.c_int64(Nt), ctypes.c_uint64(TOFSize), ctypes.c_int32(options.ndist_side), SinoP, SinoTP, SinoCP, SinoRP, SinoDP, det_w_pseudoP, ctypes.c_uint32(nPseudos), 
                       ctypes.c_double(options.TOF_width), ctypes.c_double(FWHM), ctypes.c_bool(options.verbose), ctypes.c_int32(options.nLayers), ctypes.c_float(dx), ctypes.c_float(dy), ctypes.c_float(dz),
                       ctypes.c_float(bx), ctypes.c_float(by), ctypes.c_float(bz), ctypes.c_int64(Nx), ctypes.c_int64(Ny), ctypes.c_int64(Nz), ctypes.c_bool(options.dualLayerSubmodule), 
                       ctypes.c_bool(options.useIndexBasedReconstruction), tPointP)
        if store_coordinates:
            if Nt <= 1:
                np.append(Fcoord, coord);
                if options.randoms_correction:
                    np.append(FDcoord, Rcoord);
            else:
                if options.randoms_correction:
                    print('Coordinates for delayed coincidences are not supported for list-mode data in dynamic mode!')
                uind = np.unique(tPoints)
                for uu in range(0, len(uind)):
                    np.append(Fcoord[uu], coord);
        elif options.useIndexBasedReconstruction:
            np.append(Fcoord, trIndices);
            np.append(FDcoord, axIndices);
            if options.randoms_correction:
                np.append(DtrIndex, DtrIndices);
                np.append(DaxIndex, DaxIndices);
        print('File ' + rootFile + ' loaded')
    return Sino, SinoT, SinoC, SinoR, SinoD, Fcoord, FDcoord, DtrIndex, DaxIndex