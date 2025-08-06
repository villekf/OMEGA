# -*- coding: utf-8 -*-


def loadInveonData(options, store_coordinates = False):
    import os
    import ctypes
    import numpy as np
    
    if len(options.fpath) == 0:
        import tkinter as tk
        from tkinter.filedialog import askopenfilename
        root = tk.Tk()
        root.withdraw()
        nimi = askopenfilename(title='Select Inveon list-mode datafile',filetypes=([('lst Files','*.lst')]))
        if len(nimi) == 0:
            raise ValueError("No file selected")
    else:
        if os.path.exists(options.fpath):
            nimi = options.fpath
        else:
            import tkinter as tk
            from tkinter.filedialog import askopenfilename
            root = tk.Tk()
            root.withdraw()
            nimi = askopenfilename(title='Select Inveon list-mode datafile',filetypes=([('lst Files','*.lst')]))
            if len(nimi) == 0:
                raise ValueError("No file selected")
            
    totSinos = options.TotSinos;
    if options.span == 1:
        totSinos = options.rings**2
    sinoSize = options.Ndist * options.Nang * totSinos
    Nentries = os.path.getsize(nimi) // 6
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
    Sino = np.zeros((options.Ndist, options.Nang, totSinos, Nt), dtype=np.uint16, order='F')
    if options.randoms_correction:
        SinoD = np.zeros((options.Ndist, options.Nang, totSinos, Nt), dtype=np.uint16, order='F')
    else:
        SinoD = np.zeros(1, dtype=np.uint16, order='F')
    seg = np.cumsum(options.segment_table)
    
    fPath = os.path.dirname( __file__ )
    if os.path.exists(os.path.join(fPath, '..', 'util', 'usingPyPi.py')):
        libdir = os.path.join(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')), "libs")
    else:
        libdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..'))
    if os.name == 'posix':
        libname = str(os.path.join(libdir,"inveon.so"))
    elif os.name == 'nt':
        libname = str(os.path.join(libdir,"inveon.dll"))
    else:
        libname = str(os.path.join(libdir,"inveon.so"))
    c_lib = ctypes.CDLL(libname)
    
    inStr = nimi.encode('utf-8')
    if store_coordinates or options.useIndexBasedReconstruction:
        storeL = True
    else:
        storeL = False
    if storeL:
        tPoints = np.zeros(Nentries, dtype=np.uint16)
    else:
        tPoints = np.zeros(1, dtype=np.uint16)
    tPointP = tPoints.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    segP = seg.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    SinoP = Sino.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    SinoDP = SinoD.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    partitionsP = options.partitions.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    if storeL:
        LL1 = np.zeros(Nentries, dtype=np.uint16)
        LL2 = np.zeros(Nentries, dtype=np.uint16)
        if options.randoms_correction:
            DD1 = np.zeros(Nentries, dtype=np.uint16)
            DD2 = np.zeros(Nentries, dtype=np.uint16)
        else:
            DD1 = np.zeros(1, dtype=np.uint16)
            DD2 = np.zeros(1, dtype=np.uint16)
    else:
        LL1 = np.zeros(1, dtype=np.uint16)
        LL2 = np.zeros(1, dtype=np.uint16)
        DD1 = np.zeros(1, dtype=np.uint16)
        DD2 = np.zeros(1, dtype=np.uint16)
    LP1 = LL1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    DP1 = DD1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    LP2 = LL2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    DP2 = DD2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
    val = c_lib.inveonMain(partitionsP, ctypes.c_double(alku), ctypes.c_double(loppu), ctypes.c_uint64(Nentries), ctypes.c_uint32(options.detectors), ctypes.c_bool(options.randoms_correction), ctypes.c_uint64(sinoSize), 
                   ctypes.c_bool(False), ctypes.c_uint32(options.Ndist), ctypes.c_uint32(options.Nang), ctypes.c_uint32(options.ring_difference), ctypes.c_uint32(options.span),
                   segP, ctypes.c_uint64(Nt), ctypes.c_int32(options.ndist_side), ctypes.c_bool(storeL), ctypes.c_char_p(inStr), LP1, LP2, tPointP, DP1, DP2, SinoP, SinoDP)
    if val == 0:
        inStr = nimi.encode('latin-1')
        val = c_lib.inveonMain(partitionsP, ctypes.c_double(alku), ctypes.c_double(loppu), ctypes.c_uint64(Nentries), ctypes.c_uint32(options.detectors), ctypes.c_bool(options.randoms_correction), ctypes.c_uint64(sinoSize), 
                       ctypes.c_bool(False), ctypes.c_uint32(options.Ndist), ctypes.c_uint32(options.Nang), ctypes.c_uint32(options.ring_difference), ctypes.c_uint32(options.span),
                       segP, ctypes.c_uint64(Nt), ctypes.c_int32(options.ndist_side), ctypes.c_bool(storeL), ctypes.c_char_p(inStr), LP1, LP2, tPointP, DP1, DP2, SinoP, SinoDP)
    coordinate = None
    Rcoordinate = None
    DtrIndices = None
    DaxIndices = None
    if store_coordinates:
        from omegatomo.projector.detcoord import detectorCoordinates
        if Nt > 1:
            coordinate = [None] * Nt

        LL1 = LL1[LL1 != 0]
        LL2 = LL2[LL2 != 0]
        
        ring_number1 = (LL1 - 1) // options.det_per_ring
        ring_number2 = (LL2 - 1) // options.det_per_ring
        ring_pos1 = (LL1 - 1) % options.det_per_ring
        ring_pos2 = (LL2 - 1) % options.det_per_ring
        
        x, y = detectorCoordinates(options)
        z_length = options.rings * options.cr_pz
        z = np.linspace(0, z_length, options.rings + 1, dtype=np.float32)
        z = z[1:] - z[1] / 2
        z = z - options.axial_fov / 2 + (options.axial_fov - options.cr_pz * options.rings) / 2
        x = x.astype(np.float32)
        y = y.astype(np.float32)

        if Nt > 1:
            tPoints = tPoints[tPoints != 0]
            for uu in range(Nt):
                ind1 = np.where(tPoints == uu + 1)[0][0]
                ind2 = np.where(tPoints == uu + 1)[0][-1]
                tempPos1 = ring_pos1[ind1:ind2 + 1]
                tempPos2 = ring_pos2[ind1:ind2 + 1]
                tempNum1 = ring_number1[ind1:ind2 + 1]
                tempNum2 = ring_number2[ind1:ind2 + 1]
                coordinate[uu] = np.array([x[tempPos1 - 1], y[tempPos1 - 1], z[tempNum1 - 1],
                                           x[tempPos2 - 1], y[tempPos2 - 1], z[tempNum2 - 1]], order='F')
        else:
            coordinate = np.array([x[ring_pos1], y[ring_pos1], z[ring_number1],
                                   x[ring_pos2], y[ring_pos2], z[ring_number2]], order='F')
        del LL1, LL2

        if options.randoms_correction:
            if Nt > 1:
                Rcoordinate = [None] * Nt
            
            DD1 = DD1[DD1 != 0]
            DD2 = DD2[DD2 != 0]
            
            ring_number1 = np.uint16(np.floor_divide(DD1 - 1, options.det_per_ring))
            ring_number2 = np.uint16(np.floor_divide(DD2 - 1, options.det_per_ring))
            ring_pos1 = np.uint16(np.mod(DD1 - 1, options.det_per_ring))
            ring_pos2 = np.uint16(np.mod(DD2 - 1, options.det_per_ring))
            
            if Nt > 1:
                for uu in range(Nt):
                    ind1 = np.argmax(tPoints == uu + 1)
                    ind2 = np.size(tPoints) - np.argmax(tPoints[::-1] == uu + 1) - 1
                    tempPos1 = ring_pos1[ind1:ind2 + 1]
                    tempPos2 = ring_pos2[ind1:ind2 + 1]
                    tempNum1 = ring_number1[ind1:ind2 + 1]
                    tempNum2 = ring_number2[ind1:ind2 + 1]
                    Rcoordinate[uu] = np.array([x[tempPos1], y[tempPos1], z[tempNum1], x[tempPos2], y[tempPos2], z[tempNum2]], order='F')
            else:
                Rcoordinate = np.array([x[ring_pos1], y[ring_pos1], z[ring_number1], x[ring_pos2], y[ring_pos2], z[ring_number2]], order='F')
            del DD1, DD2
    elif options.useIndexBasedReconstruction:
        LL1 = LL1[LL1 != 0]
        LL2 = LL2[LL2 != 0]
        ring_number1 = np.uint16(np.floor_divide(LL1 - 1, options.det_per_ring))
        ring_number2 = np.uint16(np.floor_divide(LL2 - 1, options.det_per_ring))
        ring_pos1 = np.uint16(np.mod(LL1 - 1, options.det_per_ring))
        ring_pos2 = np.uint16(np.mod(LL2 - 1, options.det_per_ring))
        ring_pos1 = np.reshape(ring_pos1, (-1, 1))
        ring_pos2 = np.reshape(ring_pos2, (-1, 1))
        ring_number1 = np.reshape(ring_number1, (-1, 1))
        ring_number2 = np.reshape(ring_number2, (-1, 1))
        coordinate = np.asfortranarray(np.concatenate((ring_pos1.T, ring_pos2.T), axis = 0))
        Rcoordinate = np.asfortranarray(np.concatenate((ring_number1.T, ring_number2.T), axis = 0))
        if options.randoms_correction:
            DD1 = DD1[DD1 != 0]
            DD2 = DD2[DD2 != 0]
            ring_number1 = np.uint16(np.floor_divide(DD1 - 1, options.det_per_ring))
            ring_number2 = np.uint16(np.floor_divide(DD2 - 1, options.det_per_ring))
            ring_pos1 = np.uint16(np.mod(DD1 - 1, options.det_per_ring))
            ring_pos2 = np.uint16(np.mod(DD2 - 1, options.det_per_ring))
            ring_pos1 = np.reshape(ring_pos1, (-1, 1))
            ring_pos2 = np.reshape(ring_pos2, (-1, 1))
            ring_number1 = np.reshape(ring_number1, (-1, 1))
            ring_number2 = np.reshape(ring_number2, (-1, 1))
            DtrIndices = np.asfortranarray(np.concatenate((ring_pos1.T, ring_pos2.T), axis = 0))
            DaxIndices = np.asfortranarray(np.concatenate((ring_number1.T, ring_number2.T), axis = 0))
        
                
    return Sino, SinoD, coordinate, Rcoordinate, DtrIndices, DaxIndices