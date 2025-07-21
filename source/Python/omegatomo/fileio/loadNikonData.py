# -*- coding: utf-8 -*-


def loadNikonData(options):
    import numpy as np
    import os
    import tkinter as tk
    from tkinter.filedialog import askopenfilename
    from omegatomo.fileio.loadProjectionImages import loadProjectionImages
    def select_file():
    
        if len(options.fpath) == 0:
            root = tk.Tk()
            root.withdraw()
            nimi = askopenfilename(title='Select Nikon xtekct file',filetypes=([('xtekct Files','*.xtekct')]))
            if not nimi:
                raise ValueError('No file was selected')
        else:
            if not os.path.exists(options.fpath):
                print('Specified file was not found! Please select Nikon xtekct file')
                root = tk.Tk()
                root.withdraw()
                nimi = askopenfilename(title='Select Nikon xtekct file',filetypes=([('xtekct Files','*.xtekct')]))
                if not nimi:
                    raise ValueError('No file was selected')
            else:
                nimi = options.fpath
        
        file = os.path.split(nimi)[1]
        file, ext = os.path.splitext(file)
        fpath = os.path.split(nimi)[0] + '/'
        return file, fpath, nimi
    
    file, fpath, nimi = select_file()
    
    
    with open(nimi, 'r') as fid:
        hdr = [line.strip() for line in fid.readlines()]

    if not hasattr(options, 'binning'):
        options.binning = 1
    if not hasattr(options, 'only_reconstructions'):
        options.only_reconstructions = False


    for i in range(len(hdr)):
        if not hasattr(options, 'Nx'):
            if hdr[i].split("=")[0] == 'VoxelsX':
                options.Nx = int(hdr[i].split("=")[1])
        if not hasattr(options, 'Ny'):
            if hdr[i].split("=")[0] == 'VoxelsY':
                options.Ny = int(hdr[i].split("=")[1])
        if not hasattr(options, 'Nz'):
            if hdr[i].split("=")[0] == 'VoxelsZ':
                options.Nz = int(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'DetectorPixelsY':
            options.nRowsD = int(hdr[i].split("=")[1]) // options.binning
        if hdr[i].split("=")[0] == 'DetectorPixelsX':
            options.nColsD = int(hdr[i].split("=")[1]) // options.binning
        if hdr[i].split("=")[0] == 'DetectorPixelSizeY':
            options.dPitchX = float(hdr[i].split("=")[1]) * options.binning
        if hdr[i].split("=")[0] == 'DetectorPixelSizeX':
            options.dPitchY = float(hdr[i].split("=")[1]) * options.binning
        if hdr[i].split("=")[0] == 'SrcToDetector':
            options.sourceToDetector = float(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'SrcToObject':
            options.sourceToCRot = float(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'ObjectOffsetX':
            options.oOffsetX = float(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'ObjectOffsetY':
            options.oOffsetY = float(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'DetectorOffsetY':
            options.detOffsetRow = float(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'DetectorOffsetX':
            options.detOffsetCol = float(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'MaskRadius':
            options.MaskRadius = float(hdr[i].split("=")[1])

    with open(os.path.join(fpath, file + '.ang'), 'r') as fid:
        data = [line.strip() for line in fid.readlines()]
    options.nProjections = (len(data) - 1)
    options.angles = np.zeros(options.nProjections, dtype=np.float32)
    data = data[1:]
    for kk in range(options.nProjections):
        parts = data[kk].split(":")
        options.angles[kk] = -float(parts[1].strip())

    if not hasattr(options, 'only_reconstructions') or not options.only_reconstructions:
        file_path = os.path.join(fpath, file + '_0001.tif')
        options.SinM = loadProjectionImages(options.nProjections, options.binning, file_path)
        options.SinM = np.transpose(options.SinM, (1, 0, 2))