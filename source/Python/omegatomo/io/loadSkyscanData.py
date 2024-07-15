# -*- coding: utf-8 -*-

import numpy as np
import os
import tkinter as tk
from tkinter.filedialog import askopenfilename
from omegatomo.io.loadProjectionImages import loadProjectionImages

def loadSkyscanData(options):
    def select_file():
    
        if len(options.fpath) == 0:
            root = tk.Tk()
            root.withdraw()
            nimi = askopenfilename(title='Select Skyscan log file',filetypes=([('log Files','*.log')]))
            if not nimi:
                raise ValueError('No file was selected')
        else:
            if not os.path.exists(options.fpath):
                print('Specified file was not found! Please select Nikon xtekct file')
                root = tk.Tk()
                root.withdraw()
                nimi = askopenfilename(title='Select Skyscan log file',filetypes=([('log Files','*.log')]))
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
        if hdr[i].split("=")[0] == 'Number of Columns':
            options.nRowsD = int(hdr[i].split("=")[1]) // options.binning
        if hdr[i].split("=")[0] == 'Number of Rows':
            options.nColsD = int(hdr[i].split("=")[1]) // options.binning
        if hdr[i].split("=")[0] == 'Camera Pixel Size (um)':
            options.dPitchX = float(hdr[i].split("=")[1]) * options.binning / 1000
        if hdr[i].split("=")[0] == 'Camera Pixel Size (um)':
            options.dPitchY = float(hdr[i].split("=")[1]) * options.binning / 1000
        if hdr[i].split("=")[0] == 'Camera to Source (mm)':
            options.sourceToDetector = float(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'Object to Source (mm)':
            options.sourceToCRot = float(hdr[i].split("=")[1])
        if hdr[i].split("=")[0] == 'Number of Files':
            options.nProjections = int(hdr[i].split("=")[1]) - 1

    if not hasattr(options, 'only_reconstructions') or not options.only_reconstructions:
        file_path = os.path.join(fpath, file + '0000.tif')
        options.SinM = loadProjectionImages(options.nProjections, options.binning, file_path)
        options.SinM = np.transpose(options.SinM, (1, 0, 2))