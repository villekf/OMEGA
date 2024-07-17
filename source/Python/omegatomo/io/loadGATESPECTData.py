# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:34:05 2024

@author: Ville-Veikko Wettenhovi
"""
import numpy as np
import os
import tkinter as tk
from tkinter.filedialog import askopenfilename
from omegatomo.io.loadInterfile import loadInterfile

def loadGATESPECTData(options):
    if not os.path.isfile(options.fpath):
        print("Specified file was not found! Please select first projection image")
        root = tk.Tk()
        root.withdraw()
        options.fpath = askopenfilename(filetypes=[("SPECT projection data files", "*.sin")],
                                       title="Select the first SPECT projection data file")
        if not options.fpath:
            raise FileNotFoundError("No file was selected")
    fpath, file = os.path.split(options.fpath)
    fpath = os.path.abspath(fpath)

    nimi = os.path.join(fpath, file)
    nimi = nimi[:-4] + '.hdr'
    proj = [None] * 1000
    for kk in range(1000):
        try:
            proj[kk] = loadInterfile(nimi)
        except FileNotFoundError:
            break
        if kk < 10:
            nimi = nimi[:-5] + str(kk + 1) + '.hdr'
        else:
            nimi = nimi[:-6] + str(kk + 1) + '.hdr'
    loppu = kk

    proj = [p for p in proj if p is not None]
    koko = [p.shape for p in proj]
    koko3 = max([p[2] for p in koko])
    A = np.zeros((proj[0].shape[0], proj[0].shape[1], koko3), dtype=proj[0].dtype, order='F')
    for kk in range(loppu):
        cKoko = proj[kk].shape[2]
        A[:,:,:cKoko] += np.squeeze(proj[kk])
    options.nRowsD = A.shape[0]
    options.nColsD = A.shape[1]
    options.SinM = A
    return options