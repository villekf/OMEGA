# -*- coding: utf-8 -*-

import numpy as np
import os
import tkinter as tk
from tkinter.filedialog import askopenfilename

def loadProjectionImages(nProjections, binning = 1, fpath = ''):
    from imageio import imread
    def select_file(fpath):
        file = None
    
        if len(fpath) == 0:
            root = tk.Tk()
            root.withdraw()
            nimi = askopenfilename(title='Select first projection image',filetypes=([('Image Files','*.bmp *.tif *.tiff')]))
            if not file:
                raise ValueError('No file was selected')
        else:
            if not os.path.exists(fpath):
                print('Specified file was not found! Please select first projection image')
                root = tk.Tk()
                root.withdraw()
                nimi = askopenfilename(title='Select first projection image',filetypes=([('Image Files','*.bmp *.tif *.tiff')]))
                if not file:
                    raise ValueError('No file was selected')
            else:
                nimi = fpath
        
        file = os.path.split(nimi)[1]
        fpath = os.path.split(nimi)[0] + '/'
        return file, fpath
    
    file, fpath = select_file(fpath)
    
    t = file.rfind('_')
    if t == -1:
        t = file.index('0')
        t = t - 1
    else:
        t = t
    suffix = file.rfind('.')
    nN = suffix - (t + 1)
    suffix = file[suffix+1:]
    
    pFile = file[:t+1]
    ch = ''
    for uu in range(nN - 1):
        ch = ch + '0'
    ch = ch + '0'
    lFile = fpath + pFile + ch + '.' + suffix
    zeroStart = True
    
    if not os.path.exists(lFile):
        ch = ch[:-1] + '1'
        lFile = fpath + pFile + ch + '.' + suffix
        zeroStart = False
    
    A = imread(lFile)
    projData = np.zeros([int(A.shape[0]/binning), int(A.shape[1]/binning), nProjections], dtype=A.dtype, order='F')
    summa = False
    for kk in range(nProjections):
        if zeroStart:
            ch = str(kk)
        else:
            ch = str(kk + 1)
        while len(ch) < nN:
            ch = '0' + ch
        lFile = os.path.join(fpath, pFile + ch + '.' + suffix)
        A = imread(lFile)
        if binning > 1:
            s = A.shape
            if not summa:
                B = np.sum(A.reshape(binning, -1, order='F'), axis=0)
            if (np.max(B) >= 2**16 - 1 and A.dtype == 'uint16') or summa:
                B = np.sum(A.reshape(-1, binning, order='F'), axis=1, dtype=np.uint32)
                B = np.sum(B.reshape(s[0] // binning, binning, -1, order='F'), axis=1, dtype=np.uint32)
                summa = True
                # A = np.mean(A.reshape(-1, binning), axis=1)
                # A = np.uint16(np.mean(A.reshape(-1, binning), axis=1))
            else:
                B = np.sum(B.reshape(s[0] // binning, binning, -1, order='F'), axis=1)
                if np.max(B) >= 2**16 - 1 and A.dtype == 'uint16':
                    B = np.sum(A.reshape(-1, binning, order='F'), axis=1, dtype=np.uint32)
                    B = np.sum(B.reshape(s[0] // binning, binning, -1, order='F'), axis=1, dtype=np.uint32)
                    summa = True
            B.reshape((s[0] // binning, s[1] // binning), order='F')
            A = B
        if A.dtype == 'uint32':
            if not projData.dtype == 'uint32':
                projData = np.uint32(projData)
        projData[:, :, kk] = A
    return projData