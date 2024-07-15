# -*- coding: utf-8 -*-
"""
Created on Mon May  6 11:14:26 2024

@author: Ville-Veikko Wettenhovi
"""

import os
import numpy as np

def loadSPECTInterfile(options):

    try:
        fid = open(options.fpath)
    except FileNotFoundError:
        print('Specified file was not found! Please select SPECT header file')
        from tkinter import Tk, filedialog
        Tk().withdraw()  # Hide Tkinter root window
        fpath = filedialog.askopenfilename(title='Select SPECT header file')
        if not fpath:
            raise ValueError('No file was selected')
        fid = open(os.path.join(fpath))
    
    lines = fid.readlines()
    # hdr = [line.strip() for line in lines]
    fid.close()
    
    uu = 0
    for line in lines:
        if '!number of projections' in line:
            options.nProjections = int(line.split('=')[1].strip())
        elif 'start angle' in line:
            options.startAngle = float(line.split('=')[1].strip())
        elif 'number of detector heads' in line:
            options.nHeads = int(line.split('=')[1].strip())
        elif 'extent of rotation' in line:
            extRot = float(line.split('=')[1].strip())
        elif 'Radius Per View' in line:
            if uu == 0:
                options.radiusPerProj = np.zeros(options.nProjections, dtype=np.float32)
            options.radiusPerProj[uu] = float(line.split('=')[1].strip())
            uu += 1
    options.angleIncrement = extRot / options.nProjections
    # for kk in range(options.nProjections):
    #     pattern = rf'Radius Per View \[{kk+1}\] := (\d+(?:\.\d+)?)'
    #     match = re.search(pattern, '\n'.join(hdr))
    #     if match:
    #         options.radiusPerProj[kk] = float(match.group(1))
    #     else:
    #         raise ValueError(f"Radius Per View [{kk+1}] not found in the file")