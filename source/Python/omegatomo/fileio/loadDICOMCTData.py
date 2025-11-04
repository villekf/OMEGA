import os
import numpy as np
import pydicom
import struct

def loadDICOMCTPD(path):
    """
    Automatically loads DICOM CT projection data from a directory.
    Loads both projections and the necessary variables/coordinates.
    Input: path to folder containing DICOM images
    Output: proj (numpy array), vars (dict)
    """
    files = [f for f in os.listdir(path) if f.lower().endswith('.dcm')]
    nFiles = len(files)
    xs = np.zeros(nFiles, dtype=np.float32)
    xd = np.zeros(nFiles, dtype=np.float32)
    ys = np.zeros(nFiles, dtype=np.float32)
    yd = np.zeros(nFiles, dtype=np.float32)
    zs = np.zeros(nFiles, dtype=np.float32)
    zd = np.zeros(nFiles, dtype=np.float32)
    angles = np.zeros(nFiles, dtype=np.float32)

    for kk, fname in enumerate(files):
        fpath = os.path.join(path, fname)
        info = pydicom.dcmread(fpath)
        nColsD = info.Columns
        nRowsD = info.Rows

        r = struct.unpack('f', info[0x7031, 0x1003].value)[0]
        deltar = struct.unpack('f', info[0x7033, 0x100d].value)[0]
        rho = r + deltar

        angles[kk] = struct.unpack('f', info[0x7031, 0x1001].value)[0] - np.pi/2

        deltaphi = struct.unpack('f', info[0x7033, 0x100b].value)[0]
        phi = angles[kk] + deltaphi

        deltaz = struct.unpack('f', info[0x7033, 0x100c].value)[0]
        zs[kk] = struct.unpack('f', info[0x7031, 0x1002].value)[0] + deltaz

        ys[kk] = rho * np.sin(phi)
        xs[kk] = -rho * np.cos(phi)

        data = info.pixel_array
        if kk == 0:
            proj = np.zeros((nRowsD, nColsD, nFiles), dtype=np.float32, order='F')
        
        single_data = np.asfortranarray(data.astype(np.float32) * info.RescaleSlope + info.RescaleIntercept)
        proj[:, :, kk] = np.fliplr(single_data)
        
        dPitchX = struct.unpack('f', info[0x7029, 0x1002].value)[0]
        dPitchY = struct.unpack('f', info[0x7029, 0x1006].value)[0]
        L = float(nRowsD) * dPitchX
        theta = L / r
        dtheta = theta / float(nRowsD)
        detElements = struct.unpack('2f', info[0x7031, 0x1033].value)
        rowStart = detElements[0] - float(nRowsD)/2 - 0.5
        colStart = detElements[1] - float(nColsD)/2 - 0.5
        d = struct.unpack('f', info[0x7031, 0x1031].value)[0]
        yd[kk] = r * np.sin(angles[kk]) - d * np.sin(angles[kk] + dtheta * colStart) - np.sin(rowStart * dtheta) * rowStart * dPitchX
        xd[kk] = -r * np.cos(angles[kk]) + d * np.cos(angles[kk] + dtheta * rowStart) - np.cos(rowStart * dtheta) * rowStart * dPitchX
        zd[kk] = struct.unpack('f', info[0x7031, 0x1002].value)[0] + colStart * dPitchY

    nProjections = nFiles
    vars = {
        'xs': xs,
        'ys': ys,
        'zs': zs,
        'xd': xd,
        'yd': yd,
        'zd': zd,
        'angles': angles,
        'nProjections': nProjections,
        'dPitchX': dPitchX,
        'dPitchY': dPitchY,
        'r': d,
        'sourceToCRot': r
    }
    return proj, vars