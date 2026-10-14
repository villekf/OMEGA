

def loadDICOMCTPD(path):
    """
    Automatically loads DICOM CT projection data from a directory.
    Loads both projections and the necessary variables/coordinates.
    Input: path to folder containing DICOM images
    Output: proj (numpy array), vars (dict)
    """
    import os
    import numpy as np
    import pydicom
    import struct
    files = sorted(f for f in os.listdir(path) if f.lower().endswith('.dcm'))
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

        a = struct.unpack('f', info[0x7031, 0x1001].value)[0] - np.pi/2

        deltaphi = struct.unpack('f', info[0x7033, 0x100b].value)[0]
        phi = a + deltaphi

        deltaz = struct.unpack('f', info[0x7033, 0x100c].value)[0]
        zf = struct.unpack('f', info[0x7031, 0x1002].value)[0]
        zs[kk] = zf + deltaz

        ys[kk] = rho * np.sin(phi)
        xs[kk] = -rho * np.cos(phi)

        data = info.pixel_array
        if kk == 0:
            proj = np.zeros((nRowsD, nColsD, nFiles), dtype=np.float32, order='F')
        
        single_data = np.asfortranarray(data.astype(np.float32) * info.RescaleSlope + info.RescaleIntercept)
        proj[:, :, kk] = np.fliplr(single_data)
        
        dPitchX = struct.unpack('f', info[0x7029, 0x1002].value)[0]
        dPitchY = struct.unpack('f', info[0x7029, 0x1006].value)[0]
        d = struct.unpack('f', info[0x7031, 0x1031].value)[0]
        # The detector is an arc of radius d centered at the focal center
        # (-r cos(a), r sin(a)), the channels are spaced by dPitchX / d radians
        dtheta = dPitchX / d
        detElements = struct.unpack('2f', info[0x7031, 0x1033].value)
        # Central channel offset from the geometric center of the channels
        deltaCh = detElements[0] - float(nRowsD)/2 - 0.5
        # The central element is not the middle of the detector, so the
        # detector center is rotated along the arc
        aD = a - deltaCh * dtheta
        # Axial offset of the central row (sign not verified, zero for the AAPM data)
        colStart = detElements[1] - float(nColsD)/2 - 0.5
        xd[kk] = -r * np.cos(a) + d * np.cos(aD)
        yd[kk] = r * np.sin(a) - d * np.sin(aD)
        zd[kk] = zf + colStart * dPitchY
        # The projector derives the arc center from the detector center and the
        # angle, hence the rotated angle has to be stored
        angles[kk] = aD

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