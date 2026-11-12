# -*- coding: utf-8 -*-


def _getTag(txt, tag):
    """Text content of the first <tag>...</tag> in txt (empty string if not found)"""
    import re
    m = re.search('<' + tag + '>([^<]*)</' + tag + '>', txt)
    return m.group(1) if m else ''


def loadElektaXVIData(options):
    """
    Loads Elekta XVI CBCT projection data (.his + _Frames.xml).

    Loads the raw projection images (Heimann HIS files, one per projection) and the geometry information
    (_Frames.xml) of an Elekta XVI CBCT scan, e.g. an XVI "projection export" folder. The projections are
    converted to intensities (the XVI raw values are inverted, i.e. air has the highest value after the
    conversion) and the OMEGA standard CT geometry fields are set (angles, sourceToCRot, sourceToDetector,
    detector pitch, per projection detector offsets), so no custom source/detector coordinates are needed.
    Only frames with Exposed = True are used.

    Inputs (options attributes):
        fpath = The folder that contains the .his files and the *_Frames.xml file. If empty, or not found,
            a folder selection dialog is shown.
        binning = (Optional, default 1) Binning factor. The projections are binned in both dimensions by
            averaging (mean) binning.
        sourceToCRot = (Optional) Source to center of rotation distance (mm). Default is the nominal XVI
            value, 1000 mm. A value that has already been set (> 1) is kept.
        sourceToDetector = (Optional) Source to detector distance (mm). Default is the nominal XVI value,
            1536 mm. A value that has already been set (> 1) is kept.
        dPitchX, dPitchY = (Optional) Detector pixel pitch (mm) AFTER binning, i.e. raw pitch * binning.
            Default is the nominal XVI raw pitch 0.8 mm * binning. Values that have already been set (> 0)
            are kept.
        elektaUCentreScale, elektaVCentreScale = (Optional, default 1) Scale factors that convert the
            UCentre/VCentre values (mm) of the _Frames.xml file to the detector offsets. The nominal
            convention (scale = 1) is used by default. The sign convention is fixed:
            detOffsetRow = elektaUCentreScale * UCentre and detOffsetCol = -elektaVCentreScale * VCentre
            (mm, in the OMEGA detector frame). Other values (also fitted effective values) can be used if
            the nominal ones do not fit the data well.
        elektaFlatPercentile = (Optional, default 99.9) The percentile of the binned intensity values (all
            projections) that is used as the global air (flat) value, options.flat. If options.flat has
            already been set (> 0) it is not modified.

    Outputs:
        options.SinM = The projection intensities, nRowsD x nColsD x nProjections, float32. The first
            dimension is the transaxial detector direction, the second the axial direction.
        options.nRowsD, nColsD, nProjections, angles (radians), sourceToCRot, sourceToDetector, dPitchX,
        dPitchY, detOffsetRow and detOffsetCol (nProjections, mm) and flat (scalar).

    Geometry notes: The XVI gantry angle a corresponds to the OMEGA projection angle a + pi/2 (the source is
    at SAD*(sin a, -cos a, 0) in the XVI frame, OMEGA places it at -SAD*(cos t, sin t, 0)). OMEGA's detector
    row direction is opposite to the direction of increasing column index in the .his files, hence the
    transaxial dimension of the images is reversed. The last detector row of the XVI panel is dead and is
    replaced by the previous row.

    Example:
        options.fpath = '/data/XVI_scan'
        options.binning = 2
        loadElektaXVIData(options)
    """
    import numpy as np
    import os
    import re

    def select_folder():
        if len(options.fpath) == 0:
            import tkinter as tk
            from tkinter.filedialog import askdirectory
            root = tk.Tk()
            root.withdraw()
            fpath = askdirectory(title='Select the Elekta XVI projection folder (.his files and _Frames.xml)')
            if not fpath:
                raise ValueError('No folder was selected')
        else:
            if not os.path.isdir(options.fpath):
                import tkinter as tk
                from tkinter.filedialog import askdirectory
                print('Specified folder was not found! Please select the Elekta XVI projection folder')
                root = tk.Tk()
                root.withdraw()
                fpath = askdirectory(title='Select the Elekta XVI projection folder (.his files and _Frames.xml)')
                if not fpath:
                    raise ValueError('No folder was selected')
            else:
                fpath = options.fpath
        return fpath

    fpath = select_folder()
    binning = int(getattr(options, 'binning', 1) or 1)
    only_reconstructions = bool(getattr(options, 'only_reconstructions', False))
    uScale = getattr(options, 'elektaUCentreScale', None)
    uScale = 1. if uScale is None else float(uScale)
    vScale = getattr(options, 'elektaVCentreScale', None)
    vScale = 1. if vScale is None else float(vScale)
    flatPercentile = getattr(options, 'elektaFlatPercentile', None)
    flatPercentile = 99.9 if flatPercentile is None else float(flatPercentile)

    # Frames file
    xmlFile = sorted([f for f in os.listdir(fpath) if f.endswith('_Frames.xml')])
    if len(xmlFile) == 0:
        raise ValueError('No *_Frames.xml file found in ' + fpath)
    with open(os.path.join(fpath, xmlFile[0]), 'r', encoding='utf-8-sig') as fid:
        txt = fid.read()
    try:
        W = int(_getTag(txt, 'Width'))
        H = int(_getTag(txt, 'Height'))
    except ValueError:
        W = 512
        H = 512
    if W % binning != 0 or H % binning != 0:
        raise ValueError('The binning factor has to divide the image size')
    frames = re.findall('<Frame>(.*?)</Frame>', txt, flags=re.DOTALL)
    seq = np.array([float(_getTag(fr, 'Seq')) for fr in frames])
    gantry = np.array([float(_getTag(fr, 'GantryAngle')) for fr in frames])
    uC = np.array([float(_getTag(fr, 'UCentre')) for fr in frames])
    vC = np.array([float(_getTag(fr, 'VCentre')) for fr in frames])
    exposed = np.array([_getTag(fr, 'Exposed').strip().lower() == 'true' for fr in frames])
    order = np.argsort(seq, kind='stable')
    order = order[exposed[order]]
    seq = seq[order]
    gantry = gantry[order]
    uC = uC[order]
    vC = vC[order]
    nProj = seq.size

    # Geometry
    if not (np.size(options.sourceToCRot) > 0 and np.ravel(options.sourceToCRot)[0] > 1):
        options.sourceToCRot = 1000.
    if not (np.size(options.sourceToDetector) > 0 and np.ravel(options.sourceToDetector)[0] > 1):
        options.sourceToDetector = 1536.
    if not (np.size(options.dPitchX) > 0 and np.ravel(options.dPitchX)[0] > 0):
        options.dPitchX = 0.8 * binning
    if not (np.size(options.dPitchY) > 0 and np.ravel(options.dPitchY)[0] > 0):
        options.dPitchY = 0.8 * binning
    options.nRowsD = W // binning
    options.nColsD = H // binning
    options.nProjections = nProj
    # The XVI gantry angle a is the OMEGA angle a + pi/2 (see above)
    options.angles = (np.deg2rad(gantry) + np.pi / 2).astype(np.float32)
    options.detOffsetRow = (uScale * uC).astype(np.float32)
    options.detOffsetCol = (-vScale * vC).astype(np.float32)

    # Projections
    if not only_reconstructions or np.size(options.SinM) == 0:
        files = sorted([f for f in os.listdir(fpath) if f.endswith('.his')])
        fileSeq = {}
        for f in files:
            m = re.match(r'^(\d+)\.', f)
            if m:
                fileSeq[int(m.group(1))] = f
        SinM = np.zeros((W // binning, H // binning, nProj), dtype=np.float32, order='F')
        for kk in range(nProj):
            if int(seq[kk]) not in fileSeq:
                raise ValueError('No .his file found for frame ' + str(int(seq[kk])))
            fname = os.path.join(fpath, fileSeq[int(seq[kk])])
            # HIS file: 16-bit words at byte 0 (0x7000, file ID), byte 2 (header size, 68) and byte 10
            # (image header size, 32); the image data (uint16, little-endian, rows in the axial direction)
            # follows after both (100 bytes). If this does not match the file size, the image is assumed to
            # be the last part of the file.
            fsize = os.path.getsize(fname)
            with open(fname, 'rb') as fid:
                h = np.frombuffer(fid.read(12), dtype='<u2')
                hdrBytes = int(h[1]) + int(h[5])
                if hdrBytes + W * H * 2 != fsize:
                    hdrBytes = fsize - W * H * 2
                fid.seek(hdrBytes)
                A = np.frombuffer(fid.read(W * H * 2), dtype='<u2').reshape(H, W)
            # A is (axial v, transaxial u) -> (u, v). Raw value inversion: air has the highest intensity
            A = np.float32(65535) - A.T.astype(np.float32)
            # Dead last detector row
            A[:, H - 1] = A[:, H - 2]
            # OMEGA's detector row direction is opposite to the file column order
            A = A[::-1, :]
            if binning > 1:
                A = A.reshape(W // binning, binning, H // binning, binning).sum(axis=(1, 3), dtype=np.float32) / np.float32(binning ** 2)
            SinM[:, :, kk] = A
        options.SinM = SinM
        del SinM

    # Global flat (air) value
    if not (np.size(options.flat) > 0 and np.ravel(options.flat)[0] > 0):
        a = np.sort(options.SinM.ravel())
        k = int(np.ceil(flatPercentile / 100 * a.size))
        k = min(max(k, 1), a.size)
        options.flat = np.float32(a[k - 1])
        del a
    return options
