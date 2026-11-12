# -*- coding: utf-8 -*-
"""
Quantification factor (Bq/mL) of the Siemens Biograph Vision images reconstructed with OMEGA.
"""

import os
import re
import warnings

import numpy as np

# Empirical constant (see quantificationFactorBiograph), unitless
_C = 2.19e-3
# Divisor of the block singles rate (empirical, see quantificationFactorBiograph)
_RATE_DIVISOR = 210.
# Vision block layout
_N_AX_BLOCKS = 8
_N_TR_BLOCKS = 38
_CRYST_AX = 10  # Axial crystals per block
_CRYST_TR = 20  # Transaxial crystals per block
_N_RINGS = _N_AX_BLOCKS * _CRYST_AX
_N_CRYST = _N_TR_BLOCKS * _CRYST_TR


def _readHeaderText(fpath, sidecar):
    """
    Reads the Interfile header that is embedded in a .ptd file (at the end of the file; the
    beginning of the file is checked as well) or, if sidecar is True, stored in the sidecar
    <file>.hdr. Returns an empty string if there is no header. Same logic as _readNormHeader
    in loadNormalizationBiograph.
    """
    if not os.path.isfile(fpath):
        raise FileNotFoundError('The file ' + str(fpath) + ' does not exist')
    fsz = os.path.getsize(fpath)
    nBuf = min(fsz, 4 * 2**20)
    with open(fpath, 'rb') as f:
        f.seek(fsz - nBuf)
        buf = f.read(nBuf)
        pos = buf.find(b'!INTERFILE')
        if pos < 0 and fsz > nBuf:
            f.seek(0)
            buf = f.read(nBuf)
            pos = buf.find(b'!INTERFILE')
    txt = ''
    if pos >= 0:
        txt = buf[pos:].decode('latin-1')
        # The header ends where the binary (DICOM) data continues
        z = txt.find('\x00')
        if z >= 0:
            txt = txt[:z]
    elif sidecar:
        base = os.path.splitext(fpath)[0]
        for hdr in (fpath + '.hdr', base + '.hdr'):
            if os.path.isfile(hdr):
                with open(hdr, 'r', encoding='latin-1') as f:
                    txt = f.read()
                break
    return txt


def _headerString(txt, key):
    """Value of the Interfile key (regular expression) as a string, empty if the key is not in the header"""
    m = re.search(r'^[!%]?' + key + r'\s*:=[ \t]*([^\r\n]*)', txt, re.MULTILINE)
    return '' if m is None else m.group(1).strip()


def _headerNumber(txt, key, what, fpath):
    """Numeric value of the Interfile key (regular expression), error if the key is missing or not a number"""
    s = _headerString(txt, key)
    try:
        v = float(s)
    except ValueError:
        v = float('nan')
    if len(s) == 0 or not np.isfinite(v):
        raise ValueError('The header of ' + str(fpath) + ' does not contain ' + what)
    return v


def quantificationFactorBiograph(normFile, acquisitionFile, options, halfLife=None, branchingFactor=None,
                                 frameDuration=None):
    """
    Computes the scalar K that converts an image reconstructed with OMEGA from measured Siemens
    Biograph Vision data to Bq/mL:
        image_in_Bq_per_mL = K * x,
    where x is an image that was reconstructed with the normalization of
    omegatomo.fileio.loadNormalizationBiograph (the Siemens PET_CALIBRATION .ptd file) and with
    the same attenuation, randoms and scatter corrections as the Siemens reconstruction. The
    factor is

        K = C * QF * xcal * DF * DT / (BF * T * mRef),

    where
        QF   = scanner quantification factor (Bq*s/ECAT counts) of the normalization file header,
        xcal = cross calibration factor of the normalization file header,
        BF   = branching factor of the isotope (acquisition header),
        T    = frame duration in seconds (acquisition header),
        DF   = lambda * T / (1 - exp(-lambda * T)), the decay of the isotope during the frame,
               lambda = log(2) / half-life. The reference time is the start of the frame (the
               Siemens decay correction to the scan start); decay from the time of the injection
               is NOT applied,
        DT   = 1 / L^2, the dead-time factor. L is the mean live-time fraction of the crystals,
               from the dead-time table (component 5) of the normalization file and the block
               singles rates of the acquisition header (the rate of a block divided by 210 is the
               table index of its 20 x 10 crystals; the divisor 210 is empirical and reproduces
               the Siemens dead-time factors),
        mRef = the mean efficiency of the LORs of the normalization with the LOR length factor
               (Dref / D)^n, i.e. the scale returned by loadNormalizationBiograph with
               returnScale=True, which depends on the sinogram size (Ndist) and on
               options.visionLORLengthExponent (n) and Dref = options.diameter + 2 * options.DOI,
        C    = an empirical constant, see below.
    The result multiplies the final image, it is not applied inside the reconstruction.

    The constant C is empirical (it has no physical interpretation). It was determined by
    comparing OMEGA reconstructions with the Siemens (e-soft) reconstructions in Bq/mL of three
    measured datasets (a Derenzo phantom, a torso phantom and a pig, mean of 2.178e-3, 2.218e-3
    and 2.183e-3).

    Validity of C: it was determined with non-TOF reconstructions (projector type 1, 3 iterations of 10 subsets) of
    sinograms with 520 radial bins and span 19, with the attenuation, normalization, randoms and scatter corrections.
    For the same data the factor was unchanged (within 0.2 %) by the projector type (1, 2 or 3) of non-TOF
    reconstructions, and, through mRef, within 1 % by Ndist (520 or 400), the angular mashing, the span (1 or 19) and
    options.visionLORLengthExponent (0 or 2). C has not been validated with TOF reconstructions. The TOF projectors
    have the same scale as the non-TOF ones, but TOF images need the randoms and scatter corrected per TOF bin, and
    uncorrected randoms or scatter change the image level differently with and without TOF. The factor has only been
    validated for non-TOF reconstructions.

    Siemens headers can be wrong about the isotope. The isotope name, the half-life and the
    branching factor are read from the acquisition header and can be replaced with the keyword
    arguments, e.g. a Ga-68 scan that was labelled F-18 needs halfLife=4062.6 and
    branchingFactor=0.8891.

    Only the Biograph Vision is supported (the dead-time table layout of the mCT is different).

    Parameters
    ----------
    normFile : str
        Path to the Vision normalization file, i.e. the PET_CALIBRATION .ptd file (the header
        with QF and xcal is required).
    acquisitionFile : str
        Path to the PET_EM_SINO .ptd file of the frame. The header of a sinogram file contains
        the isotope, the frame duration and the block singles rates. A PET_LISTMODE .ptd file
        contains the isotope and the duration, but not the singles rates: a warning is then
        shown and DT = 1 is used.
    options : class object
        The options object that is used for the reconstruction (the sinogram size and the
        detector geometry as in loadNormalizationBiograph).
    halfLife : float, optional
        Half-life of the isotope (s), replaces the header value.
    branchingFactor : float, optional
        Branching factor, replaces the header value.
    frameDuration : float, optional
        Frame duration (s), replaces the header value.

    Returns
    -------
    factor : float
        K. Multiply the reconstructed image with it to get Bq/mL.
    components : dict
        The keys C, QF, xcal, DF, DT, BF, T, halfLife, mRef, livetime (L, NaN if the singles
        are not available), isotope (the isotope name of the header), Dref and singlesFound.

    Raises
    ------
    ValueError
        If the files are not Vision files or the header lacks a required field.
    FileNotFoundError
        If a file does not exist.

    Examples
    --------
    >>> K, comp = quantificationFactorBiograph(options.normalization_datafile, acquisitionFile, options)
    >>> image = image * K
    >>> # Ga-68 scan that is labelled as F-18
    >>> K, comp = quantificationFactorBiograph(normFile, acquisitionFile, options, halfLife=4062.6,
    ...                                        branchingFactor=0.8891)
    """
    from .loadNormalizationBiograph import loadNormalizationBiograph
    for name, val in (('halfLife', halfLife), ('branchingFactor', branchingFactor), ('frameDuration', frameDuration)):
        if val is not None and (np.asarray(val).size != 1 or not np.isfinite(float(val)) or float(val) <= 0):
            raise ValueError(name + ' has to be a positive finite scalar')

    # Normalization file header
    txtN = _readHeaderText(normFile, True)
    if len(txtN) == 0:
        raise ValueError('No Interfile header was found in ' + str(normFile) + '. The normalization file has to be '
                         'the PET_CALIBRATION .ptd file of the Vision.')
    offsets = re.findall(r'data offset in bytes\[(\d+)\]\s*:=\s*(\d+)', txtN)
    sizes = re.findall(r'matrix size\[(\d+)\]\s*:=\s*\{([^}]*)\}', txtN)
    off = {int(k): int(v) for k, v in offsets}
    dim = {int(k): tuple(int(x) for x in v.replace(',', ' ').split()) for k, v in sizes}
    if len(off) < 5 or len(dim) < 5 or dim.get(3) != (798, _N_RINGS) or dim.get(5) != (6000, 42, _N_RINGS):
        raise ValueError(str(normFile) + ' is not a Biograph Vision normalization file with the dead-time table (five '
                         'components, crystal efficiencies 798 x 80, dead-time table 6000 x 42 x 80). Only the Vision '
                         'is supported.')
    QF = _headerNumber(txtN, r'scanner quantification factor \(Bq\*s/ECAT counts\)', 'the scanner quantification factor',
                       normFile)
    xcal = _headerNumber(txtN, r'cross calibration factor', 'the cross calibration factor', normFile)

    # Acquisition header
    txtA = _readHeaderText(acquisitionFile, False)
    if len(txtA) == 0:
        raise ValueError('No Interfile header was found in ' + str(acquisitionFile) + '. The acquisition file has to be '
                         'the PET_EM_SINO (or PET_LISTMODE) .ptd file.')
    isotope = _headerString(txtA, r'isotope name')
    if halfLife is None:
        halfLife = _headerNumber(txtA, r'isotope gamma halflife \(sec\)', "the half-life (give 'halfLife')", acquisitionFile)
    if branchingFactor is None:
        branchingFactor = _headerNumber(txtA, r'isotope branching factor', "the branching factor (give 'branchingFactor')",
                                        acquisitionFile)
    if frameDuration is None:
        frameDuration = _headerNumber(txtA, r'image duration \(sec\)', "the frame duration (give 'frameDuration')",
                                      acquisitionFile)
    halfLife = float(halfLife)
    BF = float(branchingFactor)
    T = float(frameDuration)
    if not (halfLife > 0 and BF > 0 and T > 0):
        raise ValueError('The half-life, branching factor and frame duration have to be positive')

    # Decay during the frame, referred to the start of the frame
    lam = np.log(2.) / halfLife
    DF = float(lam * T / (-np.expm1(-lam * T)))

    # Dead-time factor
    singles = None
    rateBlocks = re.findall(r'singles rate\[(\d+)\]\s*:=\s*\{([^}]*)\}', txtA)
    if len(rateBlocks) > 0:
        nAx = _headerNumber(txtA, r'number of axial blocks', 'the number of axial blocks', acquisitionFile)
        nTr = _headerNumber(txtA, r'number of transaxial blocks', 'the number of transaxial blocks', acquisitionFile)
        if int(nAx) != _N_AX_BLOCKS or int(nTr) != _N_TR_BLOCKS:
            raise ValueError('The acquisition header has %d x %d detector blocks, the Vision has %d x %d. Only the Vision '
                             'is supported.' % (int(nAx), int(nTr), _N_AX_BLOCKS, _N_TR_BLOCKS))
        nBlocks = _N_AX_BLOCKS * _N_TR_BLOCKS
        sg = []
        for _, vals in rateBlocks:
            v = np.array([float(x) for x in vals.replace(',', ' ').split()])
            # The first value is not a block (the blocks follow it)
            if v.size == nBlocks + 1:
                v = v[1:]
            if v.size != nBlocks:
                raise ValueError('The singles rate of the acquisition header has %d values, %d expected' % (v.size, nBlocks))
            sg.append(v)
        # Several sampling points: the mean rate of the frame
        singles = np.mean(np.stack(sg), axis=0)
    if singles is None:
        warnings.warn('The block singles rates were not found in the header of ' + str(acquisitionFile) + ' (a list-mode '
                      'file does not contain them), the dead-time factor is set to one. Use the PET_EM_SINO file for the '
                      'dead-time factor.')
        DT = 1.
        livetime = float('nan')
    else:
        # Block rates [axial block, transaxial block]
        rates = singles.reshape((_N_AX_BLOCKS, _N_TR_BLOCKS))
        nTab = 6000 * 42 * _N_RINGS
        tab = np.fromfile(normFile, dtype='<f4', count=nTab, offset=off[5])
        if tab.size != nTab:
            raise ValueError('The normalization file ended before the end of the dead-time table')
        tab = tab.reshape((_N_RINGS, 42, 6000)).astype(np.float64)
        # Live-time fraction of every crystal: the table index is the block rate / 210 (rounded half up) and the
        # transaxial position of the table skips the gap crystals 0 and 21; odd blocks use the positions 22-41
        ring = np.arange(_N_RINGS)
        t = np.arange(_N_CRYST)
        x = rates[ring // _CRYST_AX][:, t // _CRYST_TR] / _RATE_DIVISOR
        xi = np.clip(np.floor(x + 0.5), 0, 5999).astype(np.int64)
        pos = ((t // _CRYST_TR) % 2) * 21 + 1 + (t % _CRYST_TR)
        livetime = float(tab[ring[:, None], pos[None, :], xi].mean())
        if not (livetime > 0) or not np.isfinite(livetime):
            raise ValueError('The live-time fraction from the dead-time table is not positive')
        DT = 1. / livetime**2

    # Mean efficiency of the LORs with the length factor (Dref / D)^n
    _, mRef = loadNormalizationBiograph(normFile, options, returnScale=True)
    mRef = float(mRef)
    if not (mRef > 0):
        raise ValueError('The mean efficiency of the normalization could not be determined')

    factor = float(_C * QF * xcal * DF * DT / (BF * T * mRef))

    Dref = float(np.asarray(options.diameter).ravel()[0])
    DOI = getattr(options, 'DOI', None)
    if DOI is not None and np.asarray(DOI).size > 0:
        Dref += 2. * float(np.asarray(DOI).ravel()[0])
    components = {'C': _C, 'QF': QF, 'xcal': xcal, 'DF': DF, 'DT': DT, 'BF': BF, 'T': T, 'halfLife': halfLife,
                  'mRef': mRef, 'livetime': livetime, 'isotope': isotope, 'Dref': Dref,
                  'singlesFound': singles is not None}
    return factor, components
