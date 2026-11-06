# -*- coding: utf-8 -*-
"""
Python port of saveMetaImage.m. Saves an image in the MetaImage format.
"""

import os
import sys

from ._saveUtil import _resolveType, _castImage, _voxelSize, _asAtLeast3D, _columnMajorBytes

_ELEMENT_TYPES = {'float64': 'MET_DOUBLE', 'float32': 'MET_FLOAT', 'int8': 'MET_CHAR', 'uint8': 'MET_UCHAR', 'int16': 'MET_SHORT',
                  'uint16': 'MET_USHORT', 'int32': 'MET_INT', 'uint32': 'MET_UINT', 'int64': 'MET_LONG_LONG',
                  'uint64': 'MET_ULONG_LONG'}


def saveMetaImage(filename, img, dtype=None, options=None, pixelSize=None, windowsFormat=False, endian='n'):
    """
    Saves the input image in MetaImage format (Python port of saveMetaImage.m). The image data are saved in column-major order
    (the first dimension, x, varies the fastest) to [filename].raw and the header to [filename].mhd. By default the header uses
    the UNIX text format (LF).

    Parameters
    ----------
    filename : str
        Name (and optionally path) of the image and header files, without the file type.
    img : ndarray
        The 3D or 4D image (x, y, z[, t]). Vector images are not supported.
    dtype : str or numpy dtype, optional
        Data type of the saved image, e.g. 'float32', 'uint16' (MATLAB names such as 'single' are also accepted). By default
        the type of the input image. The image is cast to this type before saving.
    options : options object or dict, optional
        The OMEGA options object (or a dict with the same field names, or the image properties of MATLAB). The voxel size
        is taken from it (FOVa_x/Nx, FOVa_y/Ny, axial_fov/Nz, i.e. the field of view divided by the size of the image).
    pixelSize : sequence of float, optional
        Voxel size (mm) in x, y and z if no options are available. Overrides the voxel size from the options. If neither is
        given, the voxel size is 1 mm.
    windowsFormat : bool, optional
        Use the Windows text format (CR LF) in the header. Default False (LF).
    endian : str, optional
        Byte order of the image data: 'n' native (default), 'l' little endian or 'b' big endian.
    """
    img = _asAtLeast3D(img)
    dt = _resolveType(dtype, img)
    img = _castImage(img, dt)
    if endian == 'n':
        endian = 'b' if sys.byteorder == 'big' else 'l'
    if endian not in ('l', 'b'):
        raise ValueError("endian has to be 'n', 'l' or 'b'")
    voxel = list(_voxelSize(img.shape, options, pixelSize))
    # Unlike the MATLAB version, the spacing of the time dimension is also saved (as 1)
    while len(voxel) < img.ndim:
        voxel.append(1.)
    folder, base = os.path.split(filename)
    with open(filename + '.raw', 'wb') as f:
        f.write(_columnMajorBytes(img, dt, '>' if endian == 'b' else '<'))
    msb = endian == 'b'
    L = []
    L.append('ObjectType = Image')
    L.append('NDims = ' + str(img.ndim))
    L.append('DimSize = ' + ' '.join(str(n) for n in img.shape))
    L.append('BinaryData = true')
    L.append('BinaryDataByteOrderMSB = ' + ('true' if msb else 'false'))
    L.append('CompressedData = false')
    nd = img.ndim
    L.append('TransformMatrix = ' + ' '.join('1' if i == j else '0' for i in range(nd) for j in range(nd)))
    L.append('CenterOfRotation = ' + ' '.join(['0'] * nd))
    L.append('ElementType = ' + _ELEMENT_TYPES[dt.name])
    L.append('HeaderSize = 0')
    L.append('ElementSpacing = ' + ' '.join('%.15g' % v for v in voxel))
    L.append('ElementNumberOfChannels = 1')
    L.append('ElementByteOrderMSB = ' + ('True' if msb else 'False'))
    L.append('ElementDataFile = ' + base + '.raw')
    eol = '\r\n' if windowsFormat else '\n'
    with open(filename + '.mhd', 'wb') as f:
        f.write((eol.join(L) + eol).encode('ascii'))
