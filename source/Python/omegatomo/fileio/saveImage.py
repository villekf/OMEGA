# -*- coding: utf-8 -*-
"""
Python port of saveImage.m. Saves reconstructed images in raw, NIfTI, Analyze 7.5, DICOM, Interfile or MetaImage format.
"""

import hashlib
import os
import struct
import time

import numpy as np

from ._saveUtil import _get, _resolveType, _castImage, _voxelSize, _lastIteration, _columnMajorBytes

_OUTPUT_TYPES = ('raw', 'nifti', 'analyze', 'dicom', 'interfile', 'metaimage')


def _importSimpleITK():
    try:
        import SimpleITK as sitk
    except ModuleNotFoundError:
        raise ImportError('SimpleITK package not found! NIfTI and DICOM images cannot be saved. You can install SimpleITK package with "pip install SimpleITK".') from None
    return sitk


def _saveNifti(filename, img, voxel, description):
    """Saves the image as a NIfTI-1 file. As in saveImage.m, the first two dimensions are rotated by 90 degrees (rot90) before
    saving, i.e. the file contains the image with the dimensions (y, x, z[, t]) flipped as rot90 does."""
    sitk = _importSimpleITK()
    dt = _resolveType(None, img)
    img = _castImage(img, dt)
    img = np.rot90(img, 1, axes=(0, 1))
    # Column-major (file) order x fastest: reverse the axes for SimpleITK
    arr = np.ascontiguousarray(img.transpose(tuple(range(img.ndim - 1, -1, -1))))
    image = sitk.GetImageFromArray(arr, isVector=False)
    # After the rotation the first dimension is y and the second one x
    spacing = [voxel[1], voxel[0], voxel[2]] + [1.] * (img.ndim - 3)
    image.SetSpacing(spacing)
    image.SetMetaData('ITK_FileNotes', description[:79])
    sitk.WriteImage(image, filename + '.nii')


_ANALYZE_TYPES = {'uint8': (2, 8), 'int16': (4, 16), 'int32': (8, 32), 'float32': (16, 32), 'float64': (64, 64),
                  # Not part of Analyze 7.5, but used by the NIfTI/Analyze tools for MATLAB
                  'int8': (256, 8), 'uint16': (512, 16), 'uint32': (768, 32)}


def _saveAnalyze(filename, img, voxel, description):
    """Saves the image as an Analyze 7.5 pair of files (.hdr, 348 bytes, and .img), little endian. As in saveImage.m, the image
    is not rotated."""
    dt = _resolveType(None, img)
    if dt.name not in _ANALYZE_TYPES:
        raise ValueError('Data type ' + dt.name + ' is not supported by the Analyze format')
    code, bitpix = _ANALYZE_TYPES[dt.name]
    img = _castImage(img, dt)
    ndim = img.ndim
    dim = [ndim] + list(img.shape) + [1] * (7 - ndim)
    pixdim = [0., voxel[0], voxel[1], voxel[2]] + [0.] * 4
    if ndim == 4:
        pixdim[4] = 1.
    if img.dtype.kind == 'f':
        finite = img[np.isfinite(img)]
    else:
        finite = img
    if finite.size > 0:
        glmax, glmin = float(finite.max()), float(finite.min())
    else:
        glmax = glmin = 0.
    clip = lambda v: int(max(min(v, 2147483647), -2147483648))
    hdr = bytearray(348)
    struct.pack_into('<i', hdr, 0, 348)                      # sizeof_hdr
    hdr[4:14] = b'\0' * 10                                   # data_type
    hdr[14:32] = b'\0' * 18                                  # db_name
    struct.pack_into('<i', hdr, 32, 16384)                   # extents
    hdr[38:39] = b'r'                                        # regular
    struct.pack_into('<8h', hdr, 40, *dim)                   # dim
    struct.pack_into('<h', hdr, 70, code)                    # datatype
    struct.pack_into('<h', hdr, 72, bitpix)                  # bitpix
    struct.pack_into('<8f', hdr, 76, *pixdim)                # pixdim
    struct.pack_into('<f', hdr, 108, 0.)                     # vox_offset
    struct.pack_into('<f', hdr, 112, 1.)                     # roi_scale (scale factor)
    struct.pack_into('<i', hdr, 140, clip(glmax))            # glmax
    struct.pack_into('<i', hdr, 144, clip(glmin))            # glmin
    hdr[148:148 + 80] = description.encode('ascii', errors='replace')[:79].ljust(80, b'\0')   # descrip
    hdr[252] = 0                                             # orient (transverse unflipped)
    with open(filename + '.hdr', 'wb') as f:
        f.write(bytes(hdr))
    with open(filename + '.img', 'wb') as f:
        f.write(_columnMajorBytes(img, dt, '<'))


def _saveDicom(filename, img, voxel, description):
    """Saves the image as a series of 2D DICOM files, [filename]_slice[n].dcm, n = 1, 2, ... The pixel data is 16-bit. Image
    with integer values within the range of a 16-bit type are saved as such, other images are scaled to the full uint16
    range with RescaleSlope and RescaleIntercept, i.e. value = stored * RescaleSlope + RescaleIntercept (DICOM viewers apply
    this automatically), which requires the pydicom package. The image is not rotated: the first dimension (x) is the row
    direction. Note that dicomwrite of MATLAB saves a double image as uint16 by clipping the values to [0, 1]."""
    sitk = _importSimpleITK()
    if img.ndim > 3:
        raise ValueError('4D data is not supported by DICOM conversion')
    a = np.asarray(img, dtype=np.float64)
    ok = np.isfinite(a)
    a = np.where(ok, a, 0.)
    vmin, vmax = float(a.min()), float(a.max())
    integral = bool(np.all(a == np.round(a)))
    rescale = False
    if integral and vmin >= 0 and vmax <= 65535:
        stored = a.astype(np.uint16)
        slope, intercept = 1., 0.
    elif integral and vmin >= -32768 and vmax <= 32767:
        stored = a.astype(np.int16)
        slope, intercept = 1., 0.
    else:
        intercept = vmin
        slope = (vmax - vmin) / 65535. if vmax > vmin else 1.
        # DICOM decimal strings have at most 16 characters
        slope = float('%.9g' % slope)
        intercept = float('%.9g' % intercept)
        stored = np.clip(np.round((a - intercept) / slope), 0, 65535).astype(np.uint16)
        rescale = True
        try:
            import pydicom
        except ModuleNotFoundError:
            raise ImportError('pydicom package not found! DICOM images with non-integer values cannot be saved (the pixel values are saved as 16-bit integers with a rescale slope and intercept). You can install pydicom package with "pip install pydicom".') from None
    # UIDs of the study and the series (the UID root 2.25 is for UUID-derived UIDs)
    def _uid(tag):
        h = hashlib.sha256((filename + tag + str(time.time())).encode('utf-8')).hexdigest()
        return '2.25.' + str(int(h[:32], 16))
    studyUID = _uid('study')
    seriesUID = _uid('series')
    writer = sitk.ImageFileWriter()
    now = time.localtime()
    for k in range(stored.shape[2]):
        sl = np.ascontiguousarray(stored[:, :, k])          # shape (rows, columns) = (x, y)
        im = sitk.GetImageFromArray(sl)
        im.SetSpacing((voxel[1], voxel[0]))                 # (column spacing, row spacing)
        md = {
            '0008|0060': 'OT',
            '0008|0008': 'DERIVED\\SECONDARY',
            '0008|0021': time.strftime('%Y%m%d', now),
            '0008|0031': time.strftime('%H%M%S', now),
            '0008|103e': description[:64],
            '0018|0050': '%.9g' % voxel[2],
            '0028|0030': '%.9g\\%.9g' % (voxel[0], voxel[1]),
            '0020|000d': studyUID,
            '0020|000e': seriesUID,
            '0020|0011': '1',
            '0020|0013': str(k + 1),
            '0020|0032': '0\\0\\%.9g' % (k * voxel[2]),
            # Increasing column index (y) is the +y direction, increasing row index (x) is the +x direction
            '0020|0037': '0\\1\\0\\1\\0\\0',
            '0020|1041': '%.9g' % (k * voxel[2]),
            '0010|0010': '',
            '0010|0020': '',
        }
        for key, value in md.items():
            im.SetMetaData(key, value)
        fname = filename + '_slice' + str(k + 1) + '.dcm'
        writer.SetFileName(fname)
        writer.Execute(im)
        if rescale:
            # The DICOM writer of SimpleITK only supports integer rescale values, set the actual ones with pydicom
            ds = pydicom.dcmread(fname)
            ds.RescaleIntercept = '%.9g' % intercept
            ds.RescaleSlope = '%.9g' % slope
            ds.save_as(fname)


def saveImage(img, outputType='raw', filename=None, options=None, voxelSize=None, reko=None):
    """
    Saves the specified reconstructed image in the specified format (Python port of saveImage.m).

    NIfTI and DICOM require the SimpleITK package, the other formats do not require additional packages. Raw data format is
    used by default. DICOM does not support time-series of images (one file is created for each slice). Interfile, MetaImage
    and raw formats store the voxels in the same column-major order (x varies the fastest) as the files saved by the MATLAB
    version. As in MATLAB, the NIfTI image is rotated by 90 degrees in the x-y plane (rot90) and the Analyze, Interfile,
    MetaImage and raw images are not rotated.

    Parameters
    ----------
    img : ndarray
        3D, 4D or 5D image. If 4D, the fourth dimension is assumed to be the number of iterations and the last iteration is
        used. If 5D, the first four dimensions behave as with a 4D matrix, and the fifth dimension is time (time-series are
        saved if the format supports it).
    outputType : str, optional
        The output file type. 'raw' (default) saves the raw binary data, 'nifti' uses NIfTI, 'analyze' uses Analyze 7.5,
        'dicom' uses DICOM, 'interfile' uses Interfile and 'metaimage' uses MetaImage. Raw, Interfile and MetaImage use the
        data type of the input image. The data type of NIfTI and Analyze is also the same as the input image. DICOM images are
        16-bit.
    filename : str, optional
        The filename (and optionally the path) of the output file, without the file type. If omitted, the name
        [options.machine_name + '_' + options.name] is used (options required).
    options : options object or dict, optional
        The OMEGA options object. Used to save additional information in the output files, e.g. the voxel size (from FOVa_x/Nx,
        FOVa_y/Ny and axial_fov/Nz, i.e. the field of view divided by the size of the image) and the corrections. If omitted,
        the voxel size is 1 mm in x, y and z unless voxelSize is given.
    voxelSize : sequence of 3 floats, optional
        The voxel size (mm) in x, y and z. Overrides the voxel size from the options.
    reko : str, optional
        Name of the reconstruction algorithm. If given, it is added to the description and the filename ('_' + reko) and
        saved in the Interfile header.

    Returns
    -------
    list of str
        The names of the files that were created (for DICOM, one per slice; for formats with a separate header, the header
        file is listed first).
    """
    if outputType is None:
        outputType = 'raw'
    outputType = str(outputType).lower()
    if outputType not in _OUTPUT_TYPES:
        raise ValueError('Unsupported output filetype')
    if filename is None or len(str(filename)) == 0:
        machineName = _get(options, ('machine_name', 'scannerName'), '')
        name = _get(options, 'name', '')
        if not machineName and not name:
            raise ValueError('No output filename found. Use "saveImage(img, outputType, \'outputFileName\')" to set the output filename')
        filename = str(machineName) + '_' + str(name)
    filename = str(filename)
    description = 'OMEGA reconstruction'
    if reko is not None and len(str(reko)) > 0:
        description += ' with ' + str(reko)
        filename += '_' + str(reko)
    else:
        reko = None
    kuva = _lastIteration(img)
    voxel = _voxelSize(kuva.shape, options, voxelSize)
    if outputType == 'nifti':
        _saveNifti(filename, kuva, voxel, description)
        return [filename + '.nii']
    elif outputType == 'analyze':
        _saveAnalyze(filename, kuva, voxel, description)
        return [filename + '.hdr', filename + '.img']
    elif outputType == 'dicom':
        _saveDicom(filename, kuva, voxel, description)
        return [filename + '_slice' + str(k + 1) + '.dcm' for k in range(kuva.shape[2])]
    elif outputType == 'interfile':
        from .saveInterfile import saveInterfile
        saveInterfile(filename, kuva, reko, kuva.dtype, options, voxelSize)
        return [filename + '.h33', filename + '.i33']
    elif outputType == 'metaimage':
        from .saveMetaImage import saveMetaImage
        saveMetaImage(filename, kuva, None, options, voxelSize)
        return [filename + '.mhd', filename + '.raw']
    else:
        dt = _resolveType(None, kuva)
        with open(filename + '.raw', 'wb') as f:
            f.write(_columnMajorBytes(_castImage(kuva, dt), dt, '<'))
        return [filename + '.raw']
