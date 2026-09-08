# -*- coding: utf-8 -*-
"""
Copyright (C) 2024-2025 Ville-Veikko Wettenhovi, Niilo Saarlemo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
from typing import Tuple
import numpy as np
from omegatomo.projector import proj

def CTDetSource(options):
    """
    Computes CT source/detector coordinates if not already input. Projection
    angles, sourceToDetector and sourceToCRot variables have to be input to 
    options. Can also take into account offsets of the detector and/or source.
    Default assumes that the source is directly in front of the center of the
    detector.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    None.

    """
    if options.angles.size == options.nProjections:
        options.angles = options.angles.reshape(1, 1, -1)
    else:
        if options.angles.shape[0] == options.nProjections:
            options.angles = options.angles[:, 0]
        else:
            options.angles = options.angles[0, :]
            options.angles = options.angles.ravel()
    
    R = np.concatenate((np.concatenate((np.cos(options.angles), -np.sin(options.angles)), axis=1), np.concatenate((np.sin(options.angles), np.cos(options.angles)), axis=1)))
    
    if isinstance(options.sourceToCRot, np.ndarray):
        sourceCoordX = -options.sourceToCRot
    else:
        sourceCoordX = np.array(-options.sourceToCRot, dtype=np.float32)
    if isinstance(options.sourceOffsetRow, np.ndarray):
        sourceCoordY = options.sourceOffsetRow
    else:
        sourceCoordY = np.array(options.sourceOffsetRow, dtype=np.float32)
    if isinstance(options.sourceOffsetCol, np.ndarray):
        sourceCoordZ = options.sourceOffsetCol
    else:
        sourceCoordZ = np.array(options.sourceOffsetCol, dtype=np.float32)
    
    if sourceCoordX.size == 1:
        sourceCoordX = np.tile(sourceCoordX, options.nProjections)
    elif sourceCoordX.size != options.nProjections:
        sourceCoordX = np.tile(sourceCoordX, options.nProjections // sourceCoordX.size)
    if sourceCoordZ.size == 1:
        sourceCoordZ = np.tile(sourceCoordZ, options.nProjections)
    elif sourceCoordZ.size != options.nProjections:
        sourceCoordZ = np.tile(sourceCoordZ, options.nProjections // sourceCoordZ.size)
    if sourceCoordY.size == 1:
        sourceCoordY = np.tile(sourceCoordY, options.nProjections)
    elif sourceCoordY.size != options.nProjections:
        sourceCoordY = np.tile(sourceCoordY, options.nProjections // sourceCoordY.size)
    
    testi = np.concatenate((sourceCoordX, sourceCoordY)).T.reshape(2, 1, -1)
    sXY = np.squeeze(np.sum(R * np.transpose(testi, (1, 0, 2)), axis=1)).T
    
    if isinstance(options.sourceToCRot, np.ndarray):
        detCoordX = options.sourceToDetector - options.sourceToCRot
    else:
        detCoordX = np.array(options.sourceToDetector - options.sourceToCRot, dtype=np.float32, ndmin=1)
    detCoordY = 0.
    detCoordZ = 0.
    
    if isinstance(options.detOffsetCol, np.ndarray):
        if options.detOffsetCol.size > 0:
            detCoordZ = options.detOffsetCol
    else:
        detCoordZ = np.array(options.detOffsetCol, dtype=np.float32)
    if isinstance(options.detOffsetRow, np.ndarray):
        if options.detOffsetRow.size > 0:
            detCoordY = options.detOffsetRow
    else:
        detCoordY = np.array(options.detOffsetRow, dtype=np.float32)
    detCoordX = np.reshape(detCoordX, (-1, 1, 1))
    detCoordY = np.reshape(detCoordY, (-1, 1, 1))
    if detCoordX.size != sourceCoordX.size:
        detCoordX = np.tile(detCoordX, sourceCoordX.size // detCoordX.size)
    if detCoordY.size != sourceCoordY.size:
        detCoordY = np.tile(detCoordY, sourceCoordY.size // detCoordY.size)
    if detCoordZ.size != sourceCoordZ.size:
        detCoordZ = np.tile(detCoordZ, sourceCoordZ.size // detCoordZ.size)
    
    XY = np.squeeze(np.sum(R * np.concatenate((detCoordX, detCoordY), axis=1), axis=1)).T
    options.x = np.column_stack((sXY[:, 0], XY[:, 0]))
    options.y = np.column_stack((sXY[:, 1], XY[:, 1]))
    
    if sourceCoordZ.size == 1:
        options.z = np.column_stack((np.tile(sourceCoordZ, detCoordZ.size), np.reshape(detCoordZ, (-1, 1))))
    elif detCoordZ.size == 1:
        options.z = np.column_stack((np.reshape(sourceCoordZ, (-1, 1)), np.tile(detCoordZ, sourceCoordZ.size)))
    else:
        options.z = np.column_stack((sourceCoordZ, np.reshape(detCoordZ, (-1, 1))))
    if options.bedOffset.size > 0:
        options.x = np.tile(options.x, (options.bedOffset.size, 1))
        options.y = np.tile(options.y, (options.bedOffset.size, 1))
        options.z = np.add(np.tile(options.z, (options.bedOffset.size, 1)), np.tile((options.bedOffset - options.bedOffset[-1] / 2), (options.nProjections, 1)))
    if options.angles.size == options.nProjections:
        options.angles = options.angles.reshape(-1, 1)

def setCTCoordinates(options):
    """
    Sets all the necessary coordinates for CT. The functionality depends highly
    on the input parameters. If source/detector coordinates are manually input,
    they are only transformed to the correct format. Otherwise the coordinates
    themselves are computed.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Raises
    ------
    ValueError
        If invalid subset type.

    Returns
    -------
    None.

    """
    if options.subsets > 1 and options.subsetType < 8 and options.subsetType > 0:
        raise ValueError('Subset types < 8 with CT data are not yet supported in Python! Use custom detector coordinates instead or subset types >= 8.')
    if options.x.size == 0 or hasattr(options, 'coord'):
        CTDetSource(options)
        options.coord = True
        if np.size(options.z) / 2 > np.size(options.angles):
            if options.angles.shape[0] == 1:
                options.angles = np.reshape(options.angles, (-1, 1, 1))
            options.angles = np.tile(options.angles, (len(options.z) // 2 // np.size(options.angles), 1, 1))
    elif abs(options.offangle) > 0:
        R = np.array([[np.cos(options.offangle), -np.sin(options.offangle)],
                      [np.sin(options.offangle), np.cos(options.offangle)]])
        det = np.column_stack((options.x[:, 0], options.y[:, 0])).T
        sXY = np.dot(R, det).T
        det = np.column_stack((options.x[:, 1], options.y[:, 1])).T
        dXY = np.dot(R, det).T
        options.x = np.column_stack((sXY[:, 0], dXY[:, 0]))
        options.y = np.column_stack((sXY[:, 1], dXY[:, 1]))
    if options.x.size != options.nColsD * options.nRowsD * options.nProjections and options.uV.size <= 1 and not options.useHelical:
        options.uV = CTDetectorCoordinates(options.angles,options.pitchRoll)
        if np.ndim(options.uV) == 3 and options.uV.shape[2] > 1:
            options.uV = np.reshape(options.uV, (-1, 2))
    if options.flip_image:
        options.y = -options.y
        options.uV[:,1] = -options.uV[:,1]
        if options.pitchRoll.size > 0:
            options.uV[:,4] = -options.uV[:,4]
    if options.uV.size > 0:
        options.x = np.column_stack((options.x[:,0], options.y[:,0], options.z[:,0], options.x[:,1], options.y[:,1], options.z[:,1]))
        options.z = options.uV
    if options.pitchRoll.size > 0:
        options.pitch = True
        options.z[:,0] = options.z[:,0] * options.dPitchX
        options.z[:,1] = options.z[:,1] * options.dPitchX
        options.z[:,5] = options.z[:,5] * options.dPitchY
        options.z[:,3] = options.z[:,3] * options.dPitchX
        options.z[:,4] = options.z[:,4] * options.dPitchX
        options.z[:,2] = options.z[:,2] * options.dPitchY
    elif options.uV.size > 0:
        options.z[:,0] = options.z[:,0] * options.dPitchX
        options.z[:,1] = options.z[:,1] * options.dPitchX
    if options.useHelical:
        options.x = np.column_stack((options.x[:,0], options.y[:,0], options.z[:,0], options.x[:,1], options.y[:,1], options.z[:,1]))
        options.z = np.float32(options.angles)

def CTDetectorCoordinates(angles, pitchRoll = np.empty(0, dtype=np.float32)):
    """
    Computes the direction vectors for each projection based on the input data.

    Parameters
    ----------
    angles : NumPy array
        Projection angles.
    pitchRoll : NumPy array, optional
        The pitch/roll/yaw angles. The default is np.empty(0, dtype=np.float32).

    Returns
    -------
    uV : NumPy array
        The direction vectors of the panel pixels for each projection.

    """
    if pitchRoll.size == 0:
        uV = np.column_stack((-np.sin(angles), np.cos(angles)))
    else:
        pitchRoll.reshape(pitchRoll.size // 2, 2)
        uV = np.column_stack(((-np.sin(angles) * np.cos(pitchRoll[:,0] - np.cos(angles) * np.sin(pitchRoll[:,0]) * np.sin(pitchRoll[:,1]))), 
                             (np.cos(angles) * np.cos(pitchRoll[:,0]) - np.cos(angles) * np.sin(pitchRoll[:,0]) * np.sin(pitchRoll[:,1])),
                             (np.sin(pitchRoll[:,0]) * np.cos(pitchRoll[:,1])),
                             (np.sin(angles) * np.sin(pitchRoll[:,0]) - np.cos(angles) * np.cos(pitchRoll[:,0]) * np.sin(pitchRoll[:,1])),
                             (-(np.cos(angles) * np.sin(pitchRoll[:,0]) + np.sin(angles) * np.cos(pitchRoll[:,0]) * np.sin(pitchRoll[:,1]))),
                             (np.cos(pitchRoll[:,1]) * np.cos(pitchRoll[:,0]))))
    return uV


def getCoordinates(options):   
    """
    Computes PET sinogram coordinates if not already input.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    x : NumPy array
        x-direction coordinates for detector/detector pairs.
    y : NumPy array
        y-direction coordinates for detector/detector pairs.
    z : NumPy array
        z-direction coordinates for detector/detector pairs.

    """
    if isinstance(options.x, np.ndarray) and options.x.size > 1 and options.y.size > 1 and options.listmode == 0:
        if not(options.z.size == options.x.size) and options.z.size < options.nProjections:
            options.z = np.zeros(options.x.shape,dtype=np.float32)
        x = options.x
        y = options.y
        z = options.z
    else:
        if options.use_raw_data == 0:
            x,y = detectorCoordinates(options)
            if options.nLayers > 1:
                koko = x.size / 2
            else:
                x, y = sinogramCoordinates2D(options, x, y)
                
            if options.arc_correction:
                from omegatomo.util.arcCorrection import arc_correction
                x, y, options = arc_correction(options, False)
            if options.sampling > 1:
                from omegatomo.util.sampling import increaseSampling
                x, y, options = increaseSampling(options, x, y, False)
    
            # if options.arc_correction and ~options.precompute_lor
            #     [x, y, options] = arcCorrection(options, xp, yp, interpolateSinogram);
            # if options.sampling > 1 and ~options.precompute_lor
            #     [x, y, options] = increaseSampling(options, x, y, interpolateSinogram);
            z = sinogramCoordinates3D(options)
            if not options.NSinos == options.TotSinos:
                z = z[0:options.NSinos,:]
            x = np.column_stack((x[:,0], y[:,0], x[:,1], y[:,1]))
        else:
            if options.det_per_ring < options.det_w_pseudo:
                options.offangle = options.offangle / options.det_w_pseudo
                options.offangle = options.offangle * options.det_per_ring
            x, y = detectorCoordinates(options)
            # if options.sampling_raw > 1 and ~options.precompute_lor
            #     [x, y, options] = increaseSampling(options, x, y, interpolateSinogram);
            z_length = float(options.rings + np.sum(options.pseudot)) * options.cr_pz
        
            z = np.linspace(-(z_length / 2 - options.cr_pz / 2), z_length / 2 - options.cr_pz / 2, options.rings + int(np.sum(options.pseudot)),dtype=np.float32)
        
            if np.sum(options.pseudot) > 0:
                z = np.delete(z, np.where(options.pseudot))
        
            x = np.asfortranarray(np.vstack((x, y)))
    return x, y, z
    
def getCoordinatesSPECT(options: proj.projectorClass) -> Tuple[np.ndarray, np.ndarray]: # TODO this function
    n1: int = len(options.angles)
    n2: int = len(options.radiusPerProj)
    n3: int = len(options.swivelAngles)
    assert (n1 == n2) and (n2 == n3), "The amount of angles, radii and swivel angles have to be equal."

    nProjections: int = n1

    x: np.ndarray = np.zeros((6, nProjections), dtype=np.float32)
    z: np.ndarray = np.zeros((2, nProjections), dtype=np.float32)

    for ii in range(nProjections):
        r1: float = options.radiusPerProj[ii]
        r2: float = options.CORtoDetectorSurface

        alpha1: float = options.angles[ii]
        alpha2: float = options.swivelAngles[ii]

        x[0, ii] = r1 * np.cos(np.deg2rad(alpha1)) + r2 * np.cos(np.deg2rad(alpha2))
        x[1, ii] = r1 * np.sin(np.deg2rad(alpha1)) + r2 * np.sin(np.deg2rad(alpha2))
        x[2, ii] = 0

        x[3, ii] = x[0, ii] + (options.colD + 0.5 * options.colLxy) * np.cos(np.deg2rad(alpha2))
        x[4, ii] = x[1, ii] + (options.colD + 0.5 * options.colLxy) * np.sin(np.deg2rad(alpha2))
        x[5, ii] = 0

        z[0, ii] = np.cos(np.deg2rad(alpha2 + 90))
        z[1, ii] = np.sin(np.deg2rad(alpha2 + 90))

    if options.flipImageX:  # Horizontal
        x[0, :] = -x[0, :]
        x[3, :] = -x[3, :]
        z[0, :] = -z[0, :]

    if options.flipImageY: # Vertical
        x[1, :] = -x[1, :]
        x[4, :] = -x[4, :]
        z[1, :] = -z[1, :]

    return np.asfortranarray(x), np.asfortranarray(z)
    
def detectorCoordinates(options):
    """
    Transaxial PET coordinates. This function mainly just calls the below one 
    with specific options. Also, the detector space rotation is done here if
    selected.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    x : NumPy array
        x-direction coordinates for detector/detector pairs.
    y : NumPy array
        y-direction coordinates for detector/detector pairs.

    """
    cr_p = options.cr_p
    diameter = options.diameter
    if isinstance(options.cryst_per_block, np.ndarray):
        cryst_per_block = options.cryst_per_block[0].item()
    else:
        cryst_per_block = options.cryst_per_block
    blocks_per_ring = options.blocks_per_ring
    x = 0
    y = 0
    
    DOI = options.DOI
    
    transaxial_multip = options.transaxial_multip
    
    diameter = diameter + DOI * 2
    
    
    if options.det_w_pseudo > options.det_per_ring:
        
        xp,yp = computeCoordinates(options, blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, True)
        x = xp
        y = yp
        
    elif options.nLayers > 1:
        if options.cryst_per_block[0] == options.cryst_per_block[1]:
            x1,y1 = computeCoordinates(options, blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, False)
            # orig_diameter = diameter + options.crystH(end) * 2;
            diameter = diameter + DOI * 2 + options.crystH[0] * 2
            x2,y2 = computeCoordinates(options, blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, False)
            # xp = [x;xp];
            # yp = [y;yp];
        elif options.cryst_per_block[0] > options.cryst_per_block[1]:
            x,y = computeCoordinates(options, blocks_per_ring, transaxial_multip,options.cryst_per_block[-1].item(), diameter, cr_p, True)
            # orig_diameter = diameter + options.crystH(end) * 2;
            # diameter = diameter + DOI * 2 + options.crystH(end) * 2;
            # [x,y] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, false);
            # xp = [x;xp];
            # yp = [y;yp];
        else:
            x,y = computeCoordinates(options, blocks_per_ring, transaxial_multip,options.cryst_per_block[-1].item(), diameter, cr_p, False)
            # orig_diameter = diameter + options.crystH(end) * 2;
            # diameter = diameter + DOI * 2 + options.crystH(end) * 2;
            # [xp,yp] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, true);
            # xp = [xp;x];
            # yp = [yp;y];
    else:
        x,y = computeCoordinates(options, blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, False)
    
    if options.flip_image:
        x = np.flip(x)
        y = np.flip(y)
    if not options.offangle == 0:
        x = np.roll(x, int(np.round(options.offangle)))
        y = np.roll(y, int(np.round(options.offangle)))
    return x, y

def computeCoordinates(options, blocks_per_ring, transaxial_multip, cryst_per_block, diameter, cr_p, usePseudo):
    """
    Transaxial PET coordinates.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    blocks_per_ring : int
        Number of blocks per transaxial ring.
    transaxial_multip : int
        If sub-blocks are present the number of those.
    cryst_per_block : int
        The number of crystals per block.
    diameter : float
        Diameter of the bore.
    cr_p : float
        Crystal pitch/size in transaxial direction.
    usePseudo : bool
        True if pseudo detectors are present, False otherwise.

    Returns
    -------
    x : NumPy array
        x-direction coordinates for detector/detector pairs.
    y : NumPy array
        y-direction coordinates for detector/detector pairs.

    """
    import math
    cryst_per_block_orig = cryst_per_block
    if usePseudo:
        extraVar = 1
        cryst_per_block = cryst_per_block + 1
    else:
        extraVar = 0
    
    angle = np.arange(90, -271, -360 / (blocks_per_ring * transaxial_multip))  # Generate angles
    
    widths = (cryst_per_block_orig) * cr_p * np.cos(np.radians(angle[:((blocks_per_ring * transaxial_multip) // 2 + 1)]))
    widthsy = (cryst_per_block_orig) * cr_p * np.sin(np.radians(angle[:((blocks_per_ring * transaxial_multip) // 2)]))
    
    erotus = (diameter - np.sum(np.abs(widths))) / np.sum(np.cos(np.radians(angle[:((blocks_per_ring * transaxial_multip) // 2 + 1)]))) / 2
    erotusy = (diameter - np.sum(np.abs(widthsy))) / np.sum(np.abs(np.sin(np.radians(angle[:((blocks_per_ring * transaxial_multip) // 2)])))) / 2
    
    alkupistex = diameter / 2.
    alkupistey = -((cryst_per_block_orig) / 2. + 0.5) * cr_p
    
    ii = 0
    x = np.zeros(((blocks_per_ring * transaxial_multip) * (cryst_per_block)), dtype=np.float32)
    y = np.zeros(((blocks_per_ring * transaxial_multip) * (cryst_per_block)), dtype=np.float32)
    
    # Compute the detector coordinates of each detector (crystal) in each block
    # Only for the 1/4th of the ring
    for blocks in range(0, math.ceil((blocks_per_ring * transaxial_multip) / 4)):
        for crystals in range(1, cryst_per_block + 1):
            if blocks > 0 and crystals == 1:
                x[ii] = alkupistex - (cr_p * 0.5) * np.cos(np.radians(angle[blocks])) - erotus * np.cos(
                    np.radians(angle[blocks - 1])) - erotus * np.cos(np.radians(angle[blocks])) - (
                                cr_p * 0.5) * np.cos(np.radians(angle[blocks - 1]))
                y[ii] = alkupistey + (cr_p * 0.5) * np.sin(np.radians(angle[blocks])) + erotusy * np.sin(
                    np.radians(angle[blocks - 1])) + erotusy * np.sin(np.radians(angle[blocks])) + (
                                cr_p * 0.5) * np.sin(np.radians(angle[blocks - 1]))
            else:
                x[ii] = alkupistex - (cr_p) * np.cos(np.radians(angle[blocks]))
                y[ii] = alkupistey + (cr_p) * np.sin(np.radians(angle[blocks]))
    
            if crystals == cryst_per_block:
                alkupistex = x[ii] + (cr_p) * np.cos(np.radians(angle[blocks])) * extraVar
                alkupistey = y[ii] - (cr_p) * np.sin(np.radians(angle[blocks])) * extraVar
            else:
                alkupistex = x[ii]
                alkupistey = y[ii]
            ii += 1
    
    # Symmetry of the coordinates
    if (blocks_per_ring * transaxial_multip) % 4 == 0:
        blocks += 1
        x[ii] = alkupistex - (cr_p * 0.5) * np.cos(np.radians(angle[blocks])) - erotus * np.cos(
            np.radians(angle[blocks - 1])) - erotus * np.cos(np.radians(angle[blocks])) - (cr_p * 0.5) * np.cos(
            np.radians(angle[blocks - 1]))
        y[ii] = alkupistey + (cr_p * 0.5) * np.sin(np.radians(angle[blocks])) + erotusy * np.sin(
            np.radians(angle[blocks - 1])) + erotusy * np.sin(np.radians(angle[blocks])) + (cr_p * 0.5) * np.sin(
            np.radians(angle[blocks - 1]))
        alkupistex = x[ii]
        alkupistey = y[ii]
        ii += 1
    
        for ll in range(2, cryst_per_block + 1):
            x[ii] = alkupistex - (cr_p) * np.cos(np.radians(angle[blocks]))
            y[ii] = alkupistey + (cr_p) * np.sin(np.radians(angle[blocks]))
            alkupistex = x[ii]
            alkupistey = y[ii]
            ii += 1
    
        if usePseudo:
            ii -= 1
            alkupistex = x[ii] + (cr_p) * np.cos(np.radians(angle[blocks])) * extraVar
            alkupistey = y[ii] - (cr_p) * np.sin(np.radians(angle[blocks])) * extraVar
            alkublock = blocks + 1
            ii += 1
            for blocks in range(alkublock, blocks_per_ring * transaxial_multip + 1):
                for crystals in range(1, cryst_per_block + 1):
                    if blocks > 0 and crystals == 1:
                        x[ii] = alkupistex - (cr_p * 0.5) * np.cos(np.radians(angle[blocks])) - erotus * np.cos(
                            np.radians(angle[blocks - 1])) - erotus * np.cos(np.radians(angle[blocks])) - (
                                        cr_p * 0.5) * np.cos(np.radians(angle[blocks - 1]))
                        y[ii] = alkupistey + (cr_p * 0.5) * np.sin(np.radians(angle[blocks])) + erotusy * np.sin(
                            np.radians(angle[blocks - 1])) + erotusy * np.sin(np.radians(angle[blocks])) + (
                                        cr_p * 0.5) * np.sin(np.radians(angle[blocks - 1]))
                    else:
                        x[ii] = alkupistex - (cr_p) * np.cos(np.radians(angle[blocks]))
                        y[ii] = alkupistey + (cr_p) * np.sin(np.radians(angle[blocks]))
    
                    if crystals == cryst_per_block:
                        alkupistex = x[ii] + (cr_p) * np.cos(np.radians(angle[blocks])) * extraVar
                        alkupistey = y[ii] - (cr_p) * np.sin(np.radians(angle[blocks])) * extraVar
                    else:
                        alkupistex = x[ii]
                        alkupistey = y[ii]
                    ii += 1
        else:
            x[ii:ii + ((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4] = -(np.flip(x[:((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4]))
            y[ii:ii + ((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4] = (np.flip(y[:((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4]))
            x[ii + ((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4:] = np.flip(x[cryst_per_block:((blocks_per_ring * transaxial_multip) * cryst_per_block) // 2])
            y[ii + ((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4:] = -(y[cryst_per_block:((blocks_per_ring * transaxial_multip) * cryst_per_block) // 2])
    else:
        if usePseudo:
            alkublock = blocks + 1
            for blocks in range(alkublock, blocks_per_ring * transaxial_multip + 1):
                for crystals in range(1, cryst_per_block + 1):
                    if blocks > 0 and crystals == 1:
                        x[ii] = alkupistex - (cr_p * 0.5) * np.cos(np.radians(angle[blocks])) - erotus * np.cos(
                            np.radians(angle[blocks - 1])) - erotus * np.cos(np.radians(angle[blocks])) - (
                                        cr_p * 0.5) * np.cos(np.radians(angle[blocks - 1]))
                        y[ii] = alkupistey + (cr_p * 0.5) * np.sin(np.radians(angle[blocks])) + erotusy * np.sin(
                            np.radians(angle[blocks - 1])) + erotusy * np.sin(np.radians(angle[blocks])) + (
                                        cr_p * 0.5) * np.sin(np.radians(angle[blocks - 1]))
                    else:
                        x[ii] = alkupistex - (cr_p) * np.cos(np.radians(angle[blocks]))
                        y[ii] = alkupistey + (cr_p) * np.sin(np.radians(angle[blocks]))
    
                    if crystals == cryst_per_block:
                        alkupistex = x[ii] + (cr_p) * np.cos(np.radians(angle[blocks])) * extraVar
                        alkupistey = y[ii] - (cr_p) * np.sin(np.radians(angle[blocks])) * extraVar
                    else:
                        alkupistex = x[ii]
                        alkupistey = y[ii]
                    ii += 1
        else:
            x[ii:ii + ((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4 + cryst_per_block // 2] = -(np.flip(x[:((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4 + cryst_per_block // 2]))
            y[ii:ii + ((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4 + cryst_per_block // 2] = (np.flip(y[:((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4 + cryst_per_block // 2]))
            x[ii + ((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4 + cryst_per_block // 2:] = np.flip(x[cryst_per_block:((blocks_per_ring * transaxial_multip) * cryst_per_block) // 2])
            y[ii + ((blocks_per_ring * transaxial_multip) * cryst_per_block) // 4 + cryst_per_block // 2:] = -(y[cryst_per_block:((blocks_per_ring * transaxial_multip) * cryst_per_block) // 2])
    return x, y
            
def formDetectorIndices( det_w_pseudo, nLayers = 1, crystN = 0):
    """
    For raw data only. Note used at the moment.

    Parameters
    ----------
    det_w_pseudo : int
        Number of detectors per transaxial ring, including pseudos.
    nLayers : int, optional
        Number of crystal layers. The default is 1.
    crystN : int, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    L : NumPy array
        Detector indices.

    """
    L = np.zeros((det_w_pseudo*(det_w_pseudo+1)//2, 2), dtype=np.int32)
    jh = 0
    for kk in range(1, det_w_pseudo + 1):
        L[jh:(jh + det_w_pseudo - kk + 1), :] = np.column_stack((np.repeat(kk, det_w_pseudo - (kk - 1)), np.arange(kk, det_w_pseudo + 1)))
        jh += det_w_pseudo - kk + 1
    L = L[~np.all(L == 0, axis=1)]
    
    if nLayers > 1:
        temp = np.arange(crystN, det_w_pseudo + 1, crystN)
        L = L[~np.isin(L[:, 0], temp)]
        L = L[~np.isin(L[:, 1], temp)]
    return L
        
def sinogramCoordinates2D(options, x, y, nLayers = 1):
    """
    Computes the transaxial sinogram coordinates using the input detector
    coordinates.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    x : NumPy array
        x-direction coordinates for detector/detector pairs.
    y : NumPy array
        y-direction coordinates for detector/detector pairs.
    nLayers : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    x : NumPy array
        x-direction coordinates for sinogram bins.
    y : NumPy array
        y-direction coordinates for sinogram bins.

    """
    det_w_pseudo = options.det_w_pseudo
    Nang = options.Nang
    Ndist = options.Ndist
    mashing = 1
    # Determine the possible mashing factor
    if Nang < options.det_w_pseudo//2:
        mashing = options.det_w_pseudo // Nang // 2
        Nang = Nang * mashing
    
    
    ## 2D coordinates
    
    # Determine the sinogram indices for each of the LOR
    
    # Form the detector vector pair
    if isinstance(options.cryst_per_block, np.ndarray):
        L = formDetectorIndices(det_w_pseudo, nLayers, options.cryst_per_block[0].item())
    else:
        L = formDetectorIndices(det_w_pseudo, nLayers, options.cryst_per_block)
    
    L = L - 1
    
    xa = np.max(L,1)
    ya = np.min(L,1)
    
    # Angle
    j = ((xa + ya + det_w_pseudo //2 + 1) % det_w_pseudo) // 2
    
    b = j + det_w_pseudo // 2;
    
    # Distance
    i = np.abs(xa - ya - det_w_pseudo // 2)
    for kk in range(len(ya)):
        if (ya[kk] < j[kk]) or (b[kk] < xa[kk]):
            i[kk] = -i[kk]
    
    # The sinogram corners need to the swapped
    swap = np.logical_or((j * 2) < -i, i <= ((j - det_w_pseudo // 2) * 2))
    L3 = L[swap,0]
    L[swap,0] = L[swap,1]
    L[swap,1] = L3
    
    # Determine the accepted LORs (distances that are within the predefined
    # value)
    if Ndist % 2 == 0:
        accepted_lors = np.logical_and(i <= (Ndist//2 + min(0,options.ndist_side)), i >= (-Ndist//2 + max(0,options.ndist_side)))
    else:
        accepted_lors = np.logical_and(i <= Ndist//2, i >= -(Ndist//2))
    
    j = j // (det_w_pseudo // 2 // Nang)
    
    i = i[accepted_lors]
    j = j[accepted_lors]
    if np.min(i) < 0:
        i = i + np.abs(np.min(i))
    
    L = L[accepted_lors,:]
    
    xx1 = x[L[:,0]]
    yy1 = y[L[:,0]]
    xx2 = x[L[:,1]]
    yy2 = y[L[:,1]]
    
    ##
    
    x = np.zeros((Ndist, Nang),dtype=np.float32,order='F')
    y = np.zeros((Ndist, Nang),dtype=np.float32,order='F')
    x2 = np.zeros((Ndist, Nang),dtype=np.float32,order='F')
    y2 = np.zeros((Ndist, Nang),dtype=np.float32,order='F')
    np.add.at(x, (i, j), xx1)
    np.add.at(y, (i, j), yy1)
    np.add.at(x2, (i, j), xx2)
    np.add.at(y2, (i, j), yy2)

    
    # If mashing is present, combine the coordinates
    if mashing > 1:
        from skimage.measure import block_reduce
        # Compute the mean coordinates
        x = block_reduce(x, block_size=(1,mashing), func=np.mean)
        y = block_reduce(y, block_size=(1,mashing), func=np.mean)
        x2 = block_reduce(x2, block_size=(1,mashing), func=np.mean)
        y2 = block_reduce(y2, block_size=(1,mashing), func=np.mean)
        
        j = j // mashing
    
    x = np.column_stack((x.ravel('F'), x2.ravel('F')))
    y = np.column_stack((y.ravel('F'), y2.ravel('F')))
    return x, y
    
def sinogramCoordinates3D(options, layers = (1,1)):
    """
    Same as above, but for axial coordinates.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    layers : tuple, optional
        The crystal layer combinations. The default is (1,1).

    Returns
    -------
    z : Numpy array
        Axial sinogram coordinates.

    """
    import math
    layer1 = layers[0]
    layer2 = layers[1]
    cr_pz = options.cr_pz
    Nz = options.rings * 2 - 1
    
    if len(options.ringGaps) > 0 and np.sum(options.ringGaps) > 0:
        z_length = float(options.rings + 1) * cr_pz + np.sum(options.ringGaps)
    else:
        z_length = float(options.rings + 1) * cr_pz
    
    if options.nLayers > 1 and options.cryst_per_block_axial[0] != options.cryst_per_block_axial[1]:
        maxZ = z_length + cr_pz * (options.linear_multip - 1)
    else:
        maxZ = z_length
    
    if options.nLayers > 1 and options.cryst_per_block_axial[0] != options.cryst_per_block_axial[1]:
        apu = np.zeros(options.linear_multip)
        for kk in range(options.linear_multip):
            apu[kk] = kk * cr_pz
        apu = np.repeat(apu, np.sum(options.cryst_per_block_axial))
    
    # Compute the 3D coordinates
    if options.span > 1:
        z = np.linspace(cr_pz, z_length + cr_pz, options.rings + 2)
        z = z[:options.rings]
        # if len(options.ringGaps) > 0 and np.sum(options.ringGaps) > 0:
        #     z = np.linspace(cr_pz, z_length + cr_pz, options.rings + 2 + len(options.ringGaps))
        # else:
        #     z = np.linspace(cr_pz, z_length + cr_pz, options.rings + 2)
        
        if len(options.ringGaps) > 0 and np.sum(options.ringGaps) > 0:
            gaps = np.cumsum(options.ringGaps)
            for kk in range(1, int(options.rings / options.cryst_per_block_axial)):
                start_idx = options.cryst_per_block_axial * kk
                end_idx = options.cryst_per_block_axial * (kk + 1)
                z[start_idx:end_idx] = z[start_idx:end_idx] + gaps[kk - 1]
            # z = z[:options.rings + len(options.ringGaps)]
            # z = np.delete(z, options.ringGaps + 1)
        # else:
        #     z = z[:options.rings]
        
        if options.nLayers > 1 and options.cryst_per_block_axial[0] != options.cryst_per_block_axial[1]:
            z += apu
        z -= maxZ / 2
        ringsp = options.rings
        z_ring = np.zeros((options.rings, options.rings, 2))
        
        z = np.reshape(z, (-1, 1))
        # Create ring combinations
        z_ring[:,:,0] = np.tile(z, (1, options.rings))
        z_ring[:,:,1] = np.tile(z.T, (options.rings, 1))
        z_ring = z_ring.reshape(options.rings * options.rings, 2)
        kkj = np.zeros(((options.ring_difference - math.ceil(options.span / 2)) // options.span) + 1, dtype=int)
    
        for kk in range(1, ((options.ring_difference - math.ceil(options.span / 2)) // options.span) + 2):
            kkj[kk - 1] = np.ceil(options.span / 2) + options.span * (kk - 1)
        
        offset2 = np.cumsum(options.segment_table)
        
        # Perpendicular rings
        z = np.zeros((options.TotSinos, 2))
        z[:Nz+1:2, 0] = z_ring[::ringsp + 1, 0]
        z[:Nz+1:2, 1] = z_ring[::ringsp + 1, 1]
        mean_jh = np.zeros(options.TotSinos)
        mean_jh[:Nz + 1:2] = 1
        
        # Then the detectors on adjacent rings
        for jh in range(0, int(np.floor(options.span / 2))):
            apu = z_ring[((jh + 1) * ringsp)::ringsp + 1, 0]
            apu2 = z_ring[(jh + 1):((ringsp - jh + 1) * ringsp):ringsp + 1, 0]
            z[jh + 1:offset2[0] - jh:2, 0] = z[jh + 1:offset2[0] - jh:2, 0] + apu + apu2
            apu = z_ring[((jh + 1) * ringsp)::ringsp + 1, 1]
            apu2 = z_ring[(jh + 1):(ringsp - jh) * ringsp:ringsp + 1, 1]
            z[jh + 1:offset2[0] - jh:2, 1] = z[jh + 1:offset2[0] - jh:2, 1] + apu + apu2
            loc2 = np.isin(np.arange(0, Nz), np.arange(jh + 1,offset2[0] - jh, 2))
            loc = np.full(len(mean_jh), False)
            loc[:len(loc2)] = loc2
            mean_jh[loc] += 2
        
        # Lastly the rest of the detectors with the amount of combined LORs
        # specified with the span value
        for ih in range(1, int(len(options.segment_table) / 2) + 1):
            for jh in range(1, options.span + 1):
                apu = z_ring[((kkj[ih - 1] + jh - 1) * options.rings)::options.rings + 1, 0]
                z[offset2[2 * (ih - 1)] + jh - 1:offset2[(2 * ih) - 1] - jh + 1:2, 0] += apu
                apu2 = z_ring[(kkj[ih - 1] + jh - 1):((options.rings - kkj[ih - 1] - jh + 2) * options.rings):options.rings + 1, 0]
                z[offset2[(2 * ih - 1)] + jh - 1:offset2[(2 * ih)] - jh + 1:2, 0] += apu2
                apu = z_ring[((kkj[ih - 1] + jh - 1) * options.rings)::(options.rings + 1), 1]
                z[offset2[2 * (ih - 1)] + jh - 1:offset2[(2 * ih) - 1] - jh + 1:2, 1] += apu
                apu2 = z_ring[(kkj[ih - 1] + jh - 1):((options.rings - kkj[ih - 1] - jh + 1) * options.rings):options.rings + 1, 1]
                z[offset2[(2 * ih - 1)] + jh - 1:offset2[(2 * ih)] - jh + 1:2, 1] += apu2
                loc = np.isin(np.arange(1, options.TotSinos + 1), np.arange(offset2[(2 * (ih - 1))] + jh, offset2[(2 * ih - 1)] - jh + 2, 2))
                mean_jh[loc] += 1
                loc = np.isin(np.arange(1, options.TotSinos + 1), np.arange(offset2[(2 * ih - 1)] + jh, offset2[(2 * ih)] - jh + 2, 2))
                mean_jh[loc] += 1
        
        mean_jh[mean_jh == 0] = 1
        # Take the mean value of coordinates
        z[:, 0] /= mean_jh
        z[:, 1] /= mean_jh
        # z = np.fliplr(z)
        ind1 = z[:, 0] < z[:, 1]
        # z[ind1] = np.fliplr(z[ind1])
        # z[ind1] = np.fliplr(z[ind1])
    else:
        dif = cr_pz

        if options.nLayers > 1 and options.cryst_per_block_axial[0] != options.cryst_per_block_axial[1]:
            gap_start = ((float(options.rings) * cr_pz) - (options.cryst_per_block_axial[0] * options.linear_multip) * cr_pz) / options.linear_multip / 2
            gap = np.arange(gap_start, gap_start * 2 * options.linear_multip + 1e-6, gap_start * 2)

        z = np.zeros((options.rings**2, 2))
        loppu = options.rings
        ind1 = np.ones((loppu, 1))
        
        if options.nLayers > 1:
            if options.cryst_per_block_axial[0] > options.cryst_per_block_axial[-1] and layer2 == 1:
                r = options.cryst_per_block_axial[-1] * options.linear_multip
                apu = np.repeat(gap, options.cryst_per_block_axial[-1]).reshape(-1, 1)
                ind2 = np.arange(0, cr_pz * r, cr_pz).reshape(-1, 1) + apu
                insert_indices = np.setdiff1d(np.arange(len(ind2)), np.arange(options.cryst_per_block[0]-1, len(ind2), options.cryst_per_block[0]))
                ind3 = np.full((loppu, 1), np.inf)
                ind3[insert_indices] = ind2[insert_indices]
                ind2 = ind3
            elif options.cryst_per_block_axial[-1] > options.cryst_per_block_axial[0] and layer2 == 1:
                r = options.cryst_per_block_axial[0] * options.linear_multip
                apu = np.repeat(gap, options.cryst_per_block_axial[0]).reshape(-1, 1)
                ind2 = np.arange(0, cr_pz * r, cr_pz).reshape(-1, 1) + apu
                ind3 = np.full((loppu, 1), np.inf)
                insert_indices = np.setdiff1d(np.arange(len(ind3)), np.arange(options.cryst_per_block[-1]-1, len(ind3), options.cryst_per_block[-1]))
                ind3[insert_indices] = ind2[insert_indices]
                ind2 = ind3
            else:
                ind2 = np.arange(0, cr_pz * loppu, cr_pz).reshape(-1, 1)
        
            if (options.cryst_per_block_axial[0] > options.cryst_per_block_axial[-1] or options.cryst_per_block_axial[-1] > options.cryst_per_block_axial[0]) and layer1 == 1:
                apu1 = np.full((loppu, 1), np.inf)
                if options.cryst_per_block_axial[-1] > options.cryst_per_block_axial[0]:
                    apu = np.repeat(gap, options.cryst_per_block_axial[0]).reshape(-1, 1)
                    insert_indices = np.setdiff1d(np.arange(len(ind1)), np.arange(options.cryst_per_block[-1]-1, len(ind1), options.cryst_per_block[-1]))
                else:
                    apu = np.repeat(gap, options.cryst_per_block_axial[-1]).reshape(-1, 1)
                    insert_indices = np.setdiff1d(np.arange(len(ind1)), np.arange(options.cryst_per_block[0]-1, len(ind1), options.cryst_per_block[0]))
                apu1[insert_indices] = apu[insert_indices]
        else:
            ind2 = np.arange(0, cr_pz * loppu, cr_pz).reshape(-1, 1)
        
        uu = 0
        if hasattr(options, 'ringGaps') and np.sum(options.ringGaps) > 0:
            gaps = np.repeat(np.insert(np.cumsum(options.ringGaps), 0, 0), options.cryst_per_block_axial).reshape(-1, 1)
        
        yy = 0
        for t in range(1, loppu + 1):
            idx_start = (t - 1) * loppu
            idx_end = t * loppu
            if options.nLayers > 1 and options.cryst_per_block_axial[0] != options.cryst_per_block_axial[1] and layer1 == 1:
                z[idx_start:idx_end, :] = np.hstack([dif * uu * ind1 + apu1[t - 1], ind2])
                uu += 1
            else:
                z[idx_start:idx_end, :] = np.hstack([dif * (t - 1) * ind1, ind2])
                if hasattr(options, 'ringGaps') and np.sum(options.ringGaps) > 0:
                    z[idx_start:idx_end, 1] += np.squeeze(gaps[:loppu])
                    if (t % (options.cryst_per_block_axial + 1)) == 0:
                        yy += 1
                    if yy > 0:
                        z[idx_start:idx_end, 0] += yy * options.ringGaps[yy - 1]
        
            if options.nLayers > 1 and ((options.cryst_per_block_axial[0] > options.cryst_per_block_axial[-1] and layer1 == 1) or (options.cryst_per_block_axial[-1] > options.cryst_per_block_axial[0] and layer1 == 1)):
                if np.isinf(apu1[t - 1]):
                    z[idx_start:idx_end, :] = np.hstack([np.inf * ind1, ind2])
                    uu -= 1
        
        z[np.isnan(z)] = np.inf
        z = z - (maxZ / 2 - cr_pz)
    return z


def _type6_scalar(values, volume: int) -> float:
    """Read a scalar image-geometry value for one type-6 volume."""
    values = np.asarray(values).reshape(-1)
    if values.size == 0:
        raise ValueError('Missing type-6 image geometry')
    return float(values[min(volume, values.size - 1)])


def _type6_auto_filter(options: proj.projectorClass, volume: int, depth: int | None = None) -> np.ndarray:
    """Generate one rotation-projector CDRF on a volume's pixel grid."""
    nx = int(np.asarray(options.Nx).reshape(-1)[volume])
    ny = int(np.asarray(options.Ny).reshape(-1)[volume])
    nz = int(np.asarray(options.Nz).reshape(-1)[volume])
    dx = _type6_scalar(options.dx, volume)
    dy = _type6_scalar(options.dy, volume)
    dz = _type6_scalar(options.dz, volume)
    if min(dx, dy, dz) <= 0.:
        raise ValueError(f'Type-6 volume {volume} has a non-positive voxel pitch')

    sigma_z_input = np.asarray(options.sigmaZ, dtype=np.float32)
    sigma_xy_input = np.asarray(options.sigmaXY, dtype=np.float32)
    if sigma_z_input.size != 1 or sigma_xy_input.size != 1:
        raise ValueError(
            'Multi-resolution projector type 6 currently requires scalar '
            'sigmaZ and sigmaXY, or an explicit gFilter entry for each volume.'
        )

    depth = max(nx * 4, int(depth or 0))
    if float(sigma_z_input.reshape(-1)[0]) < 0.:
        distances = 0.5 * dx + np.arange(depth, dtype=np.float32) * dx
        col_l_xy = float(options.colLxy) if options.colLxy > 0. else float(options.colL)
        col_l_z = float(options.colLz) if options.colLz > 0. else float(options.colL)
        distance_from_exit = distances + float(options.cr_p) / 2.
        rg_z = np.maximum(0., 2. * float(options.colR) * distance_from_exit / col_l_z)
        rg_xy = np.maximum(0., 2. * float(options.colR) * distance_from_exit / col_l_xy)
        fwhm_z = np.sqrt(rg_z**2 + float(options.iR)**2) / dz
        fwhm_xy = np.sqrt(rg_xy**2 + float(options.iR)**2) / dy
        sigma_z = fwhm_z / (2. * np.sqrt(2. * np.log(2.)))
        sigma_xy = np.sqrt(np.maximum(fwhm_xy**2 - 1., 1e-16)) / (2. * np.sqrt(2. * np.log(2.)))
    else:
        sigma_z = np.full(depth, float(sigma_z_input.reshape(-1)[0]), dtype=np.float32)
        sigma_xy = np.full(depth, float(sigma_xy_input.reshape(-1)[0]), dtype=np.float32)

    max_i = max(nx, ny, nz)
    coordinates = np.arange(max_i // 2 - 1, -max_i // 2, -1, dtype=np.float32)
    xx, yy = np.meshgrid(coordinates, coordinates, indexing='ij')
    kernel = np.exp(
        -(xx[:, :, None]**2 / (2. * sigma_z[None, None, :]**2)
          + yy[:, :, None]**2 / (2. * sigma_xy[None, None, :]**2))
    )
    mid_slice = kernel[:, :, kernel.shape[2] // 4]
    row, col = np.where(mid_slice > 1e-6)
    if row.size == 0 or col.size == 0:
        raise ValueError(f'Type-6 PSF for volume {volume} has empty support')
    kernel = kernel[row.min():row.max() + 1, col.min():col.max() + 1, :]
    return np.asfortranarray((kernel / np.sum(kernel, axis=(0, 1), keepdims=True)).astype(np.float32))


def _type6_multires_filter_resources(options: proj.projectorClass, volumes: int, depth_planes: np.ndarray) -> list[np.ndarray]:
    """Return one CDRF resource per volume for the Python custom operator."""
    provided = options.gFilter
    if isinstance(provided, (list, tuple)) and len(provided) > 0:
        if len(provided) != volumes:
            raise ValueError(f'gFilter must contain {volumes} filters, got {len(provided)}')
        filters = [np.asarray(value, dtype=np.float32) for value in provided]
    elif np.asarray(options.gFilter).size:
        raise ValueError(
            'A single gFilter is ambiguous for multi-resolution projector type 6. '
            'Provide one filter per volume as a gFilter list, or leave gFilter empty '
            'to generate the physical CDRFs automatically.'
        )
    else:
        filters = [
            _type6_auto_filter(options, volume, int(depth_planes[volume]))
            for volume in range(volumes)
        ]
    for volume, kernel in enumerate(filters):
        if kernel.ndim != 3 or min(kernel.shape) < 1:
            raise ValueError(f'gFilter[{volume}] must be a non-empty 3-D CDRF')
        filters[volume] = np.asfortranarray(kernel.astype(np.float32))
    return filters


def _type6_multires_geometry(options: proj.projectorClass, volumes: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create volume/view panel shifts without flattening image volumes into views."""
    panel_tilt = np.asarray(options.swivelAngles, dtype=np.float32).reshape(-1) - \
        np.asarray(options.angles, dtype=np.float32).reshape(-1) + 180.
    radius = np.asarray(options.radiusPerProj, dtype=np.float32).reshape(-1)
    if panel_tilt.size != radius.size:
        raise ValueError('Type-6 swivelAngles, angles, and radiusPerProj must have matching view counts')
    linear = radius * np.sin(np.deg2rad(panel_tilt))
    blur_planes = np.empty((volumes, radius.size), dtype=np.int32)
    blur_planes2 = np.empty_like(blur_planes)
    blur_planes2_linear = np.broadcast_to(linear, (volumes, radius.size)).copy()
    for volume in range(volumes):
        fov_x = _type6_scalar(options.FOVa_x, volume)
        dx = _type6_scalar(options.dx, volume)
        dy = _type6_scalar(options.dy, volume)
        if min(fov_x, dx, dy) <= 0.:
            raise ValueError(f'Type-6 volume {volume} has invalid FOV or voxel pitch')
        depth_mm = fov_x / 2. - (radius * np.cos(np.deg2rad(panel_tilt)) - float(options.CORtoDetectorSurface))
        blur_planes[volume] = np.rint(depth_mm / dx).astype(np.int32)
        # The custom helper shifts image axis 1 (the y axis), hence dy rather
        # than dx is the physical conversion for this in-plane panel offset.
        blur_planes2[volume] = np.rint(linear / dy).astype(np.int32)
    return blur_planes, blur_planes2, blur_planes2_linear


def _type6_total_lengths(options: proj.projectorClass) -> np.ndarray:
    """Full type-6 ray lengths, independent of EFOV volume decomposition."""
    full_fov_x = float(getattr(options, 'type6FullFOVaX', _type6_scalar(options.FOVa_x, 0)))
    full_nx = int(getattr(options, 'type6FullNx', int(np.asarray(options.Nx).reshape(-1)[0])))
    if full_fov_x <= 0. or full_nx < 1:
        raise ValueError('Type-6 full image geometry is invalid')
    full_dx = full_fov_x / full_nx
    panel_tilt = np.asarray(options.swivelAngles, dtype=np.float32).reshape(-1) - \
        np.asarray(options.angles, dtype=np.float32).reshape(-1) + 180.
    radius = np.asarray(options.radiusPerProj, dtype=np.float32).reshape(-1)
    depth_mm = full_fov_x / 2. - (radius * np.cos(np.deg2rad(panel_tilt)) - float(options.CORtoDetectorSurface))
    retained = np.maximum(1, full_nx - np.maximum(np.rint(depth_mm / full_dx).astype(np.int64), 0))
    return np.asarray(retained * full_dx, dtype=np.float32)


def SPECTParameters(options: proj.projectorClass):
    if options.projector_type in [1, 11, 12, 16, 2, 21, 22, 26, 61, 62]: # Ray tracing projectors
        nRays = int(options.n_rays_transaxial * options.n_rays_axial)
        if options.rayShiftsDetector.size == 0: # Collimator modeling
            options.rayShiftsDetector = np.zeros((2*nRays, options.nRowsD, options.nColsD, options.nHeads), dtype=np.float32)
            
            if options.colFxy == 0 and options.colFz == 0:
                dx = np.linspace(-(options.nRowsD / 2 - 0.5) * options.dPitchX, (options.nRowsD / 2 - 0.5) * options.dPitchX, options.nRowsD)
                dy = np.linspace(-(options.nColsD / 2 - 0.5) * options.dPitchY, (options.nColsD / 2 - 0.5) * options.dPitchY, options.nColsD)
                
                for ii in range(options.nRowsD):
                    for jj in range(options.nColsD):
                        for kk in range(nRays):
                            options.rayShiftsDetector[2 * kk, ii, jj, :] = -dx[ii]
                            options.rayShiftsDetector[2 * kk + 1, ii, jj, :] = -dy[jj]    

        if options.rayShiftsSource.size == 0:
            options.rayShiftsSource = np.zeros((2*nRays, options.nRowsD, options.nColsD, options.nHeads), dtype=np.float32)
            
            if nRays > 1: # Multiray shifts
                tmp_x, tmp_y = np.meshgrid(
                    np.linspace(-0.5, 0.5, options.n_rays_transaxial),
                    np.linspace(-0.5, 0.5, options.n_rays_axial)
                )
                if options.colFxy == 0 and options.colFz == 0: # Pinhole collimator
                    tmp_x *= options.dPitchX
                    tmp_y *= options.dPitchY
                elif np.isinf(options.colFxy) and np.isinf(options.colFz):  # Parallel-hole collimator
                    tmp_x *= 2 * options.colR
                    tmp_y *= 2 * options.colR

                tmp_shift = np.column_stack((tmp_x.ravel(), tmp_y.ravel())).T.reshape(-1, 1, order='F')

                for kk in range(nRays):
                    options.rayShiftsSource[2 * kk, :, :, :] = tmp_shift[2 * kk]
                    options.rayShiftsSource[2 * kk + 1, :, :, :] = tmp_shift[2 * kk + 1]

        if options.projector_type in [1, 11, 12, 16, 21, 61]:
            lengthXY = options.colD + 0.5 * options.colLxy
            lengthZ = options.colD + 0.5 * options.colLz
            if lengthXY != lengthZ:
                detectorShiftZ = options.rayShiftsDetector[1::2, :, :, :]
                options.rayShiftsSource[1::2, :, :, :] = detectorShiftZ + \
                    (options.rayShiftsSource[1::2, :, :, :] - detectorShiftZ) * (lengthXY / lengthZ)

        options.rayShiftsDetector = options.rayShiftsDetector.ravel('F')
        options.rayShiftsSource = options.rayShiftsSource.ravel('F')

    if options.projector_type in [2, 12, 21, 26, 22, 62]: # Orthogonal distance ray tracer
        if options.coneOfResponseStdCoeffA < 0:
            options.coneOfResponseStdCoeffA = 2*options.colR/options.colL
        if options.coneOfResponseStdCoeffB < 0:
            options.coneOfResponseStdCoeffB = 2*options.colR/options.colL*(options.colL+options.colD+options.cr_p/2)
        if options.coneOfResponseStdCoeffC < 0:
            options.coneOfResponseStdCoeffC = options.iR

    if options.projector_type in (6, 16, 26, 61, 62, 66): # Rotation-based projector side
        volume_count = max(
            np.asarray(options.Nx).size,
            np.asarray(options.Ny).size,
            np.asarray(options.Nz).size,
            np.asarray(options.dx).size,
            np.asarray(options.dy).size,
            np.asarray(options.dz).size,
        )
        if volume_count > 1:
            # Multi-resolution image geometry is established before this
            # function runs.  Keep it indexed by (volume, view): flattening
            # it here used to mix the volume and view axes, and failed before
            # custom type-6 operators could be constructed.
            if options.angles.size == 0:
                options.angles = (
                    np.repeat(options.startAngle, (options.nProjections // options.nHeads)) +
                    np.tile(
                        np.arange(0, options.angleIncrement * (options.nProjections / options.nHeads), options.angleIncrement),
                        (options.nHeads, 1),
                    )
                )
            options.angles = np.asarray(options.angles, dtype=np.float32).ravel('F')
            options.swivelAngles = np.asarray(options.swivelAngles, dtype=np.float32).ravel('F')
            options.radiusPerProj = np.asarray(options.radiusPerProj, dtype=np.float32).ravel('F')
            if not (options.angles.size == options.swivelAngles.size == options.radiusPerProj.size):
                raise ValueError('Multi-resolution type-6 geometry needs one angle, swivel angle, and radius per view')

            (
                options.blurPlanes,
                options.blurPlanes2,
                options.blurPlanes2Linear,
            ) = _type6_multires_geometry(options, volume_count)
            # Shifted CDRFs must include every depth plane selected by this
            # volume's panel geometry, not only its local Nx planes.
            nx = np.asarray(options.Nx, dtype=np.int64).reshape(-1)[:volume_count]
            filter_depth = nx + np.max(np.abs(options.blurPlanes), axis=1)
            options.gFilter = _type6_multires_filter_resources(options, volume_count, filter_depth)
            options.type6TotalLength = _type6_total_lengths(options)
            if not options.useTotLength:
                raise ValueError(
                    'Multi-resolution projector type 6 requires useTotLength=True so one '
                    'measurement is allocated across all volumes with one common ray length.'
                )

            # Type-6 PSF/geometry is always indexed by volume after setup.
            # Native callers extract entry zero at their compatibility boundary.
            options.uu = 1
            options.ub = 1
            return

        # Accept the same one-entry gFilter list used by the custom
        # multi-resolution path, while retaining the established local
        # single-volume calculation below.
        if isinstance(options.gFilter, (list, tuple)):
            if len(options.gFilter) > 1:
                raise ValueError('Single-volume projector type 6 accepts exactly one gFilter entry')
            options.gFilter = np.empty(0, dtype=np.float32) if len(options.gFilter) == 0 else np.asarray(options.gFilter[0])
        DistanceToFirstRow = 0.5 * options.dx
        # ``addProjector`` normalizes image dimensions to one-element arrays.
        # NumPy 2 no longer accepts such an array as the length for arange.
        nx = int(np.asarray(options.Nx).reshape(-1)[0])

        Distances = DistanceToFirstRow[..., np.newaxis] + np.arange(nx * 4, dtype=np.float32) * options.dx

        if options.gFilter.size == 0:
            if options.sigmaZ < 0.:
                col_l_xy = float(options.colLxy) if options.colLxy > 0. else float(options.colL)
                col_l_z = float(options.colLz) if options.colLz > 0. else float(options.colL)
                distance_from_exit = Distances + options.cr_p / 2.
                # Anger, "Scintillation Camera with Multichannel Collimators", J Nucl Med 5:515-531 (1964).
                # Axial and transaxial responses use their respective collimator septa lengths.
                rg_z = 2. * options.colR * distance_from_exit / col_l_z
                rg_xy = 2. * options.colR * distance_from_exit / col_l_xy
                rg_z[rg_z < 0] = 0.
                rg_xy[rg_xy < 0] = 0.
                FWHMrot = 1.

                fwhm_z = np.sqrt(rg_z**2 + options.iR**2)
                fwhm_xy = np.sqrt(rg_xy**2 + options.iR**2)
                fwhm_z_pixel = fwhm_z / options.dz[0]
                fwhm_xy_pixel = fwhm_xy / options.dy[0]
                expr = fwhm_xy_pixel**2 - FWHMrot**2
                expr[expr <= 0] = 10**-16

                # Rotation interpolation already contributes one pixel of in-plane blur; retain the existing compensation there.
                options.sigmaZ = fwhm_z_pixel / (2. * np.sqrt(2. * np.log(2.)))
                options.sigmaXY = np.sqrt(expr) / (2. * np.sqrt(2. * np.log(2.)))

            maxI = max(options.Nx[0].item(), max(options.Ny[0].item(), options.Nz[0].item()))
            y = np.arange(maxI // 2 - 1, -maxI // 2, -1, dtype=np.float32).reshape((1, -1), order='F')
            x = np.arange(maxI // 2 - 1, -maxI // 2, -1, dtype=np.float32).reshape((1, -1), order='F')
            xx = np.tile(x.T, (1, x.shape[1]))
            yy = np.tile(y, (xx.shape[1], 1))

            if np.any(options.sigmaXY < 0.):
                s1 = np.tile(options.sigmaZ**2, (xx.shape[0], yy.shape[1], 1))
                options.gFilter = (1 / (2 * np.pi * s1)) * np.exp(-(xx[:, :, None]**2 + yy[:, :, None]**2) / (2 * s1))
            else:
                s1 = np.tile(options.sigmaZ, (xx.shape[0], yy.shape[1], 1))
                s2 = np.tile(options.sigmaXY, (xx.shape[0], yy.shape[1], 1))
                options.gFilter = np.exp(-(xx[:, :, None]**2 / (2 * s1**2) + yy[:, :, None]**2 / (2 * s2**2)))


            mid_slice = options.gFilter[:, :, options.gFilter.shape[2] // 4]
            rowE, colE = np.where(mid_slice > 1e-6)
            rowS = rowE.min()
            colS = colE.min()
            rowE = rowE.max()
            colE = colE.max()

            options.gFilter = options.gFilter[rowS:rowE+1, colS:colE+1, :]
            options.gFilter /= np.sum(options.gFilter, axis=(0, 1), keepdims=True)


        panelTilt = options.swivelAngles - options.angles + 180
        options.blurPlanes = np.round((options.FOVa_x / 2 - (options.radiusPerProj * np.cos(np.deg2rad(panelTilt)) - options.CORtoDetectorSurface)) / options.dx)
        # Retain the detector-panel displacement in millimetres for the MPS custom operator.  It converts this distance to fractional pixels using the translated axis' pitch for the current image volume. Keep blurPlanes2 as the integer native-projector representation.
        options.blurPlanes2Linear = options.radiusPerProj * np.sin(np.deg2rad(panelTilt))
        options.blurPlanes2 = options.blurPlanes2Linear / options.dx

        if options.angles.size == 0:
            options.angles = (np.repeat(options.startAngle, (options.nProjections // options.nHeads)) + np.tile(np.arange(0,options.angleIncrement * (options.nProjections / options.nHeads),options.angleIncrement), (options.nHeads, 1)))
        options.uu = 1
        options.ub = 1

        options.gFilter = [np.asfortranarray(options.gFilter.astype(dtype=np.float32))]
        options.angles = options.angles.ravel('F').astype(dtype=np.float32)
        options.swivelAngles = options.swivelAngles.ravel('F').astype(dtype=np.float32)
        options.radiusPerProj = options.radiusPerProj.ravel('F').astype(dtype=np.float32)
        options.blurPlanes = options.blurPlanes.ravel('F').astype(dtype=np.int32)[None, :]
        options.blurPlanes2Linear = options.blurPlanes2Linear.ravel('F').astype(dtype=np.float32)[None, :]
        options.blurPlanes2 = options.blurPlanes2.ravel('F').astype(dtype=np.int32)[None, :]
        options.type6TotalLength = _type6_total_lengths(options)
