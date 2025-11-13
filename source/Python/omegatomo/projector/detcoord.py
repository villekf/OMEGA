# -*- coding: utf-8 -*-
import numpy as np
from omegatomo.projector import proj
from typing import Tuple

def CTDetSource(options):
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
    if options.x.size != options.nColsD * options.nRowsD * options.nProjections and options.uV.size <= 1:
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
    else:
        options.z[:,0] = options.z[:,0] * options.dPitchX
        options.z[:,1] = options.z[:,1] * options.dPitchX
        
def CTDetectorCoordinates(angles, pitchRoll = np.empty(0, dtype=np.float32)):
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
    if isinstance(options.x, np.ndarray) and options.x.size > 1 and options.y.size > 1 and options.listmode == 0:
        if not(options.z.size == options.x.size) and options.z.size < options.nProjections:
            options.z = np.zeros(options.x.shape,dtype=np.float32)
        x = options.x
        y = options.y
        z = options.z
    else:
        if options.use_raw_data == 0:
            x,y = detectorCoordinates(options);
            if options.nLayers > 1:
                koko = x.size / 2
            else:
                x, y = sinogramCoordinates2D(options, x, y)
                
            if options.arc_correction:
                from omegatomo.util.arcCorrection import arc_correction
                x, y, options = arc_correction(options, False);
            if options.sampling > 1:
                from omegatomo.util.sampling import increaseSampling
                x, y, options = increaseSampling(options, x, y, False);
    
            # if options.arc_correction and ~options.precompute_lor
            #     [x, y, options] = arcCorrection(options, xp, yp, interpolateSinogram);
            # if options.sampling > 1 and ~options.precompute_lor
            #     [x, y, options] = increaseSampling(options, x, y, interpolateSinogram);
            z = sinogramCoordinates3D(options)
            if not(options.NSinos == options.TotSinos):
                z = z[0:options.NSinos,:]
            x = np.column_stack((x[:,0], y[:,0], x[:,1], y[:,1]))
        else:
            if options.det_per_ring < options.det_w_pseudo:
                options.offangle = options.offangle / options.det_w_pseudo;
                options.offangle = options.offangle * options.det_per_ring;
            x, y = detectorCoordinates(options);
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

        x[3, ii] = x[0, ii] + (options.colD + 0.5 * options.colL) * np.cos(np.deg2rad(alpha2))
        x[4, ii] = x[1, ii] + (options.colD + 0.5 * options.colL) * np.sin(np.deg2rad(alpha2))
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
    cr_p = options.cr_p
    diameter = options.diameter
    if isinstance(options.cryst_per_block, np.ndarray):
        cryst_per_block = options.cryst_per_block[0].item()
    else:
        cryst_per_block = options.cryst_per_block
    blocks_per_ring = options.blocks_per_ring
    x = 0
    y = 0
    
    DOI = options.DOI;
    
    transaxial_multip = options.transaxial_multip;
    
    diameter = diameter + DOI * 2;
    
    
    if options.det_w_pseudo > options.det_per_ring:
        
        xp,yp = computeCoordinates(options, blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, True)
        x = xp
        y = yp
        
    elif options.nLayers > 1:
        if options.cryst_per_block[0] == options.cryst_per_block[1]:
            x1,y1 = computeCoordinates(options, blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, False);
            # orig_diameter = diameter + options.crystH(end) * 2;
            diameter = diameter + DOI * 2 + options.crystH[0] * 2
            x2,y2 = computeCoordinates(options, blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, False);
            # xp = [x;xp];
            # yp = [y;yp];
        elif options.cryst_per_block[0] > options.cryst_per_block[1]:
            x,y = computeCoordinates(options, blocks_per_ring, transaxial_multip,options.cryst_per_block[-1].item(), diameter, cr_p, True);
            # orig_diameter = diameter + options.crystH(end) * 2;
            # diameter = diameter + DOI * 2 + options.crystH(end) * 2;
            # [x,y] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, false);
            # xp = [x;xp];
            # yp = [y;yp];
        else:
            x,y = computeCoordinates(options, blocks_per_ring, transaxial_multip,options.cryst_per_block[-1].item(), diameter, cr_p, False);
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
    if not(options.offangle == 0):
        x = np.roll(x, int(np.round(options.offangle)))
        y = np.roll(y, int(np.round(options.offangle)))
    return x, y

def computeCoordinates(options, blocks_per_ring, transaxial_multip, cryst_per_block, diameter, cr_p, usePseudo):
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
    
    alkupistex = (diameter / 2)
    alkupistey = -((cryst_per_block_orig) / 2 + 0.5) * cr_p
    
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
    det_w_pseudo = options.det_w_pseudo
    Nang = options.Nang
    Ndist = options.Ndist
    mashing = 1
    # Determine the possible mashing factor
    if Nang < options.det_w_pseudo//2:
        mashing = (options.det_w_pseudo // Nang // 2)
        Nang = Nang * mashing
    
    
    ## 2D coordinates
    
    # Determine the sinogram indices for each of the LOR
    
    # Form the detector vector pair
    if isinstance(options.cryst_per_block, np.ndarray):
        L = formDetectorIndices(det_w_pseudo, nLayers, options.cryst_per_block[0].item());
    else:
        L = formDetectorIndices(det_w_pseudo, nLayers, options.cryst_per_block);
    
    L = L - 1
    
    xa = np.max(L,1)
    ya = np.min(L,1)
    
    # Angle
    j = ((xa + ya + det_w_pseudo //2 + 1) % det_w_pseudo) // 2
    
    b = j + det_w_pseudo // 2;
    
    # Distance
    i = np.abs(xa - ya - det_w_pseudo // 2);
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


def SPECTParameters(options: proj.projectorClass) -> None:
    if options.projector_type in [1, 11, 12, 2, 21, 22]: # Ray tracing projectors
        if options.rayShiftsDetector.size == 0: # Collimator modeling
            options.rayShiftsDetector = np.zeros((2*options.nRays, options.nRowsD, options.nColsD, options.nProjections), dtype=np.float32)
            
            if options.colFxy == 0 and options.colFz == 0:
                dx = np.linspace(-(options.nRowsD / 2 - 0.5) * options.dPitchX, (options.nRowsD / 2 - 0.5) * options.dPitchX, options.nRowsD)
                dy = np.linspace(-(options.nColsD / 2 - 0.5) * options.dPitchY, (options.nColsD / 2 - 0.5) * options.dPitchY, options.nColsD)
                
                for ii in range(options.nRowsD):
                    for jj in range(options.nColsD):
                        for kk in range(options.nRays):
                            options.rayShiftsDetector[2 * kk, ii, jj, :] = -dx[ii]
                            options.rayShiftsDetector[2 * kk + 1, ii, jj, :] = -dy[jj]    

        if options.rayShiftsSource.size == 0:
            options.rayShiftsSource = np.zeros((2*options.nRays, options.nRowsD, options.nColsD, options.nProjections), dtype=np.float32)
            
            if options.nRays > 1: # Multiray shifts
                nRays = int(np.sqrt(options.nRays))
                tmp_x, tmp_y = np.meshgrid(np.linspace(-0.5, 0.5, nRays), np.linspace(-0.5, 0.5, nRays))
                if options.colFxy == 0 and options.colFz == 0: # Pinhole collimator
                    tmp_x *= options.dPitchX
                    tmp_y *= options.dPitchY
                elif np.isinf(options.colFxy) and np.isinf(options.colFz):  # Parallel-hole collimator
                    tmp_x *= 2 * options.colR
                    tmp_y *= 2 * options.colR

                tmp_shift = np.column_stack((tmp_x.ravel(), tmp_y.ravel())).T.reshape(-1, 1, order='F')

                for kk in range(options.nRays):
                    options.rayShiftsSource[2 * kk, :, :, :] = tmp_shift[2 * kk]
                    options.rayShiftsSource[2 * kk + 1, :, :, :] = tmp_shift[2 * kk + 1]
            
        options.rayShiftsSource = options.rayShiftsSource.ravel('F')

    if options.projector_type in [12, 21, 2, 22]: # Orthogonal distance ray tracer
        if options.coneOfResponseStdCoeffA < 0:
            options.coneOfResponseStdCoeffA = 2*options.colR/options.colL
        if options.coneOfResponseStdCoeffB < 0:
            options.coneOfResponseStdCoeffB = 2*options.colR/options.colL*(options.colL+options.colD+options.cr_p/2)
        if options.coneOfResponseStdCoeffC < 0:
            options.coneOfResponseStdCoeffC = options.iR

    if options.projector_type == 6: # Rotation-based projector
        DistanceToFirstRow = 0.5 * options.dx
        Distances = DistanceToFirstRow[..., np.newaxis] + np.arange(options.Nx * 4, dtype=np.float32) * options.dx
        Distances -= (options.colL + options.colD)  # distances to detector surface

        if options.gFilter.size == 0:
            if options.sigmaZ < 0.:
                Rg = 2. * options.colR * (options.colL + options.colD + Distances + options.cr_p / 2.) / options.colL #Anger, "Scintillation Camera with Multichannel Collimators", J Nucl Med 5:515-531 (1964)
                Rg[Rg < 0] = 0.
                FWHMrot = 1.

                FWHM = np.sqrt(Rg**2 + options.iR**2)
                FWHM_pixel = FWHM / options.dx[0]
                expr = FWHM_pixel**2 - FWHMrot**2
                expr[expr <= 0] = 10**-16

                FWHM_WithinPlane = np.sqrt(expr)
                #Parametrit CDR-mallinnukseen
                options.sigmaZ = FWHM_pixel / (2. * np.sqrt(2. * np.log(2.)))
                options.sigmaXY = FWHM_WithinPlane / (2. * np.sqrt(2. * np.log(2.)))

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
        options.blurPlanes2 = options.radiusPerProj * np.sin(np.deg2rad(panelTilt)) / options.dx

        if options.angles.size == 0:
            options.angles = (np.repeat(options.startAngle, (options.nProjections // options.nHeads)) + np.tile(np.arange(0,options.angleIncrement * (options.nProjections / options.nHeads),options.angleIncrement), (options.nHeads, 1)))
        options.uu = 1
        options.ub = 1

        options.gFilter = np.asfortranarray(options.gFilter.astype(dtype=np.float32))
        options.angles = options.angles.ravel('F').astype(dtype=np.float32)
        options.swivelAngles = options.swivelAngles.ravel('F').astype(dtype=np.float32)
        options.radiusPerProj = options.radiusPerProj.ravel('F').astype(dtype=np.float32)
        options.blurPlanes = options.blurPlanes.ravel('F').astype(dtype=np.int32)
        options.blurPlanes2 = options.blurPlanes2.ravel('F').astype(dtype=np.int32)

        return None