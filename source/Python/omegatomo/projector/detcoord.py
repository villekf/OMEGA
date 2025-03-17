# -*- coding: utf-8 -*-
import numpy as np
from omegatomo import proj
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
    if options.subsets > 1 and options.subsetType < 8:
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
                # x = np.zeros(options.Ndist * options.Nang * 4, 2);
                # y = np.zeros(options.Ndist * options.Nang * 4, 2);
                # for kk = 1 : options.nLayers
                #     if options.cryst_per_block(1) == options.cryst_per_block(2)
                #         [x1, y1] = sinogram_coordinates_2D(options, xp(1 + (kk - 1) * koko : kk * koko), yp(1 + (kk - 1) * koko : kk * koko));
                #     else
                #         if kk == 2
                #             [x1, y1] = sinogram_coordinates_2D(options, xp(1 + (kk - 1) * koko : kk * koko), yp(1 + (kk - 1) * koko : kk * koko), options.nLayers);
                #         else
                #             [x1, y1] = sinogram_coordinates_2D(options, xp(1 + (kk - 1) * koko : kk * koko), yp(1 + (kk - 1) * koko : kk * koko));
                #         end
                #     end
                #     x1(ismember(x1, [0 0], 'rows'),:) = repmat([inf inf], nnz(ismember(x1, [0 0], 'rows')), 1);
                #     y1(ismember(y1, [0 0], 'rows'),:) = repmat([inf inf], nnz(ismember(y1, [0 0], 'rows')), 1);
                #     x(1 + (kk - 1) * options.Ndist * options.Nang * 3: options.Ndist * options.Nang + (kk - 1) * options.Ndist * options.Nang * 3,:) = x1;
                #     y(1 + (kk - 1) * options.Ndist * options.Nang * 3: options.Ndist * options.Nang + (kk - 1) * options.Ndist * options.Nang * 3,:) = y1;
                # end
                # ind2 = 1 + options.Ndist * options.Nang * 3;
                # ind1 = options.Ndist * options.Nang;
                # x(1 + options.Ndist * options.Nang : options.Ndist * options.Nang * 2,:) = [x(ind2:end,1) x(1 : ind1, 2)];
                # y(1 + options.Ndist * options.Nang : options.Ndist * options.Nang * 2,:) = [y(ind2:end,1) y(1 : ind1, 2)];
                # x(1 + options.Ndist * options.Nang * 2 : options.Ndist * options.Nang * 3,:) = [x(1 : ind1, 1) x(ind2:end,2)];
                # y(1 + options.Ndist * options.Nang * 2 : options.Ndist * options.Nang * 3,:) = [y(1 : ind1, 1) y(ind2:end,2)];
            else:
                x, y = sinogramCoordinates2D(options, x, y)
    
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
    
def getCoordinatesSPECT(options: proj.projectorClass) -> Tuple[np.ndarray, np.ndarray]:
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
        alpha3: float = options.swivelAngles[ii]
        
        x[3, ii] = r1 * np.cos(np.deg2rad(alpha1)) + r2 * np.cos(np.deg2rad(alpha3))
        x[4, ii] = r1 * np.sin(np.deg2rad(alpha1)) + r2 * np.sin(np.deg2rad(alpha3))
        x[5, ii] = 0
        
        x[0, ii] = x[3, ii] + options.colL * np.cos(np.deg2rad(alpha3))
        x[1, ii] = x[4, ii] + options.colL * np.sin(np.deg2rad(alpha3))
        x[2, ii] = 0
        
        z[0, ii] = options.crXY * np.cos(np.deg2rad(alpha3 + 90))
        z[1, ii] = options.crXY * np.sin(np.deg2rad(alpha3 + 90))
    
    if options.flipImageX:
        x[0, :] = -x[0, :]
        x[3, :] = -x[3, :]
        z[0, :] = -z[0, :]
    
    if options.flipImageY:
        x[1, :] = -x[1, :]
        x[4, :] = -x[4, :]
        z[1, :] = -z[1, :]
    #return x, z
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
    
    # if nargout >= 6
    #     temppi = j*2 < -i;
    #     temppi2 = (i <= (j-det_w_pseudo/2)*2);
        
    #     temppi3 = false(det_w_pseudo);
    #     temppi3(tril(true(det_w_pseudo))) = temppi;
    #     temppi = logical(temppi3 + tril(temppi3,-1)');
        
    #     temppi3 = false(det_w_pseudo);
    #     temppi3(tril(true(det_w_pseudo))) = temppi2;
    #     temppi2 = logical(temppi3 + tril(temppi3,-1)');
        
    #     swap1 = triu(temppi);
    #     swap3 = tril(temppi);
    #     swap2 = triu(temppi2);
    #     swap4 = tril(temppi2);
    #     varargout{6} = cat(3, swap1, swap2, swap3, swap4);
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
    # Accumulate the coordinates
    # if mashing > 1:
    #     x = npg.aggregate(np.column_stack((i,j)), xx1, func='mean', size=[Ndist,Nang],fill_value=np.nan)
    #     y = npg.aggregate(np.column_stack((i,j)), yy1, func='mean', size=[Ndist,Nang],fill_value=np.nan)
    #     x2 = npg.aggregate(np.column_stack((i,j)), xx2, func='mean', size=[Ndist,Nang],fill_value=np.nan)
    #     y2 = npg.aggregate(np.column_stack((i,j)), yy2, func='mean', size=[Ndist,Nang],fill_value=np.nan)
    # else:
    np.add.at(x, (i, j), xx1)
    np.add.at(y, (i, j), yy1)
    np.add.at(x2, (i, j), xx2)
    np.add.at(y2, (i, j), yy2)

# Use np.add.at to accumulate values from xx1 at indices specified by i and j
    
    # Remove NaN values
    # Source: https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python
    # if np.sum(np.isnan(x.flatten())) > 0:
    #     xi = np.arange(0, x.shape[1])
    #     yi = np.arange(0, x.shape[0])
    #     #mask invalid values
    #     x = np.ma.masked_invalid(x)
    #     xx, yy = np.meshgrid(xi, yi)
    #     #get only the valid values
    #     x1 = xx[~x.mask]
    #     y1 = yy[~x.mask]
    #     newarr = x[~x.mask]
    #     x = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='linear')
    #     y = np.ma.masked_invalid(y)
    #     x1 = xx[~y.mask]
    #     y1 = yy[~y.mask]
    #     newarr = y[~y.mask]
    #     y = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='linear')
    #     x2 = np.ma.masked_invalid(x2)
    #     x1 = xx[~x2.mask]
    #     y1 = yy[~x2.mask]
    #     newarr = x2[~x2.mask]
    #     x2 = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='linear')
    #     y2 = np.ma.masked_invalid(y2)
    #     x1 = xx[~y2.mask]
    #     y1 = yy[~y2.mask]
    #     newarr = y2[~y2.mask]
    #     y2 = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='linear')
    
    # If mashing is present, combine the coordinates
    if mashing > 1:
        from skimage.measure import block_reduce
        # x = np.reshape(x,(Ndist,Nang));
        # x2 = np.reshape(x2,(Ndist,Nang));
        # y = np.reshape(y,(Ndist,Nang));
        # y2 = np.reshape(y2,(Ndist,Nang));
        # Compute the mean coordinates
        x = block_reduce(x, block_size=(1,mashing), func=np.mean)
        y = block_reduce(y, block_size=(1,mashing), func=np.mean)
        x2 = block_reduce(x2, block_size=(1,mashing), func=np.mean)
        y2 = block_reduce(y2, block_size=(1,mashing), func=np.mean)
        
        j = j // mashing
    
    x = np.column_stack((x.ravel('F'), x2.ravel('F')))
    y = np.column_stack((y.ravel('F'), y2.ravel('F')))
    return x, y
    
def sinogramCoordinates3D(options):
    import math
    cr_pz = options.cr_pz
    Nz = options.rings * 2 - 1
    
    if len(options.ringGaps) > 0 and np.sum(options.ringGaps) > 0:
        z_length = float(options.rings + 1 + len(options.ringGaps)) * cr_pz
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
        if len(options.ringGaps) > 0 and np.sum(options.ringGaps) > 0:
            z = np.linspace(cr_pz, z_length + cr_pz, options.rings + 2 + len(options.ringGaps))
        else:
            z = np.linspace(cr_pz, z_length + cr_pz, options.rings + 2)
        
        if len(options.ringGaps) > 0 and np.sum(options.ringGaps) > 0:
            z = z[:options.rings + len(options.ringGaps)]
            z = np.delete(z, options.ringGaps + 1)
        else:
            z = z[:options.rings]
        
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
        if len(options.ringGaps) > 0 and sum(options.ringGaps) > 0:
            z = np.zeros(((options.rings + len(options.ringGaps))**2, 2))
            z2 = np.zeros(((options.rings + len(options.ringGaps))**2, 2))
            loppu = options.rings + len(options.ringGaps)
        else:
            z = np.zeros((options.rings**2, 2))
            z2 = np.zeros((options.rings**2, 2))
            loppu = options.rings
        
        ind = np.arange(0, cr_pz * loppu, cr_pz)
        
        if options.nLayers > 1 and options.cryst_per_block_axial[0] != options.cryst_per_block_axial[1]:
            ind += apu
            layers = np.repeat(np.array([[0], [1]]), options.cryst_per_block_axial[1], axis=1)
            layers = np.repeat(np.concatenate((layers, np.array([[0]])), axis=1), options.linear_multip, axis=1)
        
        for t in range(1, loppu + 1):
            if options.nLayers > 1 and options.cryst_per_block_axial[0] != options.cryst_per_block_axial[1]:
                z[(t - 1) * loppu:t * loppu, :] = np.column_stack((dif * (t - 1) * np.ones(loppu), ind))
                z2[(t - 1) * loppu:t * loppu, :] = np.column_stack((layers[t - 1] * np.ones(loppu), layers))
            else:
                z[(t - 1) * loppu:t * loppu, :] = np.column_stack((dif * (t - 1) * np.ones(loppu), ind))
        
        if len(options.ringGaps) > 0 and sum(options.ringGaps) > 0:
            ind = np.tile(options.ringGaps, (loppu, 1)) + np.repeat(np.arange(loppu), 3) * loppu + 1
            ind += np.tile(np.arange(len(options.ringGaps)), (len(ind) // len(options.ringGaps), 1))
            z[ind.flatten() - 1, :] = []
            for kk in range(len(options.ringGaps), 0, -1):
                z[options.rings * (options.ringGaps[kk - 1] + kk - 2):options.rings * (options.ringGaps[kk - 1] + kk - 1), :] = []
        z -= (maxZ / 2 - cr_pz)
        
        if options.nLayers > 1 and options.cryst_per_block_axial[0] != options.cryst_per_block_axial[1]:
            apuZ = np.zeros(options.rings**2)
            apuZ[np.all(z2 == [1, 1], axis=1)] = 3
            apuZ[np.all(z2 == [1, 0], axis=1)] = 1
            apuZ[np.all(z2 == [0, 1], axis=1)] = 2
            z = np.column_stack((apuZ, z))
    return z


def SPECTParameters(options):
    import math
    DistanceToFirstRow = options.radiusPerProj - (options.Nx.item() / 2. - 0.5) * options.dx[0].item()
    DistanceToFirstRow = DistanceToFirstRow.reshape((-1, 1), order='F')
    Distances = np.asfortranarray(np.tile(DistanceToFirstRow,(1,options.Nx.item())) + np.tile(np.arange(0, options.Nx.item(),dtype=np.float32) * options.dx[0].item(), (DistanceToFirstRow.size,1)))
    Distances = Distances - options.colL - options.colD #these are distances to the actual detector surface
    
    if options.gFilter.size == 0:
        if options.sigmaZ < 0.:
            Rg = 2. * options.colR * (options.colL + options.colD + Distances + options.cr_p / 2.) / options.colL #Anger, "Scintillation Camera with Multichannel Collimators", J Nucl Med 5:515-531 (1964)
            Rg[Rg < 0] = 0.
            FWHMrot = 1.
    
            FWHM = np.sqrt(Rg**2 + options.iR**2)
            FWHM_pixel = FWHM / options.dx[0]
            expr = FWHM_pixel**2 - FWHMrot**2
            expr[expr <= 0] = 10**-16;
            FWHM_WithinPlane = np.sqrt(expr)
    
            #Parametrit CDR-mallinnukseen
            options.sigmaZ = FWHM_pixel / (2. * math.sqrt(2. * math.log(2.)))
            options.sigmaXY = FWHM_WithinPlane / (2. * math.sqrt(2. * math.log(2.)))
        maxI = max(options.Nx[0].item(), max(options.Ny[0].item(), options.Nz[0].item()))
        y = np.arange(maxI // 2 - 1, -maxI // 2, -1, dtype=np.float32).reshape((1, -1), order='F')
        x = np.arange(maxI // 2 - 1, -maxI // 2, -1, dtype=np.float32).reshape((1, -1), order='F')
        xx = np.tile(x.T, (1, x.shape[1]))
        yy = np.tile(y, (xx.shape[1], 1))
        if np.any(options.sigmaXY < 0.):
            options.sigmaZ = np.reshape(options.sigmaZ, (options.sigmaZ.shape[0], options.sigmaZ.shape[1], 1, 1), order='F')
            s1 = np.asfortranarray(np.tile(np.transpose(options.sigmaZ, (3, 2, 1, 0)), (xx.shape[0], yy.shape[1], 1, 1)))
            xx = np.asfortranarray(np.reshape(xx,(xx.shape[0], xx.shape[1], 1, 1), order='F'))
            yy = np.asfortranarray(np.reshape(yy,(yy.shape[0], yy.shape[1], 1, 1), order='F'))
            options.gFilter = (1. / (2. * np.pi * s1) * np.exp(-(xx**2 + yy**2) / (2. * s1)))
        else:
            options.sigmaZ = np.reshape(options.sigmaZ, (options.sigmaZ.shape[0], options.sigmaZ.shape[1], 1, 1), order='F')
            options.sigmaXY = np.reshape(options.sigmaXY, (options.sigmaXY.shape[0], options.sigmaXY.shape[1], 1, 1), order='F')
            s1 = np.asfortranarray(np.tile(np.transpose(options.sigmaZ, (3, 2, 1, 0)), (xx.shape[0], yy.shape[1], 1, 1)))
            s2 = np.asfortranarray(np.tile(np.transpose(options.sigmaXY, (3, 2, 1, 0)), (xx.shape[0], yy.shape[1], 1, 1)))
            xx = np.asfortranarray(np.reshape(xx,(xx.shape[0], xx.shape[1], 1, 1), order='F'))
            yy = np.asfortranarray(np.reshape(yy,(yy.shape[0], yy.shape[1], 1, 1), order='F'))
            options.gFilter = np.exp(-(xx**2 / (2. * s1**2) + yy**2 / (2. * s2**2)))
        ind = np.argmax((options.radiusPerProj))
        
        rowE = np.where(options.gFilter[:, :, -1, ind] > 1e-6)[0].max()
        colE = np.where(options.gFilter[:, :, -1, ind] > 1e-6)[-1].max()
        
        rowS = np.where(options.gFilter[:, :, -1, ind] > 1e-6)[0].min()
        colS = np.where(options.gFilter[:, :, -1, ind] > 1e-6)[1].min()
        # rowE,colE = np.nonzero(options.gFilter[:,:,-1,ind] > 1e-6)
        # rowE = rowE[-1]
        # colE = colE[-1]
        # [rowS,colS] = np.nonzero(options.gFilter[:,:,-1,ind] > 1e-6,1)
        # rowS = rowE[0]
        # colS = colE[0]
        options.gFilter = options.gFilter[rowS:rowE+1,colS:colE+1,:,:]
        options.gFilter = options.gFilter / np.sum(options.gFilter)
        options.gFilter = options.gFilter.astype(dtype=np.float32)
    options.blurPlanes = np.argmax(Distances>0,1)
    if options.angles.size == 0:
        options.angles = (np.repeat(options.startAngle, (options.nProjections // options.nHeads)) + np.tile(np.arange(0,options.angleIncrement * (options.nProjections / options.nHeads),options.angleIncrement), (options.nHeads, 1)))
    options.uu = 1
    options.ub = 1
    if abs(options.offangle) > 0:
        if np.max(np.abs(options.angles.flatten())) > 10. * np.pi:
            options.angles = options.angles + (options.offangle * 180./np.pi)
        else:
            options.angles = options.angles + options.offangle
    if np.max(np.abs(options.angles.flatten())) > 10. * np.pi and options.implementation == 2:
        options.angles = options.angles / 180. * np.pi
    if options.flip_image:
        options.angles = -(options.angles)
    options.angles = options.angles.ravel('F').astype(dtype=np.float32)
    options.blurPlanes = options.blurPlanes.ravel('F').astype(dtype=np.uint32)
