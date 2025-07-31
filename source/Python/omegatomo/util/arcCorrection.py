# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 11:40:44 2025

@author: Ville-Veikko Wettenhovi
"""

import numpy as np

def arc_correction(options, interpolate_sinogram):
    from omegatomo.projector.detcoord import detectorCoordinates, sinogramCoordinates2D
    xp, yp = detectorCoordinates(options)
    x_o, y_o = sinogramCoordinates2D(options, xp, yp)

    new_xp = np.zeros_like(xp)
    new_yp = np.zeros_like(yp)

    # Detector angles
    if options.blocks_per_ring % 4 == 0:
        l_angles = np.linspace(0, 90, options.blocks_per_ring // 4 + 1)
    else:
        l_angles = np.linspace(0, 180, options.blocks_per_ring // 2 + 1)
        l_angles = l_angles[:int(np.ceil(options.blocks_per_ring / 4))]
    ll = 0  # Python index

    # Shift to zero angle
    shift_val = int(options.offangle + (options.cryst_per_block + 1) / 2 if options.det_w_pseudo > options.det_per_ring else options.offangle + options.cryst_per_block / 2)
    xp = np.roll(xp, -shift_val)
    yp = np.roll(yp, -shift_val)

    # xp -= options.diameter / 2
    # yp -= options.diameter / 2

    quarter_len = len(xp) // 4
    for kk in range(quarter_len):
        if options.det_w_pseudo > options.det_per_ring:
            r = int(options.det_w_pseudo / 2)
            val_range = np.arange(kk + r - (options.cryst_per_block + 1), kk + r + (options.cryst_per_block + 1) + 1)
        else:
            r = int(options.det_w_pseudo / 2)
            val_range = np.arange(kk + r - options.cryst_per_block, kk + r + options.cryst_per_block + 1)

        angle_ref = l_angles[ll]
        dx = xp[kk] - xp[val_range]
        dy = yp[kk] - yp[val_range]
        angle1 = np.round(np.degrees(np.arctan2(dy, dx)) * 1e4) / 1e4

        if options.blocks_per_ring % 4 == 0:
            diffs = np.abs(np.abs(angle1) - angle_ref)
            ind1 = np.argmin(diffs)
        else:
            ind1 = np.argmin(np.abs(np.abs(angle1) - angle_ref))

        x2 = xp[val_range[ind1]]
        y2 = yp[val_range[ind1]]

        p = np.array([xp[kk], yp[kk]])
        q = np.array([x2, y2])
        d = q - p
        d_dot = np.dot(d, d)
        p_dot = np.dot(p, p)
        pd_dot = np.dot(2 * p, d)
        rad2 = (options.diameter / 2) ** 2

        sqrt_term = pd_dot ** 2 - 4 * d_dot * (p_dot - rad2)
        l = (-pd_dot - np.sqrt(sqrt_term)) / (2 * d_dot)

        lx = p + l * d
        new_xp[kk] = lx[0]
        new_yp[kk] = lx[1]

    # Mirror fill
    new_xp[:quarter_len] += options.diameter / 2
    new_yp[:quarter_len] += options.diameter / 2

    new_yp[quarter_len:2*quarter_len] = np.flip(new_yp[:quarter_len])
    diffi = np.diff(np.concatenate([np.flip(new_yp[:quarter_len]), [options.diameter/2]]))
    diffi = options.diameter/2 + np.cumsum(np.flip(diffi))
    new_yp[2*quarter_len:] = np.concatenate([diffi, np.flip(diffi)])

    diffi = np.diff(np.concatenate([new_xp[:quarter_len], [options.diameter/2]]))
    diffi = options.diameter/2 + np.cumsum(np.flip(diffi))
    new_xp[quarter_len:2*quarter_len] = diffi
    new_xp[2*quarter_len:] = np.flip(new_xp[:2*quarter_len])

    # Undo shift
    xp = np.roll(new_xp, shift_val)
    yp = np.roll(new_yp, shift_val)

    x, y = sinogramCoordinates2D(options, xp, yp)
    
    # Reshape LOR coordinates
    xx1 = x[:, 0].reshape(options.Ndist, options.Nang, order='F')
    xx2 = x[:, 1].reshape(options.Ndist, options.Nang, order='F')
    yy1 = y[:, 0].reshape(options.Ndist, options.Nang, order='F')
    yy2 = y[:, 1].reshape(options.Ndist, options.Nang, order='F')
    
    # Compute LOR angles (degrees), shift to avoid negatives, find J
    angle = np.degrees(np.arctan((yy1 - yy2) / (xx1 - xx2))) + 90
    angle[angle == 180] = 0
    J = np.argmin(np.mean(angle, axis=0))
    
    # Create arc-corrected coordinates from perpendicular LOR
    eka = xx1[0, J]
    vika = xx1[-1, J]
    alkux = np.linspace(eka, vika, options.Ndist)
    
    alkuy = np.sqrt((options.diameter / 2)**2 - (alkux - options.diameter / 2)**2) + options.diameter / 2
    alku = np.asfortranarray(np.vstack([alkux, alkuy]))
    alku2 = alku - options.diameter / 2
    
    alku = np.asfortranarray(np.vstack([alkux, np.abs(alkuy - options.diameter)]))
    alku1 = alku - options.diameter / 2
    
    # Compute angles
    angles = np.linspace(0, 180, options.Nang + 1)[:-1]
    angles = np.roll(angles, J + 1)
    angles = angles.reshape(1, 1, -1)
    
    # Rotation matrices
    cos_a = np.cos(np.radians(angles)).squeeze()
    sin_a = np.sin(np.radians(angles)).squeeze()
    rot_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]], order='F')  # shape: (2, 2, Nang)
    
    # Rotate all points with each matrix (broadcasting)
    new_xy1 = np.zeros((options.Ndist * options.Nang,2),order='F')
    new_xy2 = np.zeros((options.Ndist * options.Nang,2),order='F')
    for ii in range(0, rot_matrix.shape[2]):
        new_xy1[ii * options.Ndist:(ii + 1) * options.Ndist] = (rot_matrix[:,:,ii] @ alku1).T + options.diameter / 2
        new_xy2[ii * options.Ndist:(ii + 1) * options.Ndist] = (rot_matrix[:,:,ii] @ alku2).T + options.diameter / 2
    # new_xy1 = np.einsum('ijk,k->ji', rot_matrix, alku1) + options.diameter / 2
    # new_xy2 = np.einsum('ijk,k->ji', rot_matrix, alku2) + options.diameter / 2

    # Multiply each matrix with alku1 and alku2 (2 x Ndist vectors)
    # Result will be shape (options.Nang, 2) after transpose
    # new_xy1 = np.array([(R @ alku1).T for R in rot_matrix]) + options.diameter / 2
    # new_xy2 = np.array([(R @ alku2).T for R in rot_matrix]) + options.diameter / 2
    # Update x and y arrays
    x = np.column_stack([new_xy1[:, 0], new_xy2[:, 0]])
    y = np.column_stack([new_xy1[:, 1], new_xy2[:, 1]])
    
    # Flip half of the LORs to match symmetry (MATLAB-style)
    xx1 = x[:, 0].reshape(options.Ndist, options.Nang, order='F')
    xx2 = x[:, 1].reshape(options.Ndist, options.Nang, order='F')
    yy1 = y[:, 0].reshape(options.Ndist, options.Nang, order='F')
    yy2 = y[:, 1].reshape(options.Ndist, options.Nang, order='F')
    
    apu = xx1[:, :J].copy()
    xx1[:, :J] = np.flipud(xx2[:, :J])
    xx2[:, :J] = np.flipud(apu)
    
    apu = yy1[:, :J].copy()
    yy1[:, :J] = np.flipud(yy2[:, :J])
    yy2[:, :J] = np.flipud(apu)
    
    # Flatten back to 2D coordinates
    x[:, 0] = xx1.ravel('F')
    x[:, 1] = xx2.ravel('F')
    y[:, 0] = yy1.ravel('F')
    y[:, 1] = yy2.ravel('F')

    x -= options.diameter / 2
    y -= options.diameter / 2
    
    # Reshape original coordinates too
    xx1_o = x_o[:, 0].reshape(options.Ndist, options.Nang, order='F')
    xx2_o = x_o[:, 1].reshape(options.Ndist, options.Nang, order='F')
    yy1_o = y_o[:, 0].reshape(options.Ndist, options.Nang, order='F')
    yy2_o = y_o[:, 1].reshape(options.Ndist, options.Nang, order='F')
    
    if interpolate_sinogram:

        from scipy.interpolate import griddata
        if options.SinM.size == options.SinM.shape[0] or options.SinM.ndim < 4:
            options.SinM = options.SinM.reshape((options.Ndist, options.Nang, options.NSinos, -1), order='F')
    
        angles_o = np.degrees(np.arctan((yy1_o - yy2_o) / (xx1_o - xx2_o)))
        angles_o -= np.min(angles_o[:, 0])
        angles_o[angles_o < 0] += 180
        angle_o = np.concatenate([
            -angles_o[:, [1]],
            angles_o,
            np.abs(angles_o[:, [0]] - 180)
        ], axis=1)
    
        angles = np.degrees(np.arctan((yy1 - yy2) / (xx1 - xx2)))
        angles -= np.min(angles[:, 0])
        angles[angles < 0] += 180
        angle = np.concatenate([
            -angles[:, [1]],
            angles,
            np.abs(angles[:, [0]] - 180)
        ], axis=1)
    
        # x0 = options.diameter / 2
        # y0 = options.diameter / 2
        x0 = 0.
        y0 = 0.
    
        def orth_dist(x1, y1, x2, y2):
            return np.abs((y1 - y2) * x0 - (x1 - x2) * y0 + x1 * y2 - y1 * x2) / \
                   np.sqrt((y1 - y2)**2 + (x1 - x2)**2)
    
        distance_o = orth_dist(x_o[:, 0], y_o[:, 0], x_o[:, 1], y_o[:, 1]).reshape(options.Ndist, options.Nang,order='F')
        distance_o[options.Ndist // 2:] *= -1
        distance_o = np.concatenate([
            distance_o[:, [0]],
            distance_o,
            distance_o[:, [0]]
        ], axis=1)
    
        distance = orth_dist(x[:, 0], y[:, 0], x[:, 1], y[:, 1]).reshape(options.Ndist, options.Nang, order='F')
        distance[options.Ndist // 2:] *= -1
        distance = np.concatenate([
            distance[:, [0]],
            distance,
            distance[:, [0]]
        ], axis=1)
    
        for uu in range(options.SinM.shape[3]):
            sinm = options.SinM[..., uu]
            sinm_ext = np.concatenate([
                sinm[:, -1:, :],
                sinm,
                sinm[:, :1, :]
            ], axis=1)
            
            for zz in range(sinm_ext.shape[2]):
                values = sinm_ext[:, :, zz].astype(np.float32)
            
                # Flatten the grid and values
                points = np.column_stack((angle_o.ravel('F'), distance_o.ravel('F')))
                vals = values.ravel('F')
            
                # Interpolate
                interpolated = griddata(
                    points,
                    vals,
                    (angle, distance),
                    method=options.arc_interpolation,
                    fill_value=0.0
                )
            
                interpolated[np.isnan(interpolated)] = 0.0
                options.SinM[:, :, zz, uu] = interpolated[:,1:interpolated.shape[1]-1].astype(np.float32)
    
        if options.verbose:
            print("Arc correction interpolation complete")

    return x, y, options