# -*- coding: utf-8 -*-
"""
Copyright (C) 2024-2025 Ville-Veikko Wettenhovi

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


import numpy as np
    
def computePixelSize(options):
    """
    This function computes the location of the FOV as well as the voxel sizes.
    Separate values are computed for multi-resolution volumes.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    xx : NumPy array
        The coordinates of the voxel boundaries in x-direction.
    yy : NumPy array
        The coordinates of the voxel boundaries in y-direction.
    zz : NumPy array
        The coordinates of the voxel boundaries in z-direction.

    """
    FOV = np.column_stack((options.FOVa_x,options.FOVa_y, options.axial_fov)).astype(dtype=np.float32)
    etaisyys = -(FOV) / 2
    options.dx = np.zeros((FOV.shape[0]),dtype=np.float32)
    options.dy = np.zeros((FOV.shape[0]),dtype=np.float32)
    options.dz = np.zeros((FOV.shape[0]),dtype=np.float32)
    options.bx = np.zeros((FOV.shape[0]),dtype=np.float32)
    options.by = np.zeros((FOV.shape[0]),dtype=np.float32)
    options.bz = np.zeros((FOV.shape[0]),dtype=np.float32)
    for kk in range(FOV.shape[0] - 1, -1, -1):
        if isinstance(options.Nx, np.ndarray):
            xx = np.linspace(etaisyys[kk,0] + options.oOffsetX, -etaisyys[kk,0] + options.oOffsetX, options.Nx[kk].item() + 1, dtype=np.float32)
            yy = np.linspace(etaisyys[kk,1] + options.oOffsetY, -etaisyys[kk,1] + options.oOffsetY, options.Ny[kk].item() + 1, dtype=np.float32)
            zz = np.linspace(etaisyys[kk,2] + options.oOffsetZ, -etaisyys[kk,2] + options.oOffsetZ, options.Nz[kk].item() + 1, dtype=np.float32)
        else:
            xx = np.linspace(etaisyys[kk,0] + options.oOffsetX, -etaisyys[kk,0] + options.oOffsetX, options.Nx + 1, dtype=np.float32)
            yy = np.linspace(etaisyys[kk,1] + options.oOffsetY, -etaisyys[kk,1] + options.oOffsetY, options.Ny + 1, dtype=np.float32)
            zz = np.linspace(etaisyys[kk,2] + options.oOffsetZ, -etaisyys[kk,2] + options.oOffsetZ, options.Nz + 1, dtype=np.float32)
        # Distance of adjacent pixels
        options.dx[kk] = xx[1] - xx[0]
        options.dy[kk] = yy[1] - yy[0]
        options.dz[kk] = zz[1] - zz[0]
    
        # Distance of image from the origin
        if kk == 0:
            options.bx[kk] = xx[0]
            options.by[kk] = yy[0]
            options.bz[kk] = zz[0]
        else:
            if kk > 4 or (FOV.shape[0] == 5 and kk > 2):
                if kk % 2 == 0:
                    options.by[kk] = options.oOffsetY + FOV[0][1] / 2
                else:
                    options.by[kk] = options.oOffsetY - FOV[0][1] / 2 - FOV[kk,1]
                options.bx[kk] = xx[0]
                options.bz[kk] = zz[0]
            elif (kk > 2 and kk < 5) or (FOV.shape[0] == 5 and kk > 0):
                if kk % 2 == 0:
                    options.bx[kk] = etaisyys[0,0] + options.oOffsetX + FOV[0,0]
                else:
                    options.bx[kk] = etaisyys[0,0] + options.oOffsetX - FOV[kk,0]
                options.by[kk] = yy[0]
                options.bz[kk] = zz[0]
            elif kk > 0 and kk < 3:
                options.bx[kk] = xx[0]
                options.by[kk] = yy[0]
                if kk % 2 == 0:
                    options.bz[kk] = options.oOffsetZ + FOV[0,2] / 2
                else:
                    options.bz[kk] = options.oOffsetZ - FOV[0][2] / 2 - FOV[kk,2]
    return xx, yy, zz


def computePixelCenters(options, xx, yy, zz):
    """
    This function computes the coordinates of the voxel centers.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    xx : NumPy array
        The coordinates of the voxel boundaries in x-direction.
    yy : NumPy array
        The coordinates of the voxel boundaries in y-direction.
    zz : NumPy array
        The coordinates of the voxel boundaries in z-direction.

    Returns
    -------
    None.

    """
    if options.projector_type in [2, 3, 22, 33, 12, 13, 21, 31, 42, 43, 24, 34]:
        options.x_center = xx[0 : -1] + options.dx[0] / 2.
        options.y_center = yy[0 : -1] + options.dy[0] / 2.
        options.z_center = zz[0 : -1] + options.dz[0] / 2.
    elif (options.projector_type == 1 and options.TOF):
        options.x_center = np.array(xx[0].item(), ndmin = 1)
        options.y_center = np.array(yy[0].item(), ndmin = 1)
        options.z_center = np.array(zz[0].item(), ndmin = 1)
    else:
        options.x_center = np.array(xx[0].item(), ndmin = 1)
        options.y_center = np.array(yy[0].item(), ndmin = 1)
        if zz.size == 0:
            options.z_center = np.zeros(1, dtype=np.float32)
        else:
            options.z_center = np.array(zz[0].item(), ndmin = 1)
    options.x_center = options.x_center.astype(dtype=np.float32)
    options.y_center = options.y_center.astype(dtype=np.float32)
    options.z_center = options.z_center.astype(dtype=np.float32)
            
def computeVoxelVolumes(options):
    """
    This is used by projector type 3 only. Computes a look-up table for the 
    volumes of intersection. Uses an analytic approach to compute the volumes.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    None.

    """
    from scipy.integrate import quad
    from scipy.special import ellipk, ellipe
    import math
    def volumeIntersection(R, r, b):
        def integrand(x, alpha, k):
            return 1. / ((1. + alpha * x**2) * np.sqrt(1. - x**2) * np.sqrt(1. - k * x**2))
        def vecquad(alpha, k):
            return quad(integrand, 0, 1, args=(alpha, k))[0]
        A = np.maximum(r**2, (b + R)**2)
        B = np.minimum(r**2, (b + R)**2)
        C = (b - R)**2
        k = np.abs((B - C) / (A - C))
        alpha = (B - C) / C
        s = (b + R) * (b - R)
        vec = np.vectorize(vecquad)
        Gamma = vec(alpha,k)
        K = ellipk(k)
        E = ellipe(k)
        heaviside_f = np.zeros(b.shape[0],dtype=np.float32)
        heaviside_f[R - b > 0] = 1
        heaviside_f[R - b == 0] = 1/2
        V = np.zeros(b.shape[0], dtype=np.float32)
        ind = r > b + R
        if np.any(ind):
            A1 = A[ind]
            B1 = B[ind]
            C1 = C[ind]
            s1 = s[ind]
            K1 = K[ind]
            E1 = E[ind]
            Gamma1 = Gamma[ind]
            V[ind] = (4. * np.pi)/3. * r**3 * heaviside_f[ind] + 4. / (3. * np.sqrt(A1 - C1)) * (Gamma1 * (A1**2 * s1) / C1 - K1 *(A1 * s1 - ((A1 - B1) * (A1 - C1)) / 3.) - E1 * (A1 - C1) * (s1 + (4. * A1 - 2. * B1 - 2. * C1) / 3.))
        ind = r < b + R
        if np.any(ind):
            A2 = A[ind]
            B2 = B[ind]
            C2 = C[ind]
            s2 = s[ind]
            K2 = K[ind]
            E2 = E[ind]
            Gamma2 = Gamma[ind]
            V[ind] = (4. * np.pi) / 3. * r**3 * heaviside_f[ind] + ((4. / 3.) / (np.sqrt(A2 - C2))) * (Gamma2 * ((B2**2 * s2) / C2) + K2 * (s2 * (A2 - 2. * B2) + (A2 - B2) * ((3. * B2 - C2 - 2 * A2) / 3.)) + E2 * (A2 - C2) * (-s2 + (2. * A2 - 4. * B2 + 2. * C2) / 3.))
        return V
    if options.projector_type in [3, 33, 13, 31, 34, 43]:
        dp = np.max(np.column_stack((options.dx[0].item(),options.dy[0].item(),options.dz[0].item())))
        options.voxel_radius = math.sqrt(2.) * options.voxel_radius * (dp / 2.)
        options.bmax = options.tube_radius + options.voxel_radius
        b = np.linspace(0, options.bmax, 10000, dtype=np.float32)
        b = b[(options.tube_radius <= (b + options.voxel_radius))]
        b = np.unique(np.round(b*10**3)/10**3)
        V = volumeIntersection(options.tube_radius, options.voxel_radius, b)
        diffis = np.append(np.diff(V), 0)
        # diffis = np.concatenate((np.diff(V),0))
        b = b[diffis <= 0]
        options.V = np.abs(V[diffis <= 0])
        options.Vmax = (4 * np.pi) / 3 * options.voxel_radius**3
        options.bmin = np.min(b)
        options.V = options.V.astype(dtype=np.float32)
    else:
        options.V = np.zeros(1, dtype=np.float32)
    
def computeProjectorScalingValues(options):
    """
    Computes scaling values for CT-based projectors (4 and 5). For example, 
    the interpolation length in projector type 4.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    None.

    """
    options.dScaleX4 = 1. / (options.dx * options.Nx)
    options.dScaleY4 = 1. / (options.dy * options.Ny)
    options.dScaleZ4 = 1. / (options.dz * options.Nz)
    if options.projector_type in [5, 15, 45, 54, 51]:
        options.dSizeY = 1. / (options.dy * options.Ny)
        options.dSizeX = 1. / (options.dx * options.Nx)
        options.dScaleX = 1. / (options.dx * (options.Nx + 1))
        options.dScaleY = 1. / (options.dy * (options.Ny + 1))
        options.dScaleZ = 1. / (options.dz * (options.Nz + 1))
        options.dSizeZBP = (options.nColsD + 1) * options.dPitchX
        options.dSizeXBP = (options.nRowsD + 1) * options.dPitchY
        options.dSizeY = options.dSizeY.astype(dtype=np.float32)
        options.dSizeX = options.dSizeX.astype(dtype=np.float32)
        options.dScaleX = options.dScaleX.astype(dtype=np.float32)
        options.dScaleY = options.dScaleY.astype(dtype=np.float32)
        options.dScaleZ = options.dScaleZ.astype(dtype=np.float32)
    # if options.projector_type == 4 or options.projector_type == 14 or options.projector_type == 54:
    if options.projector_type in [4, 14, 41, 42, 43, 44, 45, 54]:
        options.kerroin = (options.dx * options.dy * options.dz) / (options.dPitchX * options.dPitchY * options.sourceToDetector)
        if options.dL == 0:
            if isinstance(options.FOVa_x, np.ndarray):
                if isinstance(options.Nx, np.ndarray):
                    options.dL = options.FOVa_x[0].item() / options.Nx[0].item()
                else:
                    options.dL = options.FOVa_x[0].item() / options.Nx
            elif isinstance(options.Nx, np.ndarray):
                options.dL = options.FOVa_x / options.Nx[0].item()
            else:
                options.dL = options.FOVa_x / options.Nx
        else:
            if isinstance(options.FOVa_x, np.ndarray):
                if isinstance(options.Nx, np.ndarray):
                    options.dL = options.dL * (options.FOVa_x[0].item() / options.Nx[0].item())
                else:
                    options.dL = options.dL * (options.FOVa_x[0].item() / options.Nx)
            elif isinstance(options.Nx, np.ndarray):
                options.dL = options.dL * (options.FOVa_x / options.Nx[0].item())
            else:
                options.dL = options.dL * (options.FOVa_x / options.Nx)
    else:
        options.kerroin = np.zeros(1,dtype=np.float32)
    if options.projector_type in [4, 5, 14, 15, 45, 54] and options.CT:
        options.use_64bit_atomics = False
        options.use_32bit_atomics = False
    options.dScaleX4 = options.dScaleX4.astype(dtype=np.float32)
    options.dScaleY4 = options.dScaleY4.astype(dtype=np.float32)
    options.dScaleZ4 = options.dScaleZ4.astype(dtype=np.float32)