# -*- coding: utf-8 -*-
import numpy as np
   
def indexMaker(options):
    subsets = options.subsets;
    Ndist = options.Ndist;
    Nang = options.Nang;
    NSinos = options.NSinos;
    options.nMeas = np.array(options.Ndist * options.Nang * options.NSinos, dtype=np.int64)
    if options.CT:
        tyyppi = np.uint64;
    else:
        tyyppi = np.uint32;
    if subsets > 1 and options.subsetType < 8:
        totalLength = Ndist*Nang*NSinos
        options.index = np.empty(0, dtype=tyyppi)
        options.nMeas = np.zeros((subsets, 1), dtype = np.int64)
        if options.subsetType == 4 and options.use_raw_data == 0:
            for i in range(subsets):
                osa = np.size(np.arange(i, Nang, subsets))
                index1 = np.tile(np.arange(i*Ndist, (i + 1) * Ndist, 1), (osa*NSinos)).astype(tyyppi)
                index1 = index1 + np.repeat(np.arange(0, (osa)*NSinos, 1) * Ndist*subsets,Ndist).astype(tyyppi)
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
        elif options.subsetType == 5:
            apu = np.tile((np.tile(np.arange(0, Nang * Ndist, Ndist, dtype=tyyppi).T[None,:], (Ndist, 1)) + np.arange(0, Ndist, dtype=tyyppi)[:,None])[:,:,None], (1, 1, NSinos)) + np.transpose(np.arange(0, Nang * Ndist * NSinos, Nang * Ndist, dtype=tyyppi)[:,None,None],(2, 1, 0))
            apu = np.transpose(apu, (0, 2, 1)).reshape(Ndist * NSinos, -1, order='F')
            for i in range(subsets):
                index1 = apu[i : : subsets, :]
                index1 = index1.T
                index1 = index1.ravel('F')
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
            # Take every nth (column) measurement
        elif options.subsetType == 1:
            for i in range(subsets):
                index1 = np.arange(i, totalLength, subsets).astype(tyyppi)
                [I,J,K] = options.ind2sub(index1, np.array((Nang, Ndist, NSinos)))
                index1 = options.sub2ind(np.array((Ndist, Nang, NSinos)), J, I,K)
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
            # Take every nth (row) measurement
            # Every nth measurements
        elif options.subsetType == 2:
            for i in range(subsets):
                index1 = np.arange(i, totalLength, subsets).astype(tyyppi)
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
            # Pick the measurements randomly
        elif options.subsetType == 3:
            port = totalLength // subsets
            generator = np.random.default_rng()
            apu = generator.permutation(totalLength).astype(tyyppi)
            for i in range(subsets):
                if i == subsets:
                    index1 = apu[port * i : -1]
                else:
                    index1 = apu[port * i : port * (i + 1)]
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
            # Pick the subsets based on the angles of the LORs
        elif options.subsetType == 6:
            raise ValueError('Not supported in Python version')
            # Use golden angle sampling
        elif options.subsetType == 7 or options.use_raw_data == 0:
            raise ValueError('Not supported in Python version')
        elif options.subsetType == 0:
            val = totalLength // subsets
            if totalLength % subsets > 0:
                valEnd = totalLength - val * (subsets - 1)
            else:
                valEnd = val
            options.nMeas = np.full(subsets - 1, val, dtype=np.int64)
            options.nMeas = np.append(options.nMeas, valEnd)
            options.index = np.zeros(1,dtype=tyyppi)
    elif (subsets > 1 and (options.subsetType == 8 or options.subsetType == 9 or options.subsetType == 10 or options.subsetType == 11)) or subsets == 1:
        sProjections = options.nProjections // subsets
        modi = np.mod(options.nProjections, subsets)
        uu = (modi > 0).astype(np.int64)
        ind1 = 0
        ind2 = sProjections + uu
        options.index = np.empty(0, dtype=tyyppi)
        options.nMeas = np.zeros((subsets, 1), dtype = np.int64)
        if options.subsetType == 8 and subsets > 1:
            for kk in range(subsets):
                index1 = np.arange(kk, options.nProjections, subsets).astype(tyyppi)
                options.index = np.append(options.index, index1)
                options.nMeas[kk] = np.size(index1)
        elif options.subsetType == 9 and subsets > 1:
            generator = np.random.default_rng()
            apu = generator.permutation(options.nProjections).astype(tyyppi)
            for kk in range(subsets):
                index1 = apu[ind1 : ind2]
                options.index = np.append(options.index, index1)
                options.nMeas[kk] = np.size(index1)
                modi = modi - 1
                if modi <= 0:
                    uu = 0
                ind1 = ind2
                ind2 = ind2 + (options.nProjections // subsets) + uu
        elif options.subsetType == 10 and subsets > 1:
            raise ValueError('Not supported in Python version!')
            # ga = 2.39996322972865332
            # options.angles = np.abs(options.angles)
            # angles = options.angles - np.min(options.angles)
            # anglesOrig = angles
            # maksimi = np.max(angles)
            # ga = ga * (maksimi / (2.*np.pi))
            # angle = 0.
            # for kk in range(subsets - 1):
            #     ind = np.zeros(options.angles.size // subsets,1, dtype = tyyppi)
            #     for ii in range(options.angles.size // subsets):
            #         I = np.argmin(np.abs(angles-angle))
            #         II = numpy.where(np.isin(anglesOrig,angles[I[0]]) > 0)
            #         angles[I] = []
            #         ind[ii] = II
            #         angle = angle + ga
            #         if angle > maksimi:
            #             angle = angle - maksimi
            #     index{kk} = ind
            #     options.nMeas(kk) = numel(index{kk})
            # ind = np.zeros(ceil(numel(options.angles) / subsets),1)
            # for ii = 1 : ceil(numel(options.angles) / subsets)
            #     [~,I] = min(abs(angles-angle))
            #     II = find(ismember(anglesOrig,angles(I)))
            #     angles(I) = []
            #     ind(ii) = II
            #     angle = angle + ga
            #     if angle > maksimi:
            #         angle = angle - maksimi
            # index{subsets} = ind
            # options.nMeas(subsets) = numel(index{subsets})
        elif options.subsetType == 11 and subsets > 1:
            def factor(n):
                factors = []
                divisor = 2
                
                while n > 1:
                    while n % divisor == 0:
                        factors.append(divisor)
                        n //= divisor
                    divisor += 1
                
                return factors
            n = factor(options.nProjections)
            n2 = np.array(n)
            nn = np.flipud(np.cumprod(np.flipud(n2[1:])))
            N = np.size(n)
            p = np.zeros((N, options.nProjections), dtype=tyyppi)
            for ll in range(1, N + 1):
                p1 = np.arange(n[ll - 1])
                if ll == 1:
                    p1 = np.tile(p1, options.nProjections // np.size(p1))
                else:
                    p1 = np.repeat(p1, np.prod(n[:ll - 1]))
                    p1 = np.tile(p1, options.nProjections // np.size(p1))
                p[ll - 1, :] = p1
            
            indices = np.zeros(options.nProjections, dtype=tyyppi)
            for r in range(1, options.nProjections + 1):
                tt = p[:len(nn), r - 1] * nn
                indices[r - 1] = np.sum(tt) + p[-1, r - 1]
            for kk in range(subsets):
                options.index = np.append(options.index, indices[ind1:ind2])
                options.nMeas[kk] = np.size(indices[ind1:ind2])
                modi -= 1
                if modi <= 0:
                    uu = 0
                ind1 = ind2
                ind2 = ind2 + (options.nProjections // subsets) + uu
        else:
            options.index = np.arange(0, options.nProjections).astype(tyyppi)
            options.nMeas[0] = options.index.size
        if options.subsets == 1:
            if options.CT or options.PET:
                options.nMeas[0] = options.NSinos
            else:
                options.nMeas[0] = options.Nang * options.Ndist * options.NSinos
    elif options.subsetType > 11:
        raise ValueError('Invalid subset type!')
    options.subsets = subsets
    
def ind2sub(options, index, dims):
    K = index // (dims[0] * dims[1])
    J = (index - K * (dims[0] * dims[1])) // dims[0]
    I = np.mod(index - K * (dims[0] * dims[1]), dims[0])
    return (I, J, K)
    
def sub2ind(options, dims, I, J, K):
    index = K * dims[0] * dims[1] + I + J * dims[0]
    return index


def computePixelSize(options):
    FOV = np.column_stack((options.FOVa_x,options.FOVa_y, options.axial_fov)).astype(dtype=np.float32)
    etaisyys = -(FOV) / 2
    options.dx = np.zeros((FOV.shape[0],1),dtype=np.float32);
    options.dy = np.zeros((FOV.shape[0],1),dtype=np.float32);
    options.dz = np.zeros((FOV.shape[0],1),dtype=np.float32);
    options.bx = np.zeros((FOV.shape[0],1),dtype=np.float32);
    options.by = np.zeros((FOV.shape[0],1),dtype=np.float32);
    options.bz = np.zeros((FOV.shape[0],1),dtype=np.float32);
    for kk in range(FOV.shape[0] - 1, -1, -1):
        xx = np.linspace(etaisyys[kk,0] + options.oOffsetX, -etaisyys[kk,0] + options.oOffsetX, options.Nx[kk].item() + 1, dtype=np.float32)
        yy = np.linspace(etaisyys[kk,1] + options.oOffsetY, -etaisyys[kk,1] + options.oOffsetY, options.Ny[kk].item() + 1, dtype=np.float32)
        zz = np.linspace(etaisyys[kk,2] + options.oOffsetZ, -etaisyys[kk,2] + options.oOffsetZ, options.Nz[kk].item() + 1, dtype=np.float32)
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
                    options.by[kk] = etaisyys[0,1] + options.oOffsetY + FOV[0,1];
                else:
                    options.by[kk] = etaisyys[0,1] + options.oOffsetY - FOV[kk,1];
                options.bx[kk] = xx[0];
                options.bz[kk] = zz[0];
            elif (kk > 2 and kk < 5) or (FOV.shape[0] == 5 and kk > 0):
                if kk % 2 == 0:
                    options.bx[kk] = etaisyys[0,0] + options.oOffsetX + FOV[0,0];
                else:
                    options.bx[kk] = etaisyys[0,0] + options.oOffsetX - FOV[kk,0];
                options.by[kk] = yy[0];
                options.bz[kk] = zz[0];
            elif kk > 0 and kk < 3:
                options.bx[kk] = xx[0];
                options.by[kk] = yy[0];
                if kk % 2 == 0:
                    options.bz[kk] = etaisyys[0,2] + options.oOffsetZ + FOV[0,2];
                else:
                    options.bz[kk] = etaisyys[0,2] + options.oOffsetZ - FOV[kk,2];
    return xx, yy, zz

def formSubsetIndices(options):
    if options.listmode > 0 and options.subsetType < 8 and options.subsets > 1:
        if options.subsetType > 0:
            if not (options.x.shape[0] == 6):
                options.x = np.reshape(options.x, (6, options.x.size // 6))
            options.x = options.x[:,options.index]
        options.x = options.x.ravel('F')
        options.xy_index = np.empty(0, dtype=np.uint32)
        options.z_index = np.empty(0, dtype=np.uint16)
    elif options.listmode > 0 and (options.subsetType == 8 or options.subsetType == 9) and options.subsets > 1:
        options.z = options.z[options.index,:]
        if options.uV.shape[0] == options.nProjections:
            options.uV = options.uV[options.index,:]
        options.xy_index = np.empty(0, dtype=np.uint32)
        options.z_index = np.empty(0, dtype=np.uint16)
    elif options.use_raw_data == 0 and ((options.subsets > 1 and (options.subsetType == 3 or options.subsetType == 6 or options.subsetType == 7))):
        options.xy_index = np.arange(0, options.Nang * options.Ndist, dtype=np.uint32)
        if options.span > 1:
            xy_index2 = np.tile(np.arange(0, options.Nang * options.Ndist, dtype=np.uint32), (options.NSinos - (options.rings * 2 - 1)))
            options.xy_index = np.concatenate((np.tile(options.xy_index,  (options.rings * 2 - 1)), xy_index2))
        else:
            xy_index2 = np.tile(np.arange(0, options.Nang * options.Ndist, dtype=np.uint32), (options.NSinos - options.rings))
            options.xy_index = np.concatenate((np.tile(options.xy_index, (options.rings)), xy_index2))
        options.z_index = np.arange(0, options.NSinos, dtype=np.uint16)
        options.z_index = np.repeat(options.z_index, (options.Nang * options.Ndist))
        options.z_index = options.z_index[options.index]
    
        options.xy_index = options.xy_index[options.index]
    else:
        options.xy_index = np.empty(0, dtype=np.uint32)
        options.z_index = np.empty(0, dtype=np.uint16)
    if options.sampling > 1 and not options.use_raw_data and not options.precompute_lor:
        options.Ndist = options.Ndist / options.sampling;
        
def computePixelCenters(options, xx, yy, zz):
    if (options.projector_type == 2 or options.projector_type == 3 or options.projector_type == 22 or options.projector_type == 33 or options.projector_type == 12 
        or options.projector_type == 13 or options.projector_type == 21 or options.projector_type == 31 or options.projector_type == 42 or options.projector_type == 43 
        or options.projector_type == 24 or options.projector_type == 34):
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
    from scipy.integrate import quad
    from scipy.special import ellipk, ellipe
    import math
    def volumeIntersection(R, r, b):
        def integrand(x, alpha, k):
            return 1. / ((1. + alpha * x**2) * np.sqrt(1. - x**2) * np.sqrt(1. - k * x**2))
        def vecquad(alpha, k):
            return quad(integrand, 0, 1, args=(alpha, k))[0]
        A = np.maximum(r**2, (b + R)**2);
        B = np.minimum(r**2, (b + R)**2);
        C = (b - R)**2;
        k = np.abs((B - C) / (A - C));
        alpha = (B - C) / C;
        s = (b + R) * (b - R);
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
    if (options.projector_type == 3 or options.projector_type == 33 or options.projector_type == 13 or options.projector_type == 31 or options.projector_type == 34 
        or options.projector_type == 43):
        dp = np.max(np.column_stack((options.dx[0].item(),options.dy[0].item(),options.dz[0].item())))
        options.voxel_radius = math.sqrt(2.) * options.voxel_radius * (dp / 2.)
        options.bmax = options.tube_radius + options.voxel_radius
        b = np.linspace(0, options.bmax, 10000, dtype=np.float32)
        b = b[(options.tube_radius <= (b + options.voxel_radius))]
        b = np.unique(np.round(b*10**3)/10**3);
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
    options.dScaleX4 = 1. / (options.dx * options.Nx)
    options.dScaleY4 = 1. / (options.dy * options.Ny)
    options.dScaleZ4 = 1. / (options.dz * options.Nz)
    if options.projector_type == 5 or options.projector_type == 15 or options.projector_type == 45 or options.projector_type == 54 or options.projector_type == 51:
        options.dSizeY = 1. / (options.dy * options.Ny)
        options.dSizeX = 1. / (options.dx * options.Nx)
        options.dScaleX = 1. / (options.dx * (options.Nx + 1))
        options.dScaleY = 1. / (options.dy * (options.Ny + 1))
        options.dScaleZ = 1. / (options.dz * (options.Nz + 1))
        options.dSizeZBP = (options.nColsD + 1) * options.dPitchX
        options.dSizeXBP = (options.nRowsD + 1) * options.dPitchY
        options.dSizeY = options.dSizeY.astype(dtype=np.float32)
        options.dSizeX = options.dSizeY.astype(dtype=np.float32)
        options.dScaleX = options.dScaleX.astype(dtype=np.float32)
        options.dScaleY = options.dScaleY.astype(dtype=np.float32)
        options.dScaleZ = options.dScaleZ.astype(dtype=np.float32)
    if options.projector_type == 4 or options.projector_type == 14 or options.projector_type == 54:
        options.kerroin = (options.dx * options.dy * options.dz) / (options.dPitchX * options.dPitchY * options.sourceToDetector)
        if options.dL == 0:
            options.dL = options.FOVa_x / options.Nx
        else:
            options.dL = options.dL * (options.FOVa_x / options.Nx)
    else:
        options.kerroin = np.zeros(1,dtype=np.float32);
    if ((options.projector_type == 4 or options.projector_type == 5 or options.projector_type == 14 or options.projector_type == 15 or options.projector_type == 45 
            or options.projector_type == 54) and options.CT):
        options.use_64bit_atomics = False
        options.use_32bit_atomics = False
    options.dScaleX4 = options.dScaleX4.astype(dtype=np.float32)
    options.dScaleY4 = options.dScaleY4.astype(dtype=np.float32)
    options.dScaleZ4 = options.dScaleZ4.astype(dtype=np.float32)