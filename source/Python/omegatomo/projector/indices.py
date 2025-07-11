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
            maksimi = np.size(np.arange(1, Nang, subsets))
            for i in range(subsets):
                osa = np.size(np.arange(i, Nang, subsets))
                index1 = np.tile(np.arange(i*Ndist, (i + 1) * Ndist, 1), (osa*NSinos)).astype(tyyppi)
                index1 = index1 + np.repeat(np.arange(0, (osa)*NSinos, 1) * Ndist*subsets,Ndist).astype(tyyppi)
                if Nang % subsets > 0:
                    if osa < maksimi:
                        erotus = osa - 1;
                    else:
                        erotus = Nang % subsets - subsets
                    index1 = (np.int64(index1) + np.int64(np.repeat((np.arange(0, NSinos - 1, 1)) * Ndist * erotus, Ndist * osa))).astype(tyyppi)
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
                [I,J,K] = ind2sub(index1, np.array((Nang, Ndist, NSinos)))
                index1 = sub2ind(np.array((Ndist, Nang, NSinos)), J, I,K)
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
        elif options.subsetType == 7:
            raise ValueError('Not supported in Python version')
        elif options.subsetType == 0:
            if options.listmode:
                val = totalLength // subsets
                if totalLength % subsets > 0:
                    valEnd = totalLength - val * (subsets - 1)
                else:
                    valEnd = val
            else:
                val = options.nProjections // subsets
                if totalLength % subsets > 0:
                    valEnd = options.nProjections - val * (subsets - 1)
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
            if (options.CT or options.PET or options.SPECT) and options.listmode == 0:
                options.nMeas[0] = options.NSinos
            elif options.listmode == 1:
                if options.useIndexBasedReconstruction:
                    options.nMeas[0] = options.trIndex.size // 2
                else:
                    options.nMeas[0] = options.x.size // 6
            else:
                options.nMeas[0] = options.Nang * options.Ndist * options.NSinos
    elif options.subsetType > 11:
        raise ValueError('Invalid subset type!')
    options.subsets = subsets
    
def ind2sub(index, dims):
    K = index // (dims[0] * dims[1])
    J = (index - K * (dims[0] * dims[1])) // dims[0]
    I = np.mod(index - K * (dims[0] * dims[1]), dims[0])
    return (I, J, K)
    
def sub2ind(dims, I, J, K):
    index = K * dims[0] * dims[1] + I + J * dims[0]
    return index


def formSubsetIndices(options):
    if options.listmode > 0 and options.subsetType < 8 and options.subsets > 1:
        if options.subsetType > 0:
            if options.useIndexBasedReconstruction:
                if options.trIndex.shape[0] != 2:
                    options.trIndex = np.reshape(options.trIndex, (2, -1), order='F')
                if options.Nt > 1:
                    options.trIndex = options.trIndex[:,options.index,:]
                else:
                    options.trIndex = options.trIndex[:,options.index]
                if options.axIndex.shape[0] != 2:
                    options.axIndex = np.reshape(options.axIndex, (2, -1), order='F')
                if options.Nt > 1:
                    options.axIndex = options.axIndex[:,options.index,:]
                else:
                    options.axIndex = options.axIndex[:,options.index]
            else:
                if not (options.x.shape[0] == 6):
                    options.x = np.reshape(options.x, (6, options.x.size // 6))
                if options.Nt > 1:
                    options.x = options.x[:,options.index,:]
                else:
                    options.x = options.x[:,options.index]
            if options.TOF_bins_used > 1:
                if options.Nt > 1:
                    options.TOFIndices = options.TOFIndices[options.index,:]
                else:
                    options.TOFIndices = options.TOFIndices[options.index]
        options.x = options.x.ravel('F')
        options.xy_index = np.empty(0, dtype=np.uint32)
        options.z_index = np.empty(0, dtype=np.uint16)
    elif options.listmode > 0 and (options.subsetType == 8 or options.subsetType == 9) and options.subsets > 1:
        if options.Nt > 1:
            options.z = options.z[options.index,:,:]
            if options.uV.shape[0] == options.nProjections:
                options.uV = options.uV[options.index,:,:]
        else:
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
        