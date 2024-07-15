# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:45:47 2024

@author: Ville-Veikko Wettenhovi
"""

def applyMeasPreconditioning(options, var):
    if (options.precondTypeMeas[0].item() or options.precondTypeMeas[1].item()):
        if options.useAF:
            import arrayfire as af
            if options.precondTypeMeas[1].item():
                if not hasattr(options, 'filterG'):
                    options.filterG = af.interop.np_to_af_array(options.filter0)
                if (options.subsets > 1 and (options.subsetType == 5 or options.subsetType == 4)):
                	if (options.subsetType == 4):
                		var = af.moddims(var, options.nRowsD, d1=var.elements() // options.nRowsD);
                	else:
                		var = af.moddims(var, options.nColsD, d1=var.elements() // options.nColsD);
                else:
                	var = af.moddims(var, options.nRowsD, d1=options.nColsD, d2=var.elements() // (options.nRowsD * options.nColsD));
                temp = af.fft(var, options.Nf)
                temp = temp * af.tile(options.filterG, 1, d1=temp.shape[1], d2=temp.shape[2])
                af.eval(temp)
                af.ifft_inplace(temp)
                var = af.flat(af.real(temp[:var.shape[0], :, :]))
        elif options.useTorch:
            # def reshape_fortran(x, shape):
            #     if len(x.shape) > 0:
            #         x = x.permute(*reversed(range(len(x.shape))))
            #     return x.reshape(*reversed(shape)).permute(*reversed(range(len(shape))))
            import torch
            if options.precondTypeMeas[1].item():
                if not hasattr(options, 'filterG'):
                    import numpy as np
                    options.filter0 = np.reshape(options.filter0, (1, 1, -1))
                    options.filterG = torch.tensor(options.filter0, device='cuda')
                if (options.subsets > 1 and (options.subsetType == 5 or options.subsetType == 4)):
                	if (options.subsetType == 4):
                		var = torch.reshape(var, (var.numel() // options.nRowsD, options.nRowsD));
                	else:
                		var = torch.reshape(var, (var.numel() // options.nColsD), options.nColsD);
                else:
                    # var = reshape_fortran(var, (options.nRowsD, options.nColsD, var.numel() // (options.nRowsD * options.nColsD)))
                    var = torch.reshape(var, (options.nColsD, var.numel() // (options.nRowsD * options.nColsD), options.nRowsD));
                temp = torch.fft.fft(var, n=options.Nf, dim=2)
                temp = temp * options.filterG
                temp = torch.fft.ifft(temp, dim=2)
                var = torch.ravel(torch.real(temp[:,:,:var.shape[2]]))
    return var
            
def circulantInverse(options, var):
    if options.useAF:
        import arrayfire as af
        if not hasattr(options, 'FilterG'):
            options.FilterG = af.interop.np_to_af_array(options.Ffilter)
        if (options.subsets > 1 and (options.subsetType == 5 or options.subsetType == 4)):
        	if (options.subsetType == 4):
        		var = af.moddims(var, options.nRowsD, d1=var.elements() // options.nRowsD);
        	else:
        		var = af.moddims(var, options.nColsD, d1=var.elements() // options.nColsD);
        else:
        	var = af.moddims(var, options.nRowsD, d1=options.nColsD, d2=var.elements() // (options.nRowsD * options.nColsD));
        temp = af.fft(var, options.Nf)
        temp /= af.tile(options.FilterG, 1, d1=temp.shape[1], d2=temp.shape[2])
        af.eval(temp)
        af.ifft_inplace(temp)
        var = af.flat(af.real(temp[:var.shape[0], :, :]))
    elif options.useTorch:
        import torch
        if not hasattr(options, 'FilterG'):
            import numpy as np
            options.Ffilter = np.reshape(options.Ffilter, (1, 1, -1))
            options.FilterG = torch.tensor(options.Ffilter, device='cuda')
        if (options.subsets > 1 and (options.subsetType == 5 or options.subsetType == 4)):
        	if (options.subsetType == 4):
        		var = torch.reshape(var, (options.nRowsD, var.numel() // options.nRowsD));
        	else:
        		var = torch.reshape(var, (options.nColsD, var.numel() // options.nColsD));
        else:
        	var = torch.reshape(var, (options.nColsD, var.numel() // (options.nRowsD * options.nColsD), options.nRowsD));
        temp = torch.fft.fft(var, n=options.Nf, dim=2)
        temp /= options.FilterG
        temp = torch.fft.ifft(temp, dim=2)
        var = torch.ravel(torch.real(temp[:,:,:var.shape[2]]))
    return var