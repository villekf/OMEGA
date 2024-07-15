function [m_size, xy_index_input, z_index_input, L_input, lor_input, lor2, norm_input, scatter_input] = splitInput(options, pituus, ll, xy_index, z_index, LL, lor_a, ScatterC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (options.subset_type >= 8 || options.subsets == 1) && options.listmode == 0 && options.projector_type ~= 6
    kerroin = options.nRowsD * options.nColsD;
else
    kerroin = 1;
end
% if ((options.CT || options.SPECT || options.PET) && options.listmode == 0)
%     m_size = options.xSize * options.ySize * pituus(ll + 1);
% else
if options.subsets == 1 && options.projector_type == 6
    m_size = options.nProjections;
else
    m_size = pituus(ll + 1) * kerroin - pituus(ll) * kerroin;
end
% end
if ~options.listmode && ((options.subset_type == 3 || options.subset_type == 6 || options.subset_type == 7) && options.subsets > 1)
    xy_index_input = xy_index(pituus(ll) * kerroin + 1: pituus(ll + 1) * kerroin);
    z_index_input = z_index(pituus(ll) * kerroin + 1: pituus(ll + 1) * kerroin);
else
    xy_index_input = uint32([]);
    z_index_input = uint16([]);
end
if options.use_raw_data
    L_input = LL(pituus(ll) * kerroin + 1 : pituus(ll + 1) * kerroin,:);
else
    L_input = uint16([]);
end
if options.implementation == 1
    if options.subset_type >= 8 && options.subsets > 1
        lor_input = lor_a(:,:,pituus(ll) + 1 : pituus(ll + 1));
    elseif options.subsets > 1
        lor_input = lor_a(pituus(ll) * kerroin + 1 : pituus(ll + 1) * kerroin);
    else
        lor_input = lor_a;
    end
    lor_input = reshape(lor_input, [], 1);
    if options.implementation == 1
        if exist('OCTAVE_VERSION','builtin') == 0
            lor2 = [uint64(0);cumsum(uint64(lor_input))];
        else
            lor2 = [uint64(0);cumsum(uint64(lor_input),'native')];
        end
    end
else
    lor2 = uint64([]);
    lor_input = uint16([]);
end
if options.normalization_correction
    norm_input = cast(options.normalization(pituus(ll) * kerroin + 1 : pituus(ll + 1) * kerroin), options.cType);
else
    norm_input = cast(0, options.cType);
end
if options.additionalCorrection
    scatter_input = cast(ScatterC(pituus(ll) * kerroin + 1 : pituus(ll + 1) * kerroin), options.cType);
else
    scatter_input = cast(0, options.cType);
end