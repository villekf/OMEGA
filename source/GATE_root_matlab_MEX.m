function [L1, L2, tpoints, A, int_loc, Ltrues, Lscatter, Lrandoms, trues_index, Ldelay1, Ldelay2, int_loc_delay, tpoints_delay, randoms_index, scatter_index, ...
    x1, x2, y1, y2, z1, z2, time, timeD] = GATE_root_matlab_MEX(nimi,vali,alku,loppu, detectors, blocks_per_ring, cryst_per_block, det_per_ring, source, time_intervals, ...
    options, scatter_components, store_coordinates, transaxial_multip, cryst_per_block_z, rings, large_case)
%GATE_ROOT_MATLAB_MEX Computes the ROOT data load utilizing mexhost
%   In R2019a and up, a mex-file can be computed separate from MATLAB
%   preventing any possible library conflicts. This can be achieved by
%   using mexhost and feval.
% persistent mh
% if ~(isa(mh,'matlab.mex.MexHost') && isvalid(mh))
    mh = mexhost;
% end
[L1, L2, tpoints, A, int_loc, Ltrues, Lscatter, Lrandoms, trues_index, Ldelay1, Ldelay2, int_loc_delay, tpoints_delay, randoms_index, scatter_index, x1, x2, y1, y2, z1, z2, time, timeD] = ...
    feval(mh, 'GATE_root_matlab', nimi,vali,alku,loppu, uint32(detectors), blocks_per_ring, cryst_per_block, det_per_ring, uint32(options.linear_multip), ...
    source, time_intervals, options.obtain_trues, options.store_scatter, options.store_randoms, scatter_components, options.randoms_correction, store_coordinates, ...
    transaxial_multip, cryst_per_block_z, rings, large_case);
clear mh mex
end