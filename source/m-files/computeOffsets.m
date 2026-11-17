function options = computeOffsets(options)
%COMPUTEOFFSETS Computes the indices for the specified neighborhood
%   This function computes the neighborhood indices that are required by
%   L-filter and FMH, possibly also MRP if no medfilt3 is present
N = (options.Nx)*(options.Ny)*(options.Nz);
s = [options.Nx + options.Ndx*2 options.Ny + options.Ndy*2 options.Nz + options.Ndz*2];
% Full neighborhood box, x fastest (column-major), centred at (Ndx+1, Ndy+1, Ndz+1).
% The number of offsets is exactly (2*Ndx+1)*(2*Ndy+1)*(2*Ndz+1).
[c1{1:3}] = ndgrid(1:(options.Ndx*2+1), 1:(options.Ndy*2+1), 1:(options.Ndz*2+1));
offsets = sub2ind(s, c1{1}, c1{2}, c1{3}) - sub2ind(s, options.Ndx + 1, options.Ndy + 1, options.Ndz + 1);
if options.Nz == 1
    tr_ind = sub2ind([options.Nx+options.Ndx*2 options.Ny + options.Ndy*2],mod((1:N)'-1,options.Nx)+(options.Ndx + 1),mod(floor(((1:double(N))'-1)/double(options.Nx)),...
        double(options.Ny))+(options.Ndy + 1));
else
    tr_ind = sub2ind([options.Nx+options.Ndx*2 options.Ny+options.Ndy*2 options.Nz+options.Ndz*2], mod((1:N)' - 1, options.Nx) + (options.Ndx + 1), ...
        mod(floor(((1:double(N))' - 1)/double(options.Nx)), double(options.Ny)) + (options.Ndy + 1), ...
        floor(((1:double(N))' - 1)/double(options.Nx * options.Ny)) + (options.Ndz + 1));
end
options.tr_offsets = uint32(bsxfun(@plus,tr_ind,offsets(:)'));
if options.implementation == 2
    if options.MRP || options.L || options.FMH || (options.TV && options.TVtype == 3) || options.RDP
        options.tr_offsets = options.tr_offsets - 1;
    end
    options.Ndx = uint32(options.Ndx);
    options.Ndy = uint32(options.Ndy);
    options.Ndz = uint32(options.Ndz);
end
if options.MRP
    options.medx = options.Ndx*2 + 1;
    options.medy = options.Ndy*2 + 1;
    options.medz = options.Ndz*2 + 1;
end
end

