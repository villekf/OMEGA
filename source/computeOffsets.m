function options = computeOffsets(options)
%COMPUTEOFFSETS Computes the indices for the specified neighborhood
%   This function computes the neighborhood indices that are required by
%   L-filter and FMH, possibly also MRP if no medfilt3 is present
N = (options.Nx)*(options.Ny)*(options.Nz);
s = [options.Nx + options.Ndx*2 options.Ny + options.Ndy*2 options.Nz + options.Ndz*2];
N_pad = min(3, options.Ndx + options.Ndy + options.Ndz);
[c1{1:N_pad}] = ndgrid(1:(options.Ndx*2+1));
c2(1:N_pad) = {options.Ndy+1};
if options.Ndz > options.Ndx && options.Ndz > 1
    c1{1} = cat(3, c1{1}, zeros(size(c1{1},1), size(c1{1},2), options.Ndz));
    c1{2} = cat(3, c1{2}, zeros(size(c1{2},1), size(c1{2},2), options.Ndz));
    c1{3} = cat(3, c1{3}, zeros(size(c1{3},1), size(c1{3},2), options.Ndz));
    for kk = options.Ndz - 1 : - 1 : 0
        c1{1}(:,:,end-kk) = c1{1}(:,:,end - kk - 1);
        c1{2}(:,:,end-kk) = c1{2}(:,:,end - kk - 1);
        c1{3}(:,:,end-kk) = c1{3}(:,:,end - kk - 1) + 1;
    end
    c2(end) = {options.Ndz+1};
elseif options.Ndz < options.Ndx && options.Ndz > 1
    c1{1}(:,:,end-2*(options.Ndx-options.Ndz) + 1) = [];
    c1{2}(:,:,end-2*(options.Ndx-options.Ndz) + 1) = [];
    c1{3}(:,:,end-2*(options.Ndx-options.Ndz) + 1) = [];
    c2(end) = {options.Ndz+1};
end
offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
if options.Nz == 1
    tr_ind = sub2ind([options.Nx+options.Ndx*2 options.Ny + options.Ndy*2],mod((1:N)'-1,options.Nx)+(options.Ndx + 1),mod(floor(((1:double(N))'-1)/double(options.Nx)),...
        double(options.Ny))+(options.Ndy + 1));
else
    tr_ind = sub2ind([options.Nx+options.Ndx*2 options.Ny+options.Ndy*2 options.Nz+options.Ndz*2], mod((1:N)' - 1, options.Nx) + (options.Ndx + 1), ...
        mod(floor(((1:double(N))' - 1)/double(options.Nx)), double(options.Ny)) + (options.Ndy + 1), ...
        floor(((1:double(N))' - 1)/double(options.Nx * options.Ny)) + (options.Ndz + 1));
end
options.tr_offsets = uint32(bsxfun(@plus,tr_ind,offsets(:)'));
end

