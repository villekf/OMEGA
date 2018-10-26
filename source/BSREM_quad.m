function x_quad = BSREM_quad(options)
%BSREM_quad Block Sequential Regularized Expectation Maximization algorithm
%with quadratic prior reconstruction.
%   Implements the BSREM quadratic prior reconstruction on input PET data.
%   See main_nongate.m for options-variables.
%
%   x_quad = BSREM_quad(options) returns the BSREM quadratic prior 
%   reconstructions for all iterations, including the initial value.

if iscell(options.SinM)
    Sino = options.SinM{1};
    Sino = Sino(:);
else
    Sino = options.SinM;
    Sino = Sino(:);
end

LL = [];
index = [];
pituus = [];
lor = [];

if options.use_raw_data == false && options.subsets > 1
    if options.precompute_lor || options.reconstruction_method == 3
        load([options.machine_name '_lor_pixel_count_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor','discard')
        if length(discard) ~= options.TotSinos*options.Nang*options.Ndist
            error('Error: Size mismatch between sinogram and LORs to be removed')
        end
        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
            discard = discard(1:options.NSinos*options.Nang*options.Ndist);
        end
        ind_apu = uint32(find(discard));
        port = ceil((options.Nang-options.subsets+1)/options.subsets);
        over = options.Nang - port*options.subsets;
        index = cell(options.subsets,1);
        pituus = zeros(options.subsets, 1, 'uint32');
        for i=1:options.subsets
            if over>0
                index1 = uint32(sort(sub2ind([options.Nang options.Ndist options.NSinos],repmat(repelem(i:options.subsets:(port + 1)*options.subsets,options.Ndist)',options.NSinos,1),repmat((1:options.Ndist)',(port+1)*options.NSinos,1),repelem((1:options.NSinos)',options.Ndist*(port+1),1))));
                over = over - 1;
            else
                index1 = uint32(sort(sub2ind([options.Nang options.Ndist options.NSinos],repmat(repelem(i:options.subsets:port*options.subsets,options.Ndist)',options.NSinos,1),repmat((1:options.Ndist)',port*options.NSinos,1),repelem((1:options.NSinos)',options.Ndist*port,1))));
            end
            index{i} = index1(ismember(index1, ind_apu));
            pituus(i) = int32(length(index{i}));
        end
        index = cell2mat(index);
        index = index(ismember(index, ind_apu));
        clear index1 ind_apu
    else
        port = ceil((options.Nang-options.subsets+1)/options.subsets);
        over = options.Nang - port*options.subsets;
        index = cell(options.subsets,1);
        pituus = zeros(options.subsets, 1, 'uint32');
        for i=1:options.subsets
            if over>0
                index1 = uint32(sort(sub2ind([options.Nang options.Ndist options.NSinos],repmat(repelem(i:options.subsets:(port + 1)*options.subsets,options.Ndist)',options.NSinos,1),repmat((1:options.Ndist)',(port+1)*options.NSinos,1),repelem((1:options.NSinos)',options.Ndist*(port+1),1))));
                over = over - 1;
            else
                index1 = uint32(sort(sub2ind([options.Nang options.Ndist options.NSinos],repmat(repelem(i:options.subsets:port*options.subsets,options.Ndist)',options.NSinos,1),repmat((1:options.Ndist)',port*options.NSinos,1),repelem((1:options.NSinos)',options.Ndist*port,1))));
            end
            index{i} = uint32(index1);
            pituus(i) = int32(length(index1));
        end
        clear index1
    end
elseif options.subsets > 1
    % for raw list-mode data, take the options.subsets randomly
    % last subset has all the spare indices
    if options.precompute_lor || options.reconstruction_method == 3 || options.reconstruction_method == 2
        load([options.machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw.mat'],'LL','lor')
        indices = uint32(length(LL));
        index = cell(options.subsets, 1);
        port = uint32(floor(length(LL)/options.subsets));
        if options.use_Shuffle
            apu = Shuffle(indices(end), 'index')';
        else
            apu = uint32(randperm(indices(end)))';
        end
        pituus = zeros(options.subsets, 1, 'uint32');
        for i = 1 : options.subsets
            if i == options.subsets
                index{i} = apu(port*(i-1)+1:end);
            else
                index{i} = apu(port*(i-1)+1:(port*(i)));
            end
            pituus(i) = int32(length(index{i}));
        end
        clear apu
    else
        load([options.machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'LL')
        indices = uint32(length(LL));
        index = cell(options.subsets, 1);
        port = uint32(floor(length(LL)/options.subsets));
        if options.use_Shuffle
            apu = Shuffle(indices(end), 'index')';
        else
            apu = uint32(randperm(indices(end)))';
        end
        for i = 1 : options.subsets
            if i == options.subsets
                index{i} = apu(port*(i-1)+1:end);
            else
                index{i} = apu(port*(i-1)+1:(port*(i)));
            end
        end
        clear apu
    end
end

if options.precompute_lor && options.subsets > 1
    pituus2 = [0;cumsum(pituus)];
    Sino = Sino(index);
end

epps = options.epps;
weights = options.weights;
N = options.Nx*options.Ny*options.Nz;

x_quad = zeros(options.Nx,options.Ny,options.Nz, options.Niter + 1);
x_quad(:,:,:,1) = options.x0;
x_quad = reshape(x_quad, options.Nx*options.Ny*options.Nz, options.Niter + 1);

lam = zeros(options.Niter,1);
lam(1) = options.b0;
for i=1:options.Niter
    lam(i+1) = 0.5*lam(i);
end

distX = options.FOVa/double(options.Nx);
distZ = (double(options.axial_fov)/double(options.Nz));
if isempty(options.weights)
    for jj = options.Ndz : -1 : -options.Ndz
        cc = [];
        for kk = options.Ndy : -1 : -options.Ndy
            if options.Ndz == 0 || options.Nz == 1
                apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndx*2+1) * distX)'];
            else
                apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndx*2+1) * distX)', (repelem(jj,options.Ndx*2+1) * distZ)'];
            end
            edist = [];
            for ii = 1 : length(apu)
                edist = [edist;sqrt(apu(ii,:)*apu(ii,:)')];
            end
            cc = [cc;edist];
        end
        weights = [weights;cc];
    end
    weights = 1./weights;
end
pz_pad = padding(reshape(options.x0(:),options.Nx,options.Ny,options.Nz),[options.Ndx options.Ndy options.Ndz]);
s = size(pz_pad);
N_pad = min(3, options.Ndx + options.Ndy + options.Ndz);
[c1{1:N_pad}]=ndgrid(1:(options.Ndx*2+1));
c2(1:N_pad)={options.Ndy+1};
if options.Ndz > options.Ndx && options.Ndz > 1
    apu = c1{1};
    apu2 = c1{2};
    apu3 = c1{3};
    for kk = 1 : options.Ndz
        apu(:,:,end+1) = apu(:,:,end);
        apu2(:,:,end+1) = apu2(:,:,end);
        apu3(:,:,end+1) = apu3(:,:,end) + 1;
    end
    c1{1} = apu;
    c1{2} = apu2;
    c1{3} = apu3;
    c2(end) = {options.Ndz+1};
elseif options.Ndz < options.Ndx
    apu = c1{1};
    apu2 = c1{2};
    apu3 = c1{3};
    for kk = 1 : 2*(options.Ndx-options.Ndz)
        apu(:,:,end) = [];
        apu2(:,:,end) = [];
        apu3(:,:,end) = [];
    end
    c1{1} = apu;
    c1{2} = apu2;
    c1{3} = apu3;
    c2(end) = {options.Ndz+1};
end
offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
if options.Nz == 1
    tr_ind = sub2ind([options.Nx+options.Ndx*2 options.Ny+options.Ndy*2],mod((1:N)'-1,options.Nx)+(options.Ndx + 1),mod(floor(((1:double(N))'-1)/double(options.Nx)),double(options.Ny))+(options.Ndy + 1));
else
    tr_ind = sub2ind([options.Nx+options.Ndx*2 options.Ny+options.Ndy*2 options.Nz+options.Ndz*2],mod((1:N)'-1,options.Nx)+(options.Ndx + 1),mod(floor(((1:double(N))'-1)/double(options.Nx)),double(options.Ny))+(options.Ndy + 1),floor(((1:double(N))'-1)/double(options.Nx*options.Ny))+(options.Ndz+1));
end
tr_offsets = bsxfun(@plus,tr_ind,offsets(:)');

weights = weights/sum(weights(~isinf(weights)));
weights_quad = [weights(1:floor(length(weights)/2));weights(ceil(length(weights)/2) + 1:end)];

for ii = 1 : options.Niter
    quad_apu = x_quad(:,ii);
    for kk = 1 : options.subsets
        [A] = observation_matrix_formation_nongate(options, kk, index, LL, pituus, lor);
        
        if options.precompute_lor == false
            uu = double(Sino(index{kk}));
        else
            uu = double(Sino(pituus2(kk)+1:pituus2(kk + 1)));
        end
        tStart = tic;
        quad_apu = quad_apu + lam(ii).*quad_apu./(full(sum(A,1)')+epps).*(A'*(uu./(A*quad_apu+epps))-full(sum(A,1)'));
        tElapsed = toc(tStart);
        disp(['BSREM quadratic prior sub-iteration ' num2str(kk) ' took ' num2str(tElapsed) ' seconds'])
        disp(['BSREM quadratic prior sub-iteration ' num2str(kk) ' finished'])
    end
    tStart = tic;
    if options.Nz==1
        pz_pad_bsrem = padding(reshape(quad_apu,options.Nx,options.Ny),[options.Ndx options.Ndy]);
        med = -bsxfun(@minus,pz_pad_bsrem(tr_offsets(:,isinf(weights))),pz_pad_bsrem(tr_offsets(:,[1:(find(isinf(weights))-1) (find(isinf(weights))+1):end])))*weights_quad;
    else
        pz_pad_bsrem = padding(reshape(quad_apu,options.Nx,options.Ny,options.Nz),[options.Ndx options.Ndy options.Ndz]);
        med = -bsxfun(@minus,pz_pad_bsrem(tr_offsets(:,isinf(weights))),pz_pad_bsrem(tr_offsets(:,[1:(find(isinf(weights))-1) (find(isinf(weights))+1):end])))*weights_quad;
    end
    x_quad(:,ii+1) = quad_apu.*(1 + options.beta_quad_bsrem*lam(ii).*med);
    x_quad(x_quad(:,ii+1)<0,ii+1) = 0;
    tElapsed = toc(tStart);
    disp(['BSREM quadratic prior iteration ' num2str(ii) ' took ' num2str(tElapsed) ' seconds'])
end
x_quad = reshape(x_quad,options.Nx,options.Ny,options.Nz, options.Niter + 1);

end