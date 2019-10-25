function x_mramla = MRAMLA(options)
%MRAMLA Modified Row-Action Maximum Likelihood Algorithm reconstruction.
%   Implements the MRAMLA reconstruction on input PET data.
%   See main_nongate.m for options-variables.
%
%   x_mramla = MRAMLA(options) returns the MRAMLA reconstructions for all
%   iterations, including the initial value.

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
    if options.precompute_lor || options.implementation == 3
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
    if options.precompute_lor || options.implementation == 3 || options.implementation == 2
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

x_mramla = zeros(options.Nx,options.Ny,options.Nz, options.Niter + 1);
x_mramla(:,:,:,1) = options.x0;
x_mramla = reshape(x_mramla, options.Nx*options.Ny*options.Nz, options.Niter + 1);

Amin = zeros(length(Sino),1);
pj = zeros(options.Nx*options.Ny*options.Nz,options.subsets);
for ll = 1 : options.subsets
    [A] = observation_matrix_formation_nongate(options, ll, index, LL, pituus, lor);
    pj(:,ll) = A'*ones(size(A,1),1,'double');
    [m,~] = size(A);
    rowMin = nan(m, 1);
    [I,~,S] = find(A);
    [I,K] = sort(I);
    S = S(K);
    markers = [find([1; diff(I)]); numel(I)+1];
    iRows = I(markers(1:end-1));
    for i = 1:numel(iRows)
        s = S(markers(i):(markers(i+1)-1));
        rowMin(iRows(i)) = min(s);
    end
    rowMin(isnan(rowMin)) = 0;
    if options.precompute_lor == false
        Amin(index{ll}) = rowMin;
    else
        Amin(pituus(ll)+1:pituus(ll + 1)) = rowMin;
    end
    clear I K S markers rowMin s iRows
end
D = sum(pj,2);

lam = zeros(options.Niter,1);
lam(1) = options.b0;
for i=1:options.Niter
    lam(i+1) = 0.5*lam(i);
end

if options.U == 0
    U = max(double(Sino)./Amin);
else
    U = options.U;
end
pj3 = D/options.subsets;
dU = zeros(size(x_mramla,1),1);

for ii = 1 : options.Niter
    rmM_apu = x_mramla(:,ii);
    for kk = 1 : options.subsets
        [A] = observation_matrix_formation_nongate(options, kk, index, LL, pituus, lor);
        
        if options.precompute_lor == false
            uu = double(Sino(index{kk}));
        else
            uu = double(Sino(pituus2(kk)+1:pituus2(kk + 1)));
        end
        tStart = tic;
        pp = rmM_apu<U/2;
        dU(pp) = rmM_apu(pp)./(pj3(pp)+epps);
        dU(~pp) = (U-rmM_apu(~pp))./(pj3(~pp)+epps);
        rmM_apu = rmM_apu + lam(ii).*dU.*(A'*(uu./(A*rmM_apu+epps))-full(sum(A,1)'));
        rmM_apu(rmM_apu<0) = 0;
        tElapsed = toc(tStart);
        disp(['MRAMLA sub-iteration ' num2str(kk) ' took ' num2str(tElapsed) ' seconds'])
        disp(['MRAMLA sub-iteration ' num2str(kk) ' finished'])
    end
    x_mramla(:,ii+1) = rmM_apu;
end
x_mramla = reshape(x_mramla,options.Nx,options.Ny,options.Nz, options.Niter + 1);

end