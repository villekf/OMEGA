function pz = reconstructions_main(options)
%% Main reconstruction file
% This function is used to compute various reconstructions with the
% selected method. Can be used with any sinogram or raw data.

if options.use_raw_data
    if isfield(options, 'coincidences') == 0
        if options.partitions == 1
            if options.use_ASCII
                load([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'coincidences')
            elseif options.use_LMF
                load([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'coincidences')
            elseif options.use_root
                load([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'coincidences')
            end
        else
            if options.use_ASCII
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_ASCII.mat'], 'coincidences')
            elseif options.use_LMF
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_LMF.mat'], 'coincidences')
            elseif options.use_root
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_root.mat'], 'coincidences')
            end
        end
        SinM = coincidences;
    else
        SinM = options.coincidences;
    end
    
    clear coincidences options.coincidences
else
    if options.partitions == 1 && isfield(options, 'SinM') == 0
        load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinM')
    elseif isfield(options, 'SinM') == 0
        load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'], 'SinM')
    else
        SinM = options.SinM;
        clear options.SinM
    end
end

if options.precompute_lor == false && options.reconstruction_method == 3
    error('precompute_lor must be set to true if using method 3')
end

rekot = false(15,1);
if options.mlem
    rekot(1) = true;
end
if options.osem
    rekot(2) = true;
end
if options.mramla
    rekot(3) = true;
end
if options.ramla
    rekot(4) = true;
end
if options.ecosem
    rekot(5) = true;
end
if options.cosem
    rekot(6) = true;
end
if options.acosem
    rekot(15) = true;
end
if options.mrp_osl
    rekot(7) = true;
end
if options.mrp_bsrem
    rekot(8) = true;
end
if options.quad_osl
    rekot(9) = true;
end
if options.quad_bsrem
    rekot(10) = true;
end
if options.L_osl
    rekot(11) = true;
end
if options.L_bsrem
    rekot(12) = true;
end
if options.FMH_osl
    rekot(13) = true;
end
if options.weighted_mean_osl
    rekot(14) = true;
end

Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
Niter = options.Niter;
subsets = options.subsets;
epps = options.epps;
precompute_obs_matrix = options.precompute_obs_matrix;
attenuation_correction = options.attenuation_correction;
diameter = options.diameter;
FOVa = options.FOVa;
axial_fov = options.axial_fov;
NSinos = options.NSinos;
pseudot = int32(options.pseudot);
rings = options.rings;
det_per_ring = options.det_per_ring;
machine_name = options.machine_name;
Nang = options.Nang;
Ndist = options.Ndist;
cr_pz = options.cr_pz;
use_raw_data = options.use_raw_data;
TotSinos = options.TotSinos;
attenuation_datafile = options.attenuation_datafile;
partitions = options.partitions;
verbose = options.verbose;
device = int32(options.use_device);

pz = cell(length(rekot) + 1,partitions);

N = Nx * Ny * Nz;


if options.reconstruction_method == 3 || options.reconstruction_method == 4
    if rekot(1) && rekot(2)
        if subsets == 1
            disp(['Both OSEM and MLEM set for method ' num2str(options.reconstruction_method) ', using MLEM'])
            rekot(2) = false;
        else
            disp(['Both OSEM and MLEM set for method ' num2str(options.reconstruction_method) ', using OSEM'])
            rekot(1) = false;
        end
    end
end

if options.reconstruction_method == 3
    if rekot(1)
        pz_ml = ones(N,Niter + 1);
        pz_ml(:,1) = options.x0(:);
    end
    if rekot(2)
        pz_osem = ones(N,Niter + 1);
        pz_osem(:,1) = options.x0(:);
    end
end

if options.reconstruction_method == 1
    
    Ndx = options.Ndx;
    Ndy = options.Ndy;
    Ndz = options.Ndz;
    a_L = options.a_L;
    
    weights = options.weights;
    fmh_weights = options.fmh_weights;
    weighted_weights = options.weighted_weights;
    
    
    
    if rekot(2)
        pz_osem = ones(N,Niter + 1);
        pz_osem(:,1) = options.x0(:);
        pz_osem_apu = pz_osem(:,1);
    end
    
    if rekot(3)
        pz_rmM = ones(N,Niter + 1);
        pz_rmM(:,1) = options.x0(:);
        pz_rmM_apu = pz_rmM(:,1);
    end
    
    if rekot(4)
        pz_rm = ones(N,Niter + 1);
        pz_rm(:,1) = options.x0(:);
        pz_rm_apu = pz_rm(:,1);
    end
    
    if rekot(5)
        pz_eco = ones(N,Niter + 1);
        pz_eco(:,1) = options.x0(:);
        pz_eco_apu = pz_eco(:,1);
    end
    
    if rekot(6)
        pz_cosem = ones(N,Niter + 1);
        pz_cosem(:,1) = options.x0(:);
        pz_cosem_apu = pz_cosem(:,1);
    end
    
    if rekot(15)
        pz_acosem = ones(N,Niter + 1);
        pz_acosem(:,1) = options.x0(:);
        pz_acosem_apu = pz_acosem(:,1);
    end
    
    if rekot(7)
        pz_mrp_osl = ones(N,Niter + 1);
        pz_mrp_osl(:,1) = options.x0(:);
        pz_mrp_osl_apu = pz_mrp_osl(:,1);
    end
    
    if rekot(8)
        pz_mrp_bsrem = ones(N,Niter + 1);
        pz_mrp_bsrem(:,1) = options.x0(:);
        pz_mrp_bsrem_apu = pz_mrp_bsrem(:,1);
    end
    
    if rekot(9)
        pz_quad_osl = ones(N,Niter + 1);
        pz_quad_osl(:,1) = options.x0(:);
        pz_quad_osl_apu = pz_quad_osl(:,1);
        if ~isempty(weights)
            if length(weights(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            elseif length(weights(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            end
            if ~isinf(weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))))
                weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))) = Inf;
            end
        end
    end
    
    if rekot(10)
        pz_quad_bsrem = ones(N,Niter + 1);
        pz_quad_bsrem(:,1) = options.x0(:);
        pz_quad_bsrem_apu = pz_quad_bsrem(:,1);
        if ~isempty(weights)
            if length(weights(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            elseif length(weights(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            end
            if ~isinf(weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))))
                weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))) = Inf;
            end
        end
    end
    
    if rekot(11)
        pz_L_osl = ones(N,Niter + 1);
        pz_L_osl(:,1) = options.x0(:);
        pz_L_osl_apu = pz_L_osl(:,1);
        if ~isempty(a_L)
            if length(a_L(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector a_L is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            elseif length(a_L(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector a_L is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            end
        end
    end
    
    if rekot(12)
        pz_L_bsrem = ones(N,Niter + 1);
        pz_L_bsrem(:,1) = options.x0(:);
        pz_L_bsrem_apu = pz_L_bsrem(:,1);
        if ~isempty(a_L)
            if length(a_L(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector a_L is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            elseif length(a_L(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector a_L is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            end
        end
    end
    
    if rekot(13)
        pz_fmh_osl = ones(N,Niter + 1);
        pz_fmh_osl(:,1) = options.x0(:);
        pz_fmh_osl_apu = pz_fmh_osl(:,1);
        if ~isempty(fmh_weights)
            if length(fmh_weights(:)) < floor(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))/2)*3
                error(['Weights vector fmh_weights is too small, needs to be [3, ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) '] in size'])
            elseif length(fmh_weights(:)) > floor(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))/2)*3
                error(['Weights vector fmh_weights is too large, needs to be [3, ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) '] in size'])
            end
        end
    end
    
    if rekot(14)
        pz_weighted_osl = ones(N,Niter + 1);
        pz_weighted_osl(:,1) = options.x0(:);
        pz_weighted_osl_apu = pz_weighted_osl(:,1);
        if ~isempty(weighted_weights)
            if length(weighted_weights(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector weighted_weights is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            elseif length(weighted_weights(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
                error(['Weights vector weighted_weights is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
            end
            if ~isinf(weighted_weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))))
                weighted_weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))) = Inf;
            end
        end
    end
    
end

%%
for llo = 1 : partitions
    if iscell(SinM)
        Sino = SinM{llo};
    else
        Sino = SinM;
        clear SinM
    end
    
    Sino = Sino(:);
    
    if issparse(Sino)
        Sino = (full(Sino));
    end
    
    if use_raw_data == false && NSinos ~= TotSinos
        Sino = Sino(1:NSinos*Nang*Ndist);
    end
    
    %% This computes a whole observation matrix and uses it to compute the MLEM (no on-the-fly calculations)
    if precompute_obs_matrix && options.reconstruction_method == 1
        
        
        if rekot(1)
            pituus = int32(length(Sino));
            A = observation_matrix_formation(diameter, FOVa, axial_fov, rings, pseudot, Nx, Ny, Nz, det_per_ring, cr_pz, options.use_fsparse, attenuation_correction, attenuation_datafile, options.precompute_lor, options.use_raw_data, pituus);
            
            ll = ll + 1;
            
            pz_ml = ones(N,Niter);
            pz_ml(:,1) = options.x0(:);
            D = A*ones(size(A,1),1,'double');
            Sino = double(Sino);
            
            for iter=1:Niter
                if rekot(1)
                    pz_ml(:,iter+1)=(pz_ml(:,iter)./(D+epps)).*((A*(Sino./(A'*pz_ml(:,iter)+epps))));
                    disp(['MLEM iteration ' num2str(iter) ' finished'])
                end
                if sum(rekot(2:end)) > 0
                    disp('Only MLEM is supported with precomputed observation matrix')
                end
            end
            if Nz == 1
                pz_ml = reshape(pz_ml,Nx,Ny,Niter+1);
            else
                pz_ml = reshape(pz_ml,Nx,Ny,Nz,Niter+1);
            end
        else
            error('Only MLEM is supported with precomputed observation matrix')
        end
        
        pz{1, llo} = pz_ml;
        
    else
        %% This part is used when the observation matrix is calculated on-the-fly
        
        % Compute the indices for the subsets used
        % When using sinogram data, collect the indices such that the
        % sinograms are divided evenly, if possible
        % Sinogram angles are used as the dimension where the subsets are
        % formed
        if use_raw_data == false && subsets > 1
            if options.precompute_lor || options.reconstruction_method == 3 || options.reconstruction_method == 5
                load([machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor','discard')
                if length(discard) ~= TotSinos*Nang*Ndist
                    error('Error: Size mismatch between sinogram and LORs to be removed')
                end
                if use_raw_data == false && NSinos ~= TotSinos
                    discard = discard(1:NSinos*Nang*Ndist);
                end
                lor_a = lor;
                clear lor
                ind_apu = uint32(find(discard));
                port = ceil((Nang-subsets+1)/subsets);
                over = Nang - port*subsets;
                index = cell(subsets,1);
                pituus = zeros(subsets, 1, 'uint32');
                for i=1:subsets
                    if over>0
                        index1 = uint32(sort(sub2ind([Nang Ndist NSinos],repmat(repelem(i:subsets:(port + 1)*subsets,Ndist)',NSinos,1),repmat((1:Ndist)',(port+1)*NSinos,1),repelem((1:NSinos)',Ndist*(port+1),1))));
                        over = over - 1;
                    else
                        index1 = uint32(sort(sub2ind([Nang Ndist NSinos],repmat(repelem(i:subsets:port*subsets,Ndist)',NSinos,1),repmat((1:Ndist)',port*NSinos,1),repelem((1:NSinos)',Ndist*port,1))));
                    end
                    index{i} = index1(ismember(index1, ind_apu));
                    pituus(i) = int32(length(index{i}));
                end
                index = cell2mat(index);
                index = index(ismember(index, ind_apu));
                clear index1 ind_apu
            else
                port = ceil((Nang-subsets+1)/subsets);
                over = Nang - port*subsets;
                index = cell(subsets,1);
                pituus = zeros(subsets, 1, 'uint32');
                for i=1:subsets
                    if over>0
                        index1 = uint32(sort(sub2ind([Nang Ndist NSinos],repmat(repelem(i:subsets:(port + 1)*subsets,Ndist)',NSinos,1),repmat((1:Ndist)',(port+1)*NSinos,1),repelem((1:NSinos)',Ndist*(port+1),1))));
                        over = over - 1;
                    else
                        index1 = uint32(sort(sub2ind([Nang Ndist NSinos],repmat(repelem(i:subsets:port*subsets,Ndist)',NSinos,1),repmat((1:Ndist)',port*NSinos,1),repelem((1:NSinos)',Ndist*port,1))));
                    end
                    index{i} = uint32(index1);
                    pituus(i) = int32(length(index1));
                end
                clear index1
            end
        elseif subsets > 1
            % for raw list-mode data, take the subsets randomly
            % last subset has all the spare indices
            if options.precompute_lor || options.reconstruction_method == 3 || options.reconstruction_method == 2 || options.reconstruction_method == 5
                load([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'LL','lor')
                lor_a = lor;
                clear lor
                indices = uint32(length(LL));
                index = cell(subsets, 1);
                port = uint32(floor(length(LL)/subsets));
                if options.use_Shuffle
                    apu = Shuffle(indices(end), 'index')';
                else
                    apu = uint32(randperm(indices(end)))';
                end
                pituus = zeros(subsets, 1, 'uint32');
                for i = 1 : subsets
                    if i == subsets
                        index{i} = apu(port*(i-1)+1:end);
                    else
                        index{i} = apu(port*(i-1)+1:(port*(i)));
                    end
                    pituus(i) = int32(length(index{i}));
                end
                clear apu
            else
                indices = uint32(length(Sino));
                index = cell(subsets, 1);
                port = uint32(floor(length(Sino)/subsets));
                if options.use_Shuffle
                    apu = Shuffle(indices(end), 'index')';
                else
                    apu = uint32(randperm(indices(end)))';
                end
                for i = 1 : subsets
                    if i == subsets
                        index{i} = apu(port*(i-1)+1:end);
                    else
                        index{i} = apu(port*(i-1)+1:(port*(i)));
                    end
                end
                clear apu
            end
        end
        
        
        
        %%
        
        % Diameter of the PET-device (bore) (mm)
        R=double(diameter);
        % Transaxial FOV (mm)
        FOVa=double(FOVa);
        % Axial FOV (mm)
        axial_fow = double(axial_fov);
        % Number of rings
        blocks=int32(rings + length(pseudot) - 1);
        % First ring
        block1=int32(0);
        % Pixel count in x- and y-directions
        pikselikoko=int32(Nx);
        
        NSinos = int32(NSinos);
        NSlices = int32(Nz);
        TotSinos = int32(TotSinos);
        
        if use_raw_data == false
            load([machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
            load([machine_name '_3D_coordinates_span' num2str(options.span) '_ringdiff' num2str(options.ring_difference) '_' num2str(TotSinos) '.mat'],'z')
            if NSinos ~= TotSinos
                z = z(1:NSinos,:);
            end
        else
            load([machine_name '_detector_coordinates.mat'],'x','y');
            if ~exist('LL','var')
                load([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'LL')
            end
            
            z_length = double(blocks) * cr_pz;
            z = linspace(0, z_length, blocks + 1);
            
            
        end
        if attenuation_correction
            data = load(attenuation_datafile);
            variables = fields(data);
            vaimennus = double(data.(variables{1}));
            if size(vaimennus,1) ~= Nx || size(vaimennus,2) ~= Ny || size(vaimennus,3) ~= Nz
                if size(vaimennus,1) ~= Nx*Ny*Nz
                    error("Error: Attenuation data is of different size than the reconstructed image")
                end
            end
            if size(vaimennus,2) == 1
                vaimennus = vaimennus(:,:,2*block1+1:2*blocks+1);
            else
                vaimennus = vaimennus(2*block1+1:(2*blocks+1)*Nx*Ny);
            end
            vaimennus = vaimennus(:);
            if options.reconstruction_method > 2 && options.reconstruction_method ~= 4
                vaimennus = single(vaimennus);
            end
            clear data
        else
            if options.reconstruction_method > 2 && options.reconstruction_method ~= 4
                vaimennus = single(0);
            else
                vaimennus = 0;
            end
        end
        
        if min(min(z)) == 0
            z = z + (axial_fow - max(max(z)))/2;
        end
        
        if options.reconstruction_method == 2 || options.reconstruction_method == 4
            x=single(x);
            y=single(y);
            z_det = single(z);
        else
            x=double(x);
            y=double(y);
            z_det = double(z);
        end
        
        size_x = int32(size(x,1));
        
        if options.precompute_lor && subsets > 1 || options.reconstruction_method == 2  && subsets > 1 || options.reconstruction_method == 4 && subsets > 1
            pituus = [0;cumsum(pituus)];
            if iscell(index)
                index = cell2mat(index);
            end
        end
        
        if use_raw_data
            if isempty(pseudot)
                pseudot = int32(1e5);
            else
                pseudot = pseudot - 1;
            end
        end
        
        % for the precomputed version, index vectors are needed
        if use_raw_data == false && options.precompute_lor || use_raw_data == false && options.reconstruction_method == 2 || options.reconstruction_method == 4 && use_raw_data == false
            
            
            
            if subsets > 1
                lor_a = (lor_a(index));
                Sino = Sino(index);
            else
                load([machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor','discard')
                if length(discard) ~= TotSinos*Nang*Ndist
                    error('Error: Size mismatch between sinogram and LORs to be removed')
                end
                if use_raw_data == false && NSinos ~= TotSinos
                    discard = discard(1:NSinos*Nang*Ndist);
                end
                lor_a = (lor(discard));
                Sino = Sino(discard);
                clear lor
            end
            [~, I] = sort(y, 2);
            sy = size(y);
            I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
            
            xy_index = uint32(I(:,1));
            xy_index2 = repmat(uint32(1:size_x)', NSinos - Nz, 1);
            xy_index = [repmat(xy_index, Nz, 1); xy_index2];
            [~, I] = sort(z_det, 2);
            sy = size(z_det);
            I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
            
            z_index = uint32(I(:,1));
            z_index = repelem(z_index, size_x);
            if subsets > 1
                z_index = z_index(index);
            else
                z_index = (z_index(discard));
            end
            apu = z_index > NSinos;
            z_index = z_index - 1;
            
            if subsets > 1
                xy_index = xy_index(index);
            else
                xy_index = (xy_index(discard));
            end
            xy_index(apu) = xy_index(apu) + uint32(size_x);
            xy_index = xy_index - 1;
            
            summa = zeros(subsets, 1, 'uint64');
            
            if subsets > 1
                for kk = 1 : subsets
                    summa(kk) = sum(int64(lor_a(pituus(kk)+1:pituus(kk+1))));
                    koko(kk) = -(pituus(kk))+pituus(kk+1);
                end
                vector_size = max(koko);
            else
                summa = uint64(sum(int64(lor_a)));
                pituus = int32([0;length(Sino)]);
            end
            
            
            clear discard I yt xt xy_index2 index apu
        elseif use_raw_data && options.precompute_lor || use_raw_data && options.reconstruction_method == 2  || use_raw_data && options.reconstruction_method == 4
            
            if subsets > 1
                LL = LL(index,:);
                lor_a = (lor_a(index));
                Sino = Sino(index);
            else
                load([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'lor')
                pituus = int32([0;length(Sino)]);
                lor_a = lor;
                clear lor
            end
            summa = zeros(subsets, 1, 'uint64');
%             y = circshift(y, -length(y)/4);
%             x = circshift(x, -length(x)/4);
            
            for kk = 1 : subsets
                apu = LL(pituus(kk) + 1 : pituus(kk + 1),:) - 1;
                apu2 = idivide(apu, uint16(det_per_ring));
                idx = apu2(:,1) == apu2(:,2);
                apu2 = apu(idx,:);
                ind = mod(apu2, uint16(det_per_ring)) + 1;
                yt = y(ind);
                y_i = yt(:,1) > yt(:,2);
                apu2(y_i,:) = fliplr(apu2(y_i,:));
                apu(idx,:) = apu2;
                LL(pituus(kk) + 1 : pituus(kk + 1),:) = apu + 1;
                summa(kk) = sum(int64(lor_a(pituus(kk)+1:pituus(kk+1))));
                koko(kk) = -(pituus(kk))+pituus(kk+1);
            end
            
            vector_size = max(koko);
            clear apu apu2 idx ind yt y_i index discard
            
            LL = LL';
            LL = LL(:);
        end
        
        
        % Pixels
        etaisyys=(R-FOVa)/2;
        if options.reconstruction_method == 2 || options.reconstruction_method == 4
            zz=linspace(single(0),single(axial_fow),Nz+1);
            xx = single(linspace(etaisyys,R-etaisyys,pikselikoko+1));
            yy = single(linspace(etaisyys,R-etaisyys,pikselikoko+1));
        else
            zz=linspace(double(0),double(axial_fow),Nz+1);
            xx = double(linspace(etaisyys,R-etaisyys,pikselikoko+1));
            yy = double(linspace(etaisyys,R-etaisyys,pikselikoko+1));
        end
        zz=zz(2*block1+1:2*blocks+2);
        
        % Distance of adjacent pixels
        d=diff(xx(1:2));
        dz=diff(zz(1:2));
        
        % Distance of image from the origin
        bx=xx(1);
        by=yy(1);
        bz=zz(1);
        
        % Number of pixels
        Ny=int32(Ny);
        Nx=int32(Nx);
        Nz=int32(Nz);
        
        if options.reconstruction_method == 2 || options.reconstruction_method == 4
            iij=single(0:Nx);
            jji=single(0:Ny);
            kkj=single(0:Nz);
        else
            iij=double(0:Nx);
            jji=double(0:Ny);
            kkj=double(0:Nz);
        end
        
        N=(Nx)*(Ny)*(Nz);
        det_per_ring = int32(det_per_ring);
        
        % How much memory is preallocated
        if use_raw_data == false
            ind_size = int32(NSinos/subsets*(det_per_ring)* Nx * (Ny));
        else
            ind_size = int32((det_per_ring)^2/subsets* Nx * (Ny));
        end
        
        
        zmax = max(max(z_det));
        if zmax==0
            if options.reconstruction_method == 2 || options.reconstruction_method == 4
                zmax = single(1);
            else
                zmax = double(1);
            end
        end
        %%
        if options.reconstruction_method == 1
            
            if rekot(3) || rekot(5) || rekot(6) || rekot(9) || rekot(10)  || rekot(8)  || rekot(11) || rekot(12) || rekot(15)
                pj = zeros(N,subsets);
                if rekot(15)
                    C_aco = zeros(double(N), subsets);
                end
                if rekot(6) || rekot(5)
                    C_co = zeros(double(N), subsets);
                end
                if rekot(3)
                    Amin = zeros(length(Sino),1);
                end
                
                if rekot(6) || rekot(5) || rekot(3)
                    
                    if verbose
                        disp('Prepass phase for MRAMLA, COSEM and ECOSEM started')
                    end
                    for osa_iter = 1 : subsets
                        if options.precompute_lor == false
                            if use_raw_data == false
                                [ lor, indices, alkiot] = improved_Siddon_algorithm( options.verbose, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, ...
                                    xx , NSinos, NSlices, size_x, zmax, NSinos, ind_size, vaimennus, index{osa_iter}, pituus(osa_iter), attenuation_correction);
                            else
                                L = LL(index{osa_iter},:);
                                L = L';
                                L = L(:);
                                [ lor, indices, alkiot] = improved_Siddon_algorithm_raw( options.verbose, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, ...
                                    yy, xx , NSinos, NSlices, size_x, zmax, NSinos, ind_size, vaimennus, L, pseudot, block1, blocks, det_per_ring, attenuation_correction);
                            end
                            lor = reshape(lor,[],2);
                            lor=repelem(int32((lor(:,1))),lor(:,2));
                            A_length = length(Sino(index{osa_iter}));
                            indices=indices + 1;
                            if verbose
                                tStart = tic;
                            end
                            if use_raw_data == false
                                if options.use_fsparse == false
                                    A = sparse(double(lor),double(indices),double(alkiot),A_length,double(N));
                                else
                                    A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
                                end
                            else
                                if options.use_fsparse == false
                                    A = sparse(double(lor),double(indices),double(alkiot),A_length,double(N));
                                else
                                    A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
                                end
                            end
                            clear mex indices alkiot lor
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
                            end
                            A = A';
                        else
                            lor2 = [0; cumsum(uint32(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                            if use_raw_data == false
                                [A] = improved_Siddon_algorithm_array( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, ...
                                    size_x, zmax, vaimennus, lor2, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, ...
                                    lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), summa(osa_iter), xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                                    z_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), NSinos, options.verbose);
                            else
                                [A] = improved_Siddon_algorithm_array_raw( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , size_x, zmax, ...
                                    vaimennus, lor2, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                                    summa(osa_iter), LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2), pseudot, det_per_ring, options.verbose);
                            end
                            clear lor2
                        end
                        
                        pj(:,osa_iter) = A*ones(size(A,2),1,'double');
                        if rekot(5) || rekot(6) || rekot(15)
                            if options.precompute_lor == false
                                uu = double(Sino(index{osa_iter}));
                            else
                                uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                            end
                        end
                        if rekot(15)
                            C_aco(:,osa_iter) = full(sum(spdiags(uu./(A'*options.x0(:)+epps),0,size(A,2),size(A,2))*(bsxfun(@times,A',(options.x0(:)')).^(1/options.h)))');
                        end
                        if rekot(6) || rekot(5)
                            C_co(:,osa_iter) = full(sum(spdiags(uu./(A'*options.x0(:)+epps),0,size(A,2),size(A,2))*(bsxfun(@times,A',(options.x0(:)'))))');
                        end
                        if rekot(3)
                            %%%% This particular piece of code was taken from:
                            %%%% https://se.mathworks.com/matlabcentral/answers/35309-max-min-of-sparse-matrices
                            [~,m] = size(A);
                            rowMin = nan(m, 1);
                            [~,I,S] = find(A);
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
                                Amin(index{osa_iter}) = rowMin;
                            else
                                Amin(pituus(osa_iter)+1:pituus(osa_iter + 1)) = rowMin;
                            end
                            clear I K S markers rowMin s iRows
                        end
                        clear A
                    end
                    D = sum(pj,2);
                    if verbose
                        disp('Prepass phase for MRAMLA, COSEM, ACOSEM and ECOSEM completed')
                    end
                end
                if rekot(3) || rekot(4) || rekot(8) || rekot(10)
                    lam = zeros(Niter,1);
                    lam(1) = options.b0;
                    for i=1:Niter
                        lam(i+1) = 0.5*lam(i);
                    end
                end
                if rekot(3)
                    if options.U == 0
                        U = max(double(Sino)./Amin);
                    else
                        U = options.U;
                    end
                    pj3 = D/subsets;
                    dU = zeros(size(pz_rmM,1),1);
                end
                if rekot(9) || rekot(10) || rekot(11) || rekot(12)
                    distX = FOVa/double(Nx);
                    distZ = (double(axial_fov)/double(Nz));
                    if isempty(options.weights)
                        for jj = Ndz : -1 : -Ndz
                            cc = [];
                            for kk = Ndy : -1 : -Ndy
                                if Ndz == 0 || Nz == 1
                                    apu = [((Ndx:-1:-Ndx) * distX)', (repelem(kk,Ndx*2+1) * distX)'];
                                else
                                    apu = [((Ndx:-1:-Ndx) * distX)', (repelem(kk,Ndx*2+1) * distX)', (repelem(jj,Ndx*2+1) * distZ)'];
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
                    pz_pad = padding(reshape(options.x0(:),Nx,Ny,Nz),[Ndx Ndy Ndz]);
                    s = size(pz_pad);
                    N_pad = min(3, Ndx + Ndy + Ndz);
                    [c1{1:N_pad}]=ndgrid(1:(Ndx*2+1));
                    c2(1:N_pad)={Ndy+1};
                    if Ndz > Ndx && Ndz > 1
                        apu = c1{1};
                        apu2 = c1{2};
                        apu3 = c1{3};
                        for kk = 1 : Ndz
                            apu(:,:,end+1) = apu(:,:,end);
                            apu2(:,:,end+1) = apu2(:,:,end);
                            apu3(:,:,end+1) = apu3(:,:,end) + 1;
                        end
                        c1{1} = apu;
                        c1{2} = apu2;
                        c1{3} = apu3;
                        c2(end) = {Ndz+1};
                    elseif Ndz < Ndx
                        apu = c1{1};
                        apu2 = c1{2};
                        apu3 = c1{3};
                        for kk = 1 : 2*(Ndx-Ndz)
                            apu(:,:,end) = [];
                            apu2(:,:,end) = [];
                            apu3(:,:,end) = [];
                        end
                        c1{1} = apu;
                        c1{2} = apu2;
                        c1{3} = apu3;
                        c2(end) = {Ndz+1};
                    end
                    offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
                    if Nz == 1
                        tr_ind = sub2ind([Nx+Ndx*2 Ny+Ndy*2],mod((1:N)'-1,Nx)+(Ndx + 1),mod(floor(((1:double(N))'-1)/double(Nx)),double(Ny))+(Ndy + 1));
                    else
                        tr_ind = sub2ind([Nx+Ndx*2 Ny+Ndy*2 Nz+Ndz*2],mod((1:N)'-1,Nx)+(Ndx + 1),mod(floor(((1:double(N))'-1)/double(Nx)),double(Ny))+(Ndy + 1),floor(((1:double(N))'-1)/double(Nx*Ny))+(Ndz+1));
                    end
                    tr_offsets = bsxfun(@plus,tr_ind,offsets(:)');
                    if verbose
                        disp('Prepass phase for quadratic prior and L-filter completed')
                    end
                end
                if rekot(9)
                    pz_pad_osl = pz_pad;
                end
                if rekot(9) || rekot(10)
                    if isempty(options.weights)
                        weights = weights/sum(weights(~isinf(weights)));
                        weights_quad = [weights(1:floor(length(weights)/2));weights(ceil(length(weights)/2) + 1:end)];
                    end
                end
                if rekot(11)
                    pz_pad_L_osl  = pz_pad;
                end
                if rekot(11) || rekot(12)
                    if isempty(a_L)
                        dd = (1:find(isinf(weights)))';
                        dd = dd/sum(dd);
                        a_L = [dd;flip(dd(1:end-1))];
                    end
                end
                if rekot(13)
                    if isempty(options.fmh_weights)
                        for jjj = 1:find(isinf(weights))-1
                            fmh_weights = [fmh_weights weights(jjj:find(isinf(weights))-jjj:end + 1 -jjj)];
                        end
                        fmh_weights(isinf(fmh_weights)) = max(max(fmh_weights(~isinf(fmh_weights))))*options.fmh_center_weight;
                        fmh_weights = fmh_weights./sum(fmh_weights,1);
                    end
                    pz_pad_fmh = pz_pad;
                end
                if rekot(14)
                    if isempty(options.weighted_weights)
                        weighted_weights = weights;
                        weighted_weights(isinf(weighted_weights)) = max(weighted_weights(~isinf(weighted_weights)))*options.weighted_center_weight;
                    end
                    pz_pad_weighted = pz_pad;
                end
            end
        end
        
        %%
        
        if options.reconstruction_method == 1
            for iter=1:Niter
                if rekot(2)
                    pz_osem_apu = pz_osem(:,iter);
                end
                if rekot(3)
                    pz_rmM_apu = pz_rmM(:,iter);
                end
                if rekot(4)
                    pz_rm_apu = pz_rm(:,iter);
                end
                if rekot(5)
                    pz_eco_apu = pz_eco(:,iter);
                end
                if rekot(6)
                    pz_cosem_apu = pz_cosem(:,iter);
                end
                if rekot(15)
                    pz_acosem_apu = pz_acosem(:,iter);
                end
                if rekot(7)
                    pz_mrp_osl_apu = pz_mrp_osl(:,iter);
                end
                if rekot(8)
                    pz_mrp_bsrem_apu = pz_mrp_bsrem(:,iter);
                end
                if rekot(9)
                    pz_quad_osl_apu = pz_quad_osl(:,iter);
                end
                if rekot(10)
                    pz_quad_bsrem_apu = pz_quad_bsrem(:,iter);
                end
                if rekot(11)
                    pz_L_osl_apu = pz_L_osl(:,iter);
                end
                if rekot(12)
                    pz_L_bsrem_apu = pz_L_bsrem(:,iter);
                end
                if rekot(13)
                    pz_fmh_osl_apu = pz_fmh_osl(:,iter);
                end
                if rekot(14)
                    pz_weighted_osl_apu = pz_weighted_osl(:,iter);
                end
                
                for osa_iter = 1 : subsets
                    if options.precompute_lor == false
                        if use_raw_data == false
                            [ lor, indices, alkiot] = improved_Siddon_algorithm( options.verbose, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx ,...
                                NSinos, NSlices, size_x, zmax, NSinos, ind_size, vaimennus, index{osa_iter}, pituus(osa_iter), attenuation_correction);
                        else
                            L = LL(index{osa_iter},:);
                            L = L';
                            L = L(:);
                            [ lor, indices, alkiot] = improved_Siddon_algorithm_raw( options.verbose, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx ,...
                                NSinos, NSlices, size_x, zmax, NSinos, ind_size, vaimennus, L, pseudot, block1, blocks, det_per_ring, attenuation_correction);
                        end
                        lor = reshape(lor,[],2);
                        lor=repelem(int32((lor(:,1))),lor(:,2));
                        uu = double(Sino(index{osa_iter}));
                        
                        A_length = length(uu);
                        indices=indices + 1;
                        if verbose
                            tStart = tic;
                        end
                        if options.use_fsparse == false
                            A = sparse(double(lor),double(indices),double(alkiot), A_length, double(N));
                        else
                            A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
                        end
                        clear indices alkiot lor
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
                        end
                        A = A';
                    else
                        lor2 = [0; cumsum(uint32(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                        if use_raw_data == false
                            [A] = improved_Siddon_algorithm_array( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, size_x, zmax, vaimennus, ...
                                lor2, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), summa(osa_iter), ...
                                xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), z_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), NSinos, options.verbose);
                        else
                            [A] = improved_Siddon_algorithm_array_raw( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , size_x, zmax, vaimennus, lor2, ...
                                pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), summa(osa_iter), ...
                                LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2), pseudot, det_per_ring, options.verbose);
                        end
                        uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                        clear lor2
                    end
                    
                    Summ = full(sum(A,2));
                    if rekot(2)
                        if verbose
                            tStart = tic;
                        end
                        pz_osem_apu=(pz_osem_apu./(Summ+epps)).*((A*(uu./(A'*pz_osem_apu+epps))));
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['OSEM sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(3)
                        if verbose
                            tStart = tic;
                        end
                        pp = pz_rmM_apu<U/2;
                        dU(pp) = pz_rmM_apu(pp)./(pj3(pp)+epps);
                        dU(~pp) = (U-pz_rmM_apu(~pp))./(pj3(~pp)+epps);
                        pz_rmM_apu = pz_rmM_apu + lam(iter).*dU.*(A*(uu./(A'*pz_rmM_apu+epps))-Summ);
                        pz_rmM_apu(pz_rmM_apu<0) = 0;
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(6)
                        if verbose
                            tStart = tic;
                        end
                        C_co(:,osa_iter) = full(sum(spdiags(uu./(A'*pz_cosem_apu+epps),0,size(A,2),size(A,2))*(bsxfun(@times,A',pz_cosem_apu')))');
                        apu = (sum(C_co,2)./D);
                        pz_cosem_apu = (apu)*sum(uu)/sum(A'*apu+epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['COSEM sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(15)
                        if verbose
                            tStart = tic;
                        end
                        C_aco(:,osa_iter) = full(sum(spdiags(uu./(A'*pz_acosem_apu+epps),0,size(A,2),size(A,2))*(bsxfun(@times,A',pz_acosem_apu').^(1/options.h)))');
                        apu = (sum(C_aco,2)./D).^options.h;
                        pz_acosem_apu = (apu)*sum(uu)/sum(A'*apu+epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(5)
                        if verbose
                            tStart = tic;
                        end
                        pz_eco_apuw = pz_eco_apu;
                        alpha_eco = 1;
                        if rekot(2) == false
                            pz_osem_apu=(pz_osem_apu./(Summ+epps)).*((A*(uu./(A'*pz_osem_apu+epps))));
                        end
                        if rekot(6) == false
                            C_co(:,osa_iter) = full(sum(spdiags(uu./(A'*pz_cosem_apu+epps),0,size(A,2),size(A,2))*(bsxfun(@times,A',pz_cosem_apu')))');
                            apu = (sum(C_co,2)./D);
                            pz_cosem_apu = (apu)*sum(uu)/sum(A'*apu+epps);
                        end
                        pz_eco_apu = alpha_eco*pz_osem_apu+(1-alpha_eco)*pz_cosem_apu;
                        while (alpha_eco>0.0096 && (sum(D.*(-pz_cosem_apu.*log(pz_eco_apuw+epps)+pz_eco_apuw))<sum(D.*(-pz_cosem_apu.*log(pz_eco_apu+epps)+pz_eco_apu))))
                            alpha_eco = alpha_eco*0.9;
                            pz_eco_apu = alpha_eco*pz_osem_apu+(1-alpha_eco)*pz_cosem_apu;
                        end
                        if alpha_eco<=0.0096
                            pz_eco_apu = pz_cosem_apu;
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ECOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ECOSEM sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(4)
                        if verbose
                            tStart = tic;
                        end
                        pz_rm_apu = pz_rm_apu + lam(iter).*pz_rm_apu./(Summ+epps).*(A*(uu./(A'*pz_rm_apu+epps))-Summ);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['RAMLA sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(7)
                        if verbose
                            tStart = tic;
                        end
                        if Nz==1
                            med = medfilt2(reshape(pz_mrp_osl_apu,Nx,Ny),[options.medx options.medy], 'symmetric');
                        else
                            med = medfilt3(reshape(pz_mrp_osl_apu,Nx,Ny,Nz), [options.medx options.medy options.medz]);
                        end
                        med = med(:);
                        pz_mrp_osl_apu=((pz_mrp_osl_apu)./(Summ+options.beta_mrp_osl*(pz_mrp_osl_apu-med)./(med+epps))).*((A*(uu./(A'*pz_mrp_osl_apu+epps))));
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(8)
                        if verbose
                            tStart = tic;
                        end
                        pz_mrp_bsrem_apu = pz_mrp_bsrem_apu + lam(iter).*pz_mrp_bsrem_apu./(Summ+epps).*(A*(uu./(A'*pz_mrp_bsrem_apu+epps))-Summ);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(9)
                        if verbose
                            tStart = tic;
                        end
                        if Nz==1
                            med = bsxfun(@minus,pz_pad_osl(tr_offsets(:,isinf(weights))),pz_pad_osl(tr_offsets(:,[1:(find(isinf(weights))-1) (find(isinf(weights))+1):end])))*weights_quad;
                        else
                            med = bsxfun(@minus,pz_pad_osl(tr_offsets(:,isinf(weights))),pz_pad_osl(tr_offsets(:,[1:(find(isinf(weights))-1) (find(isinf(weights))+1):end])))*weights_quad;
                        end
                        pz_quad_osl_apu=(pz_quad_osl_apu./(Summ+options.beta_quad_osl*med)).*((A*(uu./(A'*pz_quad_osl_apu+epps))));
                        if Nz==1
                            pz_pad_osl = padding(reshape(pz_quad_osl_apu,Nx,Ny),[Ndx Ndy]);
                        else
                            pz_pad_osl = padding(reshape(pz_quad_osl_apu,Nx,Ny,Nz),[Ndx Ndy Ndz]);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(10)
                        if verbose
                            tStart = tic;
                        end
                        pz_quad_bsrem_apu = pz_quad_bsrem_apu + lam(iter).*pz_quad_bsrem_apu./(Summ+epps).*(A*(uu./(A'*pz_quad_bsrem_apu+epps))-Summ);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(11)
                        if verbose
                            tStart = tic;
                        end
                        med = sort(pz_pad_L_osl(tr_offsets),2)*a_L;
                        pz_L_osl_apu=((pz_L_osl_apu)./(Summ+options.beta_L_osl*(pz_L_osl_apu-med)./(med+epps))).*((A*(uu./(A'*pz_L_osl_apu+epps))));
                        if Nz==1
                            pz_pad_L_osl = padding(reshape(pz_L_osl_apu,Nx,Ny),[Ndx Ndy]);
                        else
                            pz_pad_L_osl = padding(reshape(pz_L_osl_apu,Nx,Ny,Nz),[Ndx Ndy Ndz]);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(12)
                        if verbose
                            tStart = tic;
                        end
                        pz_L_bsrem_apu = pz_L_bsrem_apu + lam(iter).*pz_L_bsrem_apu./(Summ+epps).*(A*(uu./(A'*pz_L_bsrem_apu+epps))-Summ);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(13)
                        if verbose
                            tStart = tic;
                        end
                        med = zeros(N, find(isinf(weights)));
                        for ii = 1 : find(isinf(weights)) - 1
                            med(:,ii) = pz_pad_fmh(tr_offsets(:,ii:(find(isinf(weights)) - ii):end + 1 - ii))*fmh_weights(:,ii);
                        end
                        med(:,end) = pz_pad_fmh(tr_offsets(:,isinf(weights)));
                        med = median(med,2);
                        pz_fmh_osl_apu=((pz_fmh_osl_apu)./(Summ+options.beta_fmh_osl*(pz_fmh_osl_apu-med)./(med+epps))).*((A*(uu./(A'*pz_fmh_osl_apu+epps))));
                        if Nz==1
                            pz_pad_fmh = padding(reshape(pz_fmh_osl_apu,Nx,Ny),[Ndx Ndy]);
                        else
                            pz_pad_fmh = padding(reshape(pz_fmh_osl_apu,Nx,Ny,Nz),[Ndx Ndy Ndz]);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    if rekot(14)
                        if verbose
                            tStart = tic;
                        end
                        med = pz_pad_weighted(tr_offsets)*weighted_weights;
                        pz_weighted_osl_apu=((pz_weighted_osl_apu)./(Summ+options.beta_weighted_osl*(pz_weighted_osl_apu-med)./(med+epps))).*((A*(uu./(A'*pz_weighted_osl_apu+epps))));
                        if Nz==1
                            pz_pad_weighted = padding(reshape(pz_weighted_osl_apu,Nx,Ny),[Ndx Ndy]);
                        else
                            pz_pad_weighted = padding(reshape(pz_weighted_osl_apu,Nx,Ny,Nz),[Ndx Ndy Ndz]);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['OSL weighted mean sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    clear A
                end
                if rekot(2)
                    pz_osem(:, iter + 1) = pz_osem_apu;
                end
                if rekot(3)
                    pz_rmM(:, iter + 1) = pz_rmM_apu;
                end
                if rekot(4)
                    pz_rm(:, iter + 1) = pz_rm_apu;
                end
                if rekot(5)
                    pz_eco(:, iter + 1) = pz_eco_apu;
                end
                if rekot(6)
                    pz_cosem(:, iter + 1) = pz_cosem_apu;
                end
                if rekot(15)
                    pz_acosem(:, iter + 1) = pz_acosem_apu;
                end
                if rekot(7)
                    pz_mrp_osl(:, iter + 1) = pz_mrp_osl_apu;
                end
                if rekot(8)
                    if verbose
                        tStart = tic;
                    end
                    if Nz==1
                        med = medfilt2(reshape(pz_mrp_bsrem_apu,Nx,Ny),[options.medx options.medy],'symmetric');
                    else
                        med = medfilt3(reshape(pz_mrp_bsrem_apu,Nx,Ny,Nz),[options.medx options.medy options.medz],'symmetric');
                    end
                    pz_mrp_bsrem(:,iter+1) = pz_mrp_bsrem_apu - options.beta_mrp_bsrem*lam(iter).*(pz_mrp_bsrem_apu-med(:)).*pz_mrp_bsrem_apu./(med(:)+epps);
                    pz_mrp_bsrem(pz_mrp_bsrem(:,iter+1)<0,iter+1) = 0;
                    if verbose
                        tElapsed = toc(tStart);
                        disp(['BSREM MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                    else
                        disp(['BSREM MRP iteration ' num2str(iter) ' finished'])
                    end
                end
                if rekot(9)
                    pz_quad_osl(:, iter + 1) = pz_quad_osl_apu;
                end
                if rekot(10)
                    if verbose
                        tStart = tic;
                    end
                    if Nz==1
                        pz_pad_bsrem = padding(reshape(pz_quad_bsrem_apu,Nx,Ny),[Ndx Ndy]);
                        med = -bsxfun(@minus,pz_pad_bsrem(tr_offsets(:,isinf(weights))),pz_pad_bsrem(tr_offsets(:,[1:(find(isinf(weights))-1) (find(isinf(weights))+1):end])))*weights_quad;
                    else
                        pz_pad_bsrem = padding(reshape(pz_quad_bsrem_apu,Nx,Ny,Nz),[Ndx Ndy Ndz]);
                        med = -bsxfun(@minus,pz_pad_bsrem(tr_offsets(:,isinf(weights))),pz_pad_bsrem(tr_offsets(:,[1:(find(isinf(weights))-1) (find(isinf(weights))+1):end])))*weights_quad;
                    end
                    pz_quad_bsrem(:,iter+1) = pz_quad_bsrem_apu.*(1 + options.beta_quad_bsrem*lam(iter).*med);
                    %                 pz_quad_bsrem(:,iter+1) = pz_quad_bsrem_apu - options.beta_quad_bsrem*lam(iter).*(pz_quad_bsrem_apu-med(:)).*pz_quad_bsrem_apu./(med(:)+epps);
                    pz_quad_bsrem(pz_quad_bsrem(:,iter+1)<0,iter+1) = 0;
                    if verbose
                        tElapsed = toc(tStart);
                        disp(['BSREM quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                    else
                        disp(['BSREM quadratic iteration ' num2str(iter) ' finished'])
                    end
                end
                if rekot(11)
                    pz_L_osl(:, iter + 1) = pz_L_osl_apu;
                end
                if rekot(12)
                    if verbose
                        tStart = tic;
                    end
                    if Nz==1
                        pz_pad_L_bsrem = padding(reshape(pz_L_bsrem_apu,Nx,Ny),[Ndx Ndy]);
                        med = sort(pz_pad_L_bsrem(tr_offsets),2)*a_L;
                    else
                        pz_pad_L_bsrem = padding(reshape(pz_L_bsrem_apu,Nx,Ny,Nz),[Ndx Ndy Ndz]);
                        med = sort(pz_pad_L_bsrem(tr_offsets),2)*a_L;
                    end
                    pz_L_bsrem(:,iter+1) = pz_L_bsrem_apu - options.beta_L_bsrem*lam(iter).*(pz_L_bsrem_apu-med(:)).*pz_L_bsrem_apu./med(:);
                    pz_L_bsrem(pz_L_bsrem(:,iter+1)<0,iter+1) = epps;
                    if verbose
                        tElapsed = toc(tStart);
                        disp(['BSREM L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                    else
                        disp(['BSREM L-filter iteration ' num2str(iter) ' finished'])
                    end
                end
                if rekot(13)
                    pz_fmh_osl(:, iter + 1) = pz_fmh_osl_apu;
                end
                if rekot(14)
                    pz_weighted_osl(:, iter + 1) = pz_weighted_osl_apu;
                end
                disp(['Iteration ' num2str(iter) ' finished'])
            end
        elseif options.reconstruction_method == 2
            if use_raw_data == false
                kernel_path = which('siddon_kernel.cl');
                [pz_osem] = improved_Siddon_openCL( kernel_path, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, single(NSlices), ...
                size_x, zmax, NSinos, vaimennus, pituus, int32(attenuation_correction), int32(Niter), int32(subsets), rekot, single(epps), ...
                single(full(Sino)), single(options.x0(:)), lor_a, summa, xy_index, z_index, options.verbose, device);
            else
                kernel_path = which('siddon_kernel_raw.cl');
                [pz_osem] = improved_Siddon_openCL_raw( kernel_path, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, single(NSlices), ...
                    size_x, zmax, NSinos, vector_size, vaimennus, int32(1), pituus, int32(attenuation_correction), int32(Niter), int32(subsets), rekot, single(epps), ...
                    single(full(Sino)), single(options.x0(:)), lor_a, summa, LL, pseudot, det_per_ring, options.verbose, device);
            end
        elseif options.reconstruction_method == 3
            for iter = 1 : Niter
                if rekot(1)
                    maksimi = int32(max(lor_a));
                    if use_raw_data == false
                        [Summ, rhs] = sequential_reconstruction( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, ...
                                size_x, zmax, vaimennus, length(Sino), attenuation_correction, lor_a, xy_index, z_index, NSinos, maksimi, epps, full(Sino), ...
                                pz_ml(:,iter), options.verbose);
                    else
                        [Summ, rhs] = sequential_reconstruction_raw( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , size_x, zmax, ...
                            vaimennus, length(Sino), attenuation_correction, lor_a, LL, pseudot, det_per_ring, maksimi, epps, full(Sino), pz_ml(:,iter),...
                            options.verbose);
                    end
                    if verbose
                        tStart = tic;
                    end
                    pz_ml(:,iter + 1) = pz_ml(:,iter)./(Summ + epps) .* rhs;
                    if verbose
                        tElapsed = toc(tStart);
                        disp(['MLEM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                    else
                        disp(['MLEM iteration ' num2str(iter) ' finished'])
                    end
                end
                if rekot(2)
                    pz_osem_apu = pz_osem(:,iter);
                    for osa_iter = 1 : subsets
                        maksimi = int32(max(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                        uu = double(full(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                        if use_raw_data == false
                            [Summ, rhs] = sequential_reconstruction( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, ...
                                size_x, zmax, vaimennus, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, ...
                                lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                                z_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), NSinos, maksimi, epps, uu, pz_osem_apu, options.verbose);
                        else
                            [Summ, rhs] = sequential_reconstruction_raw( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , size_x, zmax, ...
                                vaimennus, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                                LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2), pseudot, det_per_ring, maksimi, epps, uu, pz_osem_apu,...
                                options.verbose);
                        end
                        if verbose
                            tStart = tic;
                        end
                        pz_osem_apu = pz_osem_apu./(Summ + epps) .* rhs;
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['OSEM sub-iteration ' num2str(osa_iter) ' finished'])
                        end
                    end
                    pz_osem(:,iter + 1) = pz_osem_apu;
                end
            end
        elseif options.reconstruction_method == 4
            maksimi = zeros(subsets,1,'int32');
            for osa_iter = 1 : subsets
                maksimi(osa_iter) = int32(max(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))));
            end
            if use_raw_data == false
                kernel_path = which('siddon_kernel_matrixfree_GPU.cl');
                [pz] = improved_Siddon_openCL_matrixfree_GPU( kernel_path, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, ...
                    single(NSlices), size_x, zmax, NSinos, vaimennus, pituus, int32(attenuation_correction), int32(Niter), int32(subsets), ...
                    rekot, single(epps), single(full(Sino)), single(options.x0(:)), lor_a, xy_index, z_index, maksimi, options.verbose, device);
            else
                kernel_path = which('siddon_kernel_matrixfree_GPU_raw.cl');
                [pz] = improved_Siddon_openCL_matrixfree_GPU_raw( kernel_path, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, single(NSlices), ...
                    size_x, zmax, NSinos, vaimennus, pituus, int32(attenuation_correction), int32(Niter), int32(subsets), rekot, single(epps), ...
                    single(full(Sino)), single(options.x0(:)), lor_a, LL, pseudot, det_per_ring, maksimi, options.verbose, device);
            end
        else
            error('Error: Unsupported reconstruction method.');
        end
        clear mex
        if Nz == 1  && options.reconstruction_method ~= 4
            if rekot(1)
                pz_ml = reshape(pz_ml,Nx,Ny,Niter+1);
            end
            if rekot(2)
                pz_osem = reshape(pz_osem,Nx,Ny,Niter+1);
            end
            if options.reconstruction_method == 1
                if rekot(3)
                    pz_rmM = reshape(pz_rmM,Nx,Ny,Niter+1);
                end
                if rekot(4)
                    pz_rm = reshape(pz_rm,Nx,Ny,Niter+1);
                end
                if rekot(5)
                    pz_eco = reshape(pz_eco,Nx,Ny,Niter+1);
                end
                if rekot(6)
                    pz_cosem = reshape(pz_cosem,Nx,Ny,Niter+1);
                end
                if rekot(15)
                    pz_acosem = reshape(pz_acosem,Nx,Ny,Niter+1);
                end
                if rekot(7)
                    pz_mrp_osl = reshape(pz_mrp_osl,Nx,Ny,Niter+1);
                end
                if rekot(8)
                    pz_mrp_bsrem = reshape(pz_mrp_bsrem,Nx,Ny,Niter+1);
                end
                if rekot(9)
                    pz_quad_osl = reshape(pz_quad_osl,Nx,Ny,Niter+1);
                end
                if rekot(10)
                    pz_quad_bsrem = reshape(pz_quad_bsrem,Nx,Ny,Niter+1);
                end
                if rekot(11)
                    pz_L_osl = reshape(pz_L_osl,Nx,Ny,Niter+1);
                end
                if rekot(12)
                    pz_L_bsrem = reshape(pz_L_bsrem,Nx,Ny,Niter+1);
                end
                if rekot(13)
                    pz_fmh_osl = reshape(pz_fmh_osl,Nx,Ny,Niter+1);
                end
                if rekot(14)
                    pz_weighted_osl = reshape(pz_weighted_osl,Nx,Ny,Niter+1);
                end
            end
        elseif options.reconstruction_method ~= 4
            if rekot(1)
                pz_ml = reshape(pz_ml,Nx,Ny,Nz,Niter+1);
            end
            if rekot(2)
                pz_osem = reshape(pz_osem,Nx,Ny,Nz,Niter+1);
            end
            if options.reconstruction_method == 1
                if rekot(3)
                    pz_rmM = reshape(pz_rmM,Nx,Ny,Nz,Niter+1);
                end
                if rekot(4)
                    pz_rm = reshape(pz_rm,Nx,Ny,Nz,Niter+1);
                end
                if rekot(5)
                    pz_eco = reshape(pz_eco,Nx,Ny,Nz,Niter+1);
                end
                if rekot(6)
                    pz_cosem = reshape(pz_cosem,Nx,Ny,Nz,Niter+1);
                end
                if rekot(15)
                    pz_acosem = reshape(pz_acosem,Nx,Ny,Nz,Niter+1);
                end
                if rekot(7)
                    pz_mrp_osl = reshape(pz_mrp_osl,Nx,Ny,Nz,Niter+1);
                end
                if rekot(8)
                    pz_mrp_bsrem = reshape(pz_mrp_bsrem,Nx,Ny,Nz,Niter+1);
                end
                if rekot(9)
                    pz_quad_osl = reshape(pz_quad_osl,Nx,Ny,Nz,Niter+1);
                end
                if rekot(10)
                    pz_quad_bsrem = reshape(pz_quad_bsrem,Nx,Ny,Nz,Niter+1);
                end
                if rekot(11)
                    pz_L_osl = reshape(pz_L_osl,Nx,Ny,Nz,Niter+1);
                end
                if rekot(12)
                    pz_L_bsrem = reshape(pz_L_bsrem,Nx,Ny,Nz,Niter+1);
                end
                if rekot(13)
                    pz_fmh_osl = reshape(pz_fmh_osl,Nx,Ny,Nz,Niter+1);
                end
                if rekot(14)
                    pz_weighted_osl = reshape(pz_weighted_osl,Nx,Ny,Nz,Niter+1);
                end
            end
        end
        
    end
    gg = 1;
    if options.reconstruction_method == 3
        if rekot(gg)
            pz{gg, llo} =  pz_ml;
        end
    end
    gg = gg + 1;
    if rekot(gg) && options.reconstruction_method ~= 4
        pz{gg, llo} =  pz_osem;
    end
    if options.reconstruction_method == 1
        gg = gg + 1;
        if rekot(gg)
            pz{gg, llo} =  pz_rmM;
        end
        gg = gg + 1;
        if rekot(gg)
            pz{gg, llo} =  pz_rm;
        end
        gg = gg + 1;
        if rekot(gg)
            pz{gg, llo} =  pz_eco;
        end
        gg = gg + 1;
        if rekot(gg)
            pz{gg, llo} =  pz_cosem;
        end
        gg = gg + 1;
        if rekot(15)
            pz{gg, llo} =  pz_acosem;
        end
        gg = gg + 1;
        if rekot(gg-1)
            pz{gg, llo} =  pz_mrp_osl;
        end
        gg = gg + 1;
        if rekot(gg-1)
            pz{gg, llo} =  pz_mrp_bsrem;
        end
        gg = gg + 1;
        if rekot(gg-1)
            pz{gg, llo} =  pz_quad_osl;
        end
        gg = gg + 1;
        if rekot(gg-1)
            pz{gg, llo} =  pz_quad_bsrem;
        end
        gg = gg + 1;
        if rekot(gg-1)
            pz{gg, llo} =  pz_L_osl;
        end
        gg = gg + 1;
        if rekot(gg-1)
            pz{gg, llo} = pz_L_bsrem;
        end
        gg = gg + 1;
        if rekot(gg-1)
            pz{gg, llo} = pz_fmh_osl;
        end
        gg = gg + 1;
        if rekot(gg-1)
            pz{gg, llo} = pz_weighted_osl;
        end
    end
end

if options.reconstruction_method == 1
    image_properties.Reconstruction_method = 'MATLAB';
elseif options.reconstruction_method == 2
    image_properties.Reconstruction_method = 'OpenCL';
elseif options.reconstruction_method == 3
    image_properties.Reconstruction_method = 'Sequential matrix-free';
elseif options.reconstruction_method == 4
    image_properties.Reconstruction_method = 'Matrix-free OpenCL';
end
image_properties.Nx = Nx;
image_properties.Ny = Ny;
image_properties.Nz = Nz;
image_properties.Nang = Nang;
image_properties.Ndist = Ndist;
image_properties.NSinos = NSinos;
image_properties.Niter = Niter;
image_properties.subsets = subsets;
image_properties.FOV = FOVa;
image_properties.axial_FOV = axial_fov;
image_properties.raw_data = use_raw_data;
image_properties.name = options.name;
image_properties.machine_name = machine_name;
image_properties.n_time_steps = partitions;
image_properties.attenuation = attenuation_correction;
if options.reconstruction_method == 1
    if rekot(15)
        image_properties.h = options.h;
    end
    if rekot(3) || rekot(4) || rekot(8) || rekot(10)
        image_properties.b0 = options.b0;
    end
    if rekot(3)
        image_properties.U = U;
    end
    if rekot(7)
        image_properties.beta_mrp_osl = options.beta_mrp_osl;
    end
    if rekot(8)
        image_properties.beta_mrp_bsrem = options.beta_mrp_bsrem;
    end
    if rekot(8) || rekot(7)
        image_properties.medx = options.medx;
        image_properties.medy = options.medy;
        image_properties.medz = options.medz;
    end
    if rekot(9)
        image_properties.beta_quad_osl = options.beta_quad_osl;
    end
    if rekot(10)
        image_properties.beta_quad_bsrem = options.beta_quad_bsrem;
    end
    if rekot(9) || rekot(10) || rekot(11) || rekot(12)
        image_properties.weights = weights;
    end
    if rekot(11)
        image_properties.beta_L_osl = options.beta_L_osl;
    end
    if rekot(12)
        image_properties.beta_L_bsrem = options.beta_L_bsrem;
    end
    if rekot(9) || rekot(10)
        image_properties.quadratic_prior_weights = weights_quad;
    end
    if rekot(11) || rekot(12)
        image_properties.L_weights = a_L;
    end
    if rekot(9) || rekot(10) || rekot(11) || rekot(12)
        image_properties.Ndx = Ndx;
        image_properties.Ndy = Ndy;
        image_properties.Ndz = Ndz;
    end
    if rekot(13)
        image_properties.fmh_center_weight = options.fmh_center_weight;
        image_properties.fmh_weights = fmh_weights;
        image_properties.beta_fmh_osl = options.beta_fmh_osl;
    end
    if rekot(14)
        image_properties.weighted_mean_center_weight = options.weighted_center_weight;
        image_properties.beta_weighted_mean_osl = options.beta_weighted_osl;
        image_properties.weighted_mean_weights = weighted_weights;
    end
end

pz{end,1} = image_properties;

end
