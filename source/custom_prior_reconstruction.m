function options = custom_prior_reconstruction(options, t, iter, osa_iter)
%% Main reconstruction file for the custom prior file
% This function is used to compute various reconstructions with the
% selected methods and the custom prior gradient.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi, Samuli Summala
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.OS_bool && options.MLEM_bool
    error('Custom prior reconstruction is not supported when BOTH an OS-type algorithm and a non-OS algorithm are selected')
end

if iscell(options.SinM)
    Sino = options.SinM{t};
else
    Sino = options.SinM;
end

Sino = Sino(:);

if issparse(Sino)
    Sino = (full(Sino));
end

if ~isfield(options, 'attenuation_phase')
    options.attenuation_phase = false;
end
if ~isfield(options,'TOF_bins')
    options.TOF_bins = 1;
end
if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'compute_sensitivity_image')
    options.compute_sensitivity_image = false;
end
TOF = options.TOF_bins > 1 && options.projector_type == 1;

if t == 1 && osa_iter == 1 && iter == 1
    options.x00 = options.x0;
end

[gaussK, options] = PSFKernel(options);

if exist('feature','builtin') == 5
    nCores = uint32(feature('numcores'));
else
    nCores = uint32(1);
end
    
Nx = uint32(options.Nx);
Ny = uint32(options.Ny);
Nz = uint32(options.Nz);

%%
if options.implementation == 1
    
    iij = double(0:options.Nx);
    jji = double(0:options.Ny);
    kkj = double(0:options.Nz);
    if ~options.use_raw_data
        if isempty(options.pseudot)
            options.pseudot = int32(0);
        end
    end
    if options.randoms_correction
        if iscell(options.SinDelayed)
            SinD = options.SinDelayed{llo}(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
        else
            SinD = options.SinDelayed(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
        end
        if issparse(SinD)
            SinD = (full(SinD));
        end
        SinD = SinD(:);
    else
        SinD = 0;
    end
    
    if options.normalization_correction
        normalization = options.normalization(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
    else
        normalization = 0;
    end
    
    if options.scatter_correction && ~options.subtract_scatter
        if options.implementation == 1
            scatter_input = options.ScatterC(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
        else
            if iscell(options.SinDelayed)
                options.ScatterFB{1} = {single(options.ScatterC{1}(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)))};
            else
                options.ScatterFB{1} = {single(options.ScatterC(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)))};
            end
        end
    else
        scatter_input = 0;
    end
    
    uu = double(Sino(options.pituus(osa_iter) + 1 : options.pituus(osa_iter + 1)));
    koko = length(uu);
    [A,~, Summ] = computeImplementation1(options,options.use_raw_data,options.randoms_correction, options.pituus,osa_iter, options.normalization_correction,...
        options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, options.x, options.y, options.z_det, options.xx, options.yy, options.size_x, options.NSinos, options.NSlices, options.zmax, options.attenuation_correction, options.pseudot, options.det_per_ring, ...
        TOF, options.sigma_x, options.TOFCenter, options.dec, nCores, options.ind_size, options.block1, options.blocks, options.index, iij, jji, kkj, options.LL, options.N, options.summa, options.lor_a, options.xy_index, options.z_index, ...
        options.x_center, options.y_center, options.z_center, options.bmin, options.bmax, options.Vmax, options.V, options.lor_orth, gaussK,options.is_transposed, scatter_input, normalization, SinD, koko);
    [options.im_vectors,options.C_co,options.C_aco,options.C_osl] = computeEstimatesImp1(options.im_vectors, options, A, uu, Summ, SinD, options.is_transposed, gaussK, iter, osa_iter, options.C_co, options.C_aco,options.C_osl,...
        options.randoms_correction, options.N, options.Ndx, options.Ndy, options.Ndz, options.D);
    clear A
else
    %%
    if options.implementation == 3
        error('Implementation 3 is not supported with custom priors!')
    end
    options = double_to_single(options);
    rekot = reko_maker(options);
    
    if t == 1 && iter == 1 && osa_iter == 1
        options.im_vectors = initialize_im_vectors(options.im_vectors, iter, options);
    end
%     if t > 1 && osa_iter == 1 && iter == 1
        options.x0 = options.x00;
%     else
%         options.x0 = updateInitialValue(options.im_vectors, options);
%     end
    if (options.RBI || options.MBSREM || options.OSL_RBI || options.OSL_COSEM > 0) && t == 1 && iter == 1 && osa_iter == 1
        options.D = zeros(options.N,1,'single');
    end
    varMAP = recNames(2);
    
    oo = 1;
    for kk = 1 : numel(varMAP)
        if isfield(options, ['grad_' varMAP{kk}]) && options.(varMAP{kk})
            options.(['grad_' varMAP{kk}]) = single(options.(['grad_' varMAP{kk}])(:));
            options.(['beta_custom_' varMAP{kk}]) = single(options.(['beta_custom_' varMAP{kk}]));
            options.(['custom_' varMAP{kk} '_apu']) = options.im_vectors.(['custom_' varMAP{kk} '_apu']);
            options.varApu{oo} = ['custom_' varMAP{kk} '_apu'];
            options.varGrad{oo} = ['grad_' varMAP{kk}];
            options.varBeta{oo} = ['beta_custom_' varMAP{kk}];
            oo = oo + 1;
        end
    end
    
    fn = fieldnames(options.im_vectors);
    options.varTot = fn(~cellfun('isempty',strfind(fn,'apu')));
    
    
    %     if (options.COSEM || options.ECOSEM) && t == 1 && osa_iter == 1 && iter == 1
    %         options.C_co = zeros(options.N, options.subsets, 'single');
    %     end
    %     if options.ACOSEM && t == 1 && osa_iter == 1 && iter == 1
    %         options.C_aco = zeros(options.N, options.subsets, 'single');
    %     end
    if any(options.OSL_COSEM) && t == 1 && osa_iter == 1 && iter == 1
        options.C_osl = zeros(options.N, options.subsets, 'single');
    end
    
    if options.partitions == 1
        if ~iscell(options.SinM)
            options.SinM = {options.SinM};
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            apu = single(full(options.SinDelayed));
            options.SinDelayed = cell(1,1);
            options.SinDelayed{1} = apu;
        end
        clear apu
    end
    if issparse(options.SinM{1})
        for kk = 1 : length(options.SinM)
            options.SinM{kk} = single(full(options.SinM{kk}));
        end
    elseif isa(options.SinM{1}, 'uint16')
        for kk = 1 : length(options.SinM)
            options.SinM{kk} = single(options.SinM{kk});
        end
    end
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction ...
            && ~options.reconstruct_trues && ~options.reconstruct_scatter
        options.randSize = uint64(diff(pituus));
        if issparse(options.SinDelayed{1})
            for kk = 1 : length(options.SinDelayed)
                options.SinDelayed{kk} = single(full(options.SinDelayed{kk}));
            end
        end
    else
        options.randSize = uint64(1);
    end
    options.tt = uint32(t - 1);
    options.iter = uint32(iter - 1);
    options.osa_iter = uint32(osa_iter - 1);
%     rekot = [rekot; false; false; false; false];
    
    [pz] = computeImplementation23(options, Ny, Nx, Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx, uint32(options.NSinos), options.NSlices, options.size_x, options.zmax, ...
        options.LL, options.pseudot, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, options.dec, uint32(options.use_device), options.use_raw_data, options.normalization, options.pituus, options.attenuation_correction, ...
        options.normalization_correction, uint32(iter), options.subsets, options.epps, options.lor_a, options.xy_index, options.z_index, options.x_center, options.y_center, options.z_center, options.SinDelayed, ...
        options.SinM, options.bmin, options.bmax, options.Vmax, options.V, gaussK, -1, rekot);
    
    options.im_vectors = transfer_im_vectors(options.im_vectors, pz, options, iter);
    %     if (options.COSEM || options.ECOSEM)
    %         options.C_co = pz{end-3};
    %     end
    %     if options.ACOSEM
    %         options.C_aco = pz{end-2};
    %     end
    if any(options.OSL_COSEM) && osa_iter == 1 && iter == 1 && t == 1
        options.C_osl = pz{end-5};
    end
    if (options.MRAMLA || options.MBSREM || options.OSL_RBI || options.RBI || options.COSEM || options.ECOSEM...
            || options.ACOSEM || any(options.OSL_COSEM)) && options.MBSREM_prepass && osa_iter == 1 && iter == 1 && t == 1
        options.D = pz{end-4}(:);
    end
    if options.MBSREM && osa_iter == 1 && iter == 1 && t == 1
        options.epsilon_mramla = pz{end-3};
        options.U = pz{end - 2};
    end
    pz{end} = 0;
    
end
