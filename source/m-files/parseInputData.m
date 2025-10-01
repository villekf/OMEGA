function options = parseInputData(options, index)
%PARSEINPUTDATA Perform subset data division
%   Since certain subset types require the input measurement and correction
%   data to be ordered in a specific way, this function makes sure that the
%   data is correctly ordered. Not all subset types require this reordering
if numel(options.partitions) > 1
    partitions = numel(options.partitions);
else
    partitions = options.partitions;
end
if options.subsets > 1 && options.subset_type > 0
    if ~options.largeDim
        if partitions > 1
            for ff = 1 : partitions
                if ~options.use_raw_data
                    if iscell(options.SinM)
                        temp = options.SinM{ff};
                    else
                        if options.listmode == 0
                            if options.TOF
                                temp = options.SinM(:,:,:,:,ff);
                            else
                                temp = options.SinM(:,:,:,ff);
                            end
                        else
                            temp = options.SinM(:,ff);
                        end
                    end
                    if options.NSinos ~= options.TotSinos && options.listmode == 0
                        temp = temp(:,:,1:options.NSinos,:);
                    end
                else
                    temp = single(full(options.SinM{ff}));
                end
                if options.TOF && options.listmode == 0
                    if options.subset_type >= 8
                        temp = temp(:,:,index,:);
                    else
                        temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                        temp = temp(index,:);
                    end
                else
                    if options.subset_type >= 8
                        temp = temp(:,:,index);
                    else
                        temp = temp(index);
                    end
                end
                if iscell(options.SinM)
                    options.SinM{ff} = temp(:);
                else
                    if options.listmode == 0
                        if options.TOF
                            options.SinM(:,:,:,:,ff) = temp;
                        else
                            options.SinM(:,:,:,ff) = temp;
                        end
                    else
                        options.SinM(:,ff) = temp;
                    end
                end
            end
            clear temp
        else
            if ~options.use_raw_data
                if options.NSinos ~= options.TotSinos && options.listmode == 0
                    if iscell(options.SinM)
                        options.SinM{1} = options.SinM{1}(:,:,1:options.NSinos,:);
                    else
                        options.SinM = options.SinM(:,:,1:options.NSinos,:);
                    end
                end
            else
                options.SinM = single(full(options.SinM{1}));
            end
            if options.subset_type >= 8
                if iscell(options.SinM)
                    options.SinM{1} = reshape(options.SinM{1}, options.nRowsD, options.nColsD, options.nProjections, options.TOF_bins);
                else
                    options.SinM = reshape(options.SinM, options.nRowsD, options.nColsD, options.nProjections, options.TOF_bins);
                end
            end
            if options.TOF && options.listmode == 0
                if options.subset_type >= 8
                    if iscell(options.SinM)
                        options.SinM{1} = options.SinM{1}(:,:,index,:);
                    else
                        options.SinM = options.SinM(:,:,index,:);
                    end
                else
                    if iscell(options.SinM)
                        options.SinM{1} = reshape(options.SinM{1}, numel(options.SinM) / options.TOF_bins, options.TOF_bins);
                        options.SinM{1} = options.SinM{1}(index,:);
                    else
                        options.SinM = reshape(options.SinM, numel(options.SinM) / options.TOF_bins, options.TOF_bins);
                        options.SinM = options.SinM(index,:);
                    end
                end
            else
                if options.subset_type >= 8
                    if iscell(options.SinM)
                        options.SinM{1} = options.SinM{1}(:,:,index);
                    else
                        options.SinM = options.SinM(:,:,index);
                    end
                elseif options.subset_type > 0
                    if iscell(options.SinM)
                        options.SinM{1} = options.SinM{1}(index);
                    else
                        options.SinM = options.SinM(index);
                    end
                end
            end
            if iscell(options.SinM)
                options.SinM{1} = options.SinM{1}(:);
            else
                options.SinM = options.SinM(:);
            end
        end
    end
    if options.normalization_correction && options.corrections_during_reconstruction
        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
            options.normalization = options.normalization(1:options.NSinos*options.Ndist*options.Nang);
        end
        if options.subset_type >= 8
            options.normalization = reshape(options.normalization, options.Ndist, options.Nang, []);
            options.normalization = options.normalization(:,:,index,:);
            options.normalization = options.normalization(:);
        else
            options.normalization = options.normalization(index);
        end
    end
    if options.additionalCorrection && isfield(options,'corrVector') && numel(options.corrVector) > 1
        if options.subset_type >= 8
            if iscell(options.corrVector)
                for kk = 1 : numel(options.corrVector)
                    options.corrVector{kk} = reshape(options.corrVector{kk}, options.Ndist, options.Nang, []);
                    options.corrVector{kk} = options.corrVector{kk}(:,:,index,:);
                    options.corrVector{kk} = options.corrVector{kk}(:);
                end
            else
                options.corrVector = reshape(options.corrVector, options.Ndist, options.Nang, []);
                options.corrVector = options.corrVector(:,:,index,:);
                options.corrVector = options.corrVector(:);
            end
        else
            if iscell(options.corrVector)
                for kk = 1 : numel(options.corrVector)
                    options.corrVector{kk} = options.corrVector{kk}(index);
                end
            else
                options.corrVector = options.corrVector(index);
            end
        end
    end
    if ~options.largeDim
        if (options.randoms_correction || (options.scatter_correction && options.subtract_scatter)) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if partitions > 1
                for ff = 1 : partitions
                    if ~options.use_raw_data
                        if iscell(options.SinDelayed)
                            temp = options.SinDelayed{ff};
                        else
                            temp = options.SinDelayed(:,:,:,ff);
                        end
                        if options.NSinos ~= options.TotSinos
                            temp = temp(:,:,1:options.NSinos);
                        end
                    else
                        temp = single(full(options.SinDelayed{ff}));
                    end
                    if options.subset_type >= 8
                        temp = reshape(temp, options.nRowsD, options.nColsD, []);
                        temp = temp(:,:,index);
                        % temp = temp(:);
                    else
                        temp = temp(index);
                    end
                    if iscell(options.SinDelayed)
                        options.SinDelayed{ff} = temp;
                    else
                        options.SinDelayed(:,:,:,ff) = temp;
                    end
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    if options.subset_type >= 8
                        options.SinDelayed{1} = reshape(options.SinDelayed{1}, options.nRowsD, options.nColsD, []);
                        options.SinDelayed{1} = options.SinDelayed{1}(:,:,index);
                        options.SinDelayed{1} = options.SinDelayed{1}(:);
                    else
                        options.SinDelayed{1} = options.SinDelayed{1}(index);
                    end
                else
                    if options.subset_type >= 8
                        options.SinDelayed = reshape(options.SinDelayed, options.nRowsD, options.nColsD, []);
                        options.SinDelayed = options.SinDelayed(:,:,index);
                        options.SinDelayed = options.SinDelayed(:);
                    else
                        options.SinDelayed = options.SinDelayed(index);
                    end
                end
            end
        end
    end
    if options.attenuation_correction && ~options.CT_attenuation
        if numel(options.vaimennus) ~= options.Nx(1) * options.Ny(1) * options.Nz(1) || (options.nRowsD == options.Nx(1) && options.nColsD == options.Ny(1))
            if options.subset_type >= 8
                options.vaimennus = reshape(options.vaimennus, options.nRowsD, options.nColsD, []);
                options.vaimennus = options.vaimennus(:,:,index);
                options.vaimennus = options.vaimennus(:);
            else
                options.vaimennus = options.vaimennus(:);
                options.vaimennus = options.vaimennus(index);
            end
        else
            options.CT_attenuation = true;
        end
    end
    if options.useMaskFP && options.maskFPZ > 1 && options.subset_type >= 8
        options.maskFP = uint8(options.maskFP(:,:,index));
    elseif options.useMaskFP && options.subset_type == 3
        error('Forward projection mask is not supported with subset type 3!')
    end
end