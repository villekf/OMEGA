function [normalization_correction, randoms_correction, options] = set_up_corrections(options, rings, RandProp, ScatterProp)
%% SET UP CORRECTIONS
% This function sets up (loads) the necessary corrections that are applied.
% Included are attenuation correction, normalization correction and randoms
% and/or scatter correction. The variables used for correction are saved in
% the options-struct that also contains the input variables except for the
% number of rings in the current machine.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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

folder = fileparts(which('set_up_corrections.m'));
folder = [folder(1:end-6), 'mat-files/'];
folder = strrep(folder, '\','/');

if ~isfield(options,'reconstruct_trues')
    options.reconstruct_trues = false;
end
if ~isfield(options,'reconstruct_scatter')
    options.reconstruct_scatter = false;
end
if ~isfield(options,'TOF_bins') || options.TOF_bins == 0
    options.TOF_bins = 1;
end
if ~isfield(options,'use_machine')
    options.use_machine = 0;
end

if ~options.use_raw_data && options.randoms_correction
    if options.partitions == 1
        load_string = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
            num2str(options.TotSinos) '_span' num2str(options.span)];
        if options.use_machine == 0
            sinoFile = [load_string '.mat'];
        elseif options.use_machine == 1
            sinoFile = [load_string '_listmode.mat'];
        elseif options.use_machine == 2
            sinoFile = [load_string '_machine_sinogram.mat'];
        elseif options.use_machine == 3
            sinoFile = [load_string '_listmode_sinogram.mat'];
        end
    else
        load_string = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_' ...
            num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
            num2str(options.span)];
        load_string2 = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
            num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
            num2str(options.span)];
        if options.use_machine == 0
            sinoFile = [load_string '.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '.mat'];
            end
        elseif options.use_machine == 1
            sinoFile = [load_string '_listmode.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_listmode.mat'];
            end
        elseif options.use_machine == 2
            sinoFile = [load_string '_machine_sinogram.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_machine_sinogram.mat'];
            end
        elseif options.use_machine == 3
            sinoFile = [load_string '_listmode_sinogram.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_listmode_sinogram.mat'];
            end
        end
    end
end


% options.use_psf = false;
block1 = 0;
if options.attenuation_correction
    if ~isfield(options,'vaimennus')
        data = load(options.attenuation_datafile);
        variables = fieldnames(data);
        options.vaimennus = double(data.(variables{1}));
        if size(options.vaimennus,1) ~= options.Nx || size(options.vaimennus,2) ~= options.Ny || size(options.vaimennus,3) ~= options.Nz
            if size(options.vaimennus,1) ~= options.Nx*options.Ny*options.Nz
                error('Error: Attenuation data is of different size than the reconstructed image')
            end
        end
        if rings > 0
            if size(options.vaimennus,3) == 1
                options.vaimennus = options.vaimennus(2*block1+1:(2*rings+1)*options.Nx*options.Ny);
            else
                options.vaimennus = options.vaimennus(:,:,2*block1+1:2*rings+1);
            end
        end
        options.vaimennus = options.vaimennus(:) / 10;
        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
            options.vaimennus = single(options.vaimennus);
        end
        clear data
    end
else
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
        options.vaimennus = single(0);
    else
        options.vaimennus = 0;
    end
end

if options.scatter_correction && options.normalize_scatter
    if ~isfield(options,'normalization')
        if options.use_user_normalization
            [file, fpath] = uigetfile({'*.nrm;*.mat'},'Select normalization datafile');
            if isequal(file, 0)
                error('No file was selected')
            end
            if any(strfind(file, '.nrm'))
                fid = fopen(fullfile(fpath, file));
                options.normalization = fread(fid, inf, 'single=>single',0,'l');
                fclose(fid);
                options.InveonNorm = true;
            else
                data = load(fullfile(fpath, file));
                variables = fieldnames(data);
                options.normalization = data.(variables{1});
                clear data
                if numel(options.normalization) ~= numel(options.SinM{1})
                    error('Size mismatch between the current data and the normalization data file')
                end
            end
            options.normalization = options.normalization(:);
            %         if any(strfind(file, '.nrm'))
            %             normalization = 1 ./ normalization;
            %         end
        else
            if options.use_raw_data
                norm_file = [folder options.machine_name '_normalization_listmode.mat'];
            else
                norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(options.Nang) '_span' num2str(options.span) '.mat'];
            end
            if exist(norm_file, 'file') == 2
                options.normalization = loadStructFromFile(norm_file,'normalization');
            else
                error('Normalization correction selected, but no normalization data found')
                %             normalization_coefficients(options);
                %             normalization = loadStructFromFile(norm_file,'normalization');
            end
            options.normalization = options.normalization(:);
        end
        if ~options.use_raw_data && options.NSinos ~= options.TotSinos
            options.normalization = options.normalization(1 : options.Ndist * options.Nang *options.NSinos);
        end
        if options.sampling > 1 && options.corrections_during_reconstruction
            options.normalization = reshape(options.normalization, options.Ndist, options.Nang, options.NSinos);
            options.normalization = interpolateSinog(options.normalization, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
            options.normalization = options.normalization(:);
        end
    end
end

% Apply randoms and scatter (if available) correction
% Apply randoms smoothing/variance reduction
if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction && ~options.reconstruct_trues...
        && ~options.reconstruct_scatter
    if ~options.randoms_correction && options.scatter_correction && ~options.subtract_scatter
        randoms_correction = false;
    else
        randoms_correction = true;
    end
    r_exist = isfield(options,'SinDelayed');
    if r_exist && isfield(options,'SinDelayed') && ~iscell(options.SinDelayed) && numel(options.SinDelayed) == 1 && randoms_correction
        r_exist = false;
    elseif r_exist && isfield(options,'SinDelayed') && iscell(options.SinDelayed) && numel(options.SinDelayed{1}) == 1 && randoms_correction
        r_exist = false;
    end
    if options.scatter_correction && ~options.subtract_scatter
        options.scatter = true;
    else
        options.scatter = false;
    end
    if exist('RandProp','var') == 0 || isempty(RandProp)
        RandProp.smoothing = false;
        RandProp.variance_reduction = false;
    end
    if exist('ScatterProp','var') == 0 || isempty(ScatterProp)
        ScatterProp.smoothing = false;
        ScatterProp.variance_reduction = false;
        ScatterProp.normalization = false;
    end
    if r_exist && randoms_correction
        if (~options.variance_reduction && RandProp.variance_reduction) || (~options.randoms_smoothing && RandProp.smoothing)
            if options.partitions == 1
                options.SinDelayed = loadStructFromFile(sinoFile,'raw_SinDelayed');
            else
                options.SinDelayed = loadStructFromFile(sinoFile, 'raw_SinDelayed');
            end
            RandProp.variance_reduction = false;
            RandProp.smoothing = false;
        end
    end
    if options.randoms_correction && ~r_exist
        options = loadDelayedData(options);
        r_exist = true;
    end
    if options.scatter_correction && (~isfield(options,'ScatterC') || (ScatterProp.normalization && ~options.normalize_scatter) || ...
            (ScatterProp.smoothing && ~options.scatter_smoothing) || (ScatterProp.variance_reduction && ~options.scatter_variance_reduction))
        if ScatterProp.normalization && ~options.normalize_scatter
            warning('Scatter data is normalized, but no normalized scatter selected. Input non-normalized scatter data!')
        end
        if ScatterProp.smoothing && ~options.scatter_smoothing
            warning('Scatter data is smoothed, but no scatter smoothing selected. Input non-smoothed scatter data!')
        end
        if ScatterProp.variance_reduction && ~options.scatter_variance_reduction
            warning('Scatter data is variance corrected, but no scatter variance correction selected. Input non-variance corrected scatter data!')
        end
        options = loadScatterData(options);
    end
    if r_exist
        if iscell(options.SinDelayed)
            for kk = 1 : options.partitions
                % SinM{kk} = SinM{kk} + SinDelayed{kk};
                % SinDelayed{kk} = 2 * SinDelayed{kk};
                if options.variance_reduction && ~RandProp.variance_reduction
                    options.SinDelayed{kk} = Randoms_variance_reduction(single(options.SinDelayed{kk}), options);
                end
                if options.randoms_smoothing && ~RandProp.smoothing
                    options.SinDelayed{kk} = randoms_smoothing(single(options.SinDelayed{kk}), options);
                end
                options.SinDelayed{kk} = options.SinDelayed{kk}(:);
                if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                    options.SinDelayed{kk} = options.SinDelayed{kk}(1:options.NSinos*options.Ndist*options.Nang);
                end
                if options.scatter_correction
                    if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                        if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data
                            options.ScatterC = permute(options.ScatterC,[2 1 3]);
                        end
                        if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                            if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                                options.ScatterC = single(options.ScatterC) .* reshape(options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                            else
                                options.ScatterC = single(options.ScatterC) .* reshape(1 ./ options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                            end
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                        end
                    elseif iscell(options.ScatterC)
                        if sum(size(options.ScatterC{1})) > 1
                            if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data
                                options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3]);
                            end
                        end
                        if options.normalization_correction && options.normalize_scatter
                            if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                            else
                                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(1 ./ options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                            end
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                        end
                    end
                    if options.subtract_scatter
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                            options.SinDelayed{kk} = single(options.SinDelayed{kk}) / options.TOF_bins + single(options.ScatterC{kk}(:));
                        else
                            options.SinDelayed{kk} = double(options.SinDelayed{kk}) / options.TOF_bins + double(options.ScatterC{kk}(:));
                        end
                    end
                elseif ~options.scatter_correction && options.TOF_bins > 1
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                        options.SinDelayed{kk} = single(options.SinDelayed{kk} / options.TOF_bins);
                    else
                        options.SinDelayed{kk} = double(options.SinDelayed{kk} / options.TOF_bins);
                    end
                end
            end
        else
            if options.variance_reduction && ~RandProp.variance_reduction
                options.SinDelayed = Randoms_variance_reduction(single(options.SinDelayed), options);
            end
            if options.randoms_smoothing && ~RandProp.smoothing
                options.SinDelayed = randoms_smoothing(single(options.SinDelayed), options);
            end
            options.SinDelayed = options.SinDelayed(:);
            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                options.SinDelayed = options.SinDelayed(1:options.NSinos*options.Ndist*options.Nang);
            end
            % Sino = Sino + SinDelayed;
            % SinDelayed = 2 * SinDelayed;
            if options.scatter_correction && isfield(options,'ScatterC')
                if isempty(ScatterProp)
                    ScatterProp.smoothing = false;
                    ScatterProp.variance_reduction = false;
                    ScatterProp.normalization = false;
                end
                if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                    if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data
                        options.ScatterC = permute(options.ScatterC,[2 1 3]);
                    end
                    if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                        if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                            options.ScatterC = single(options.ScatterC) .* reshape(options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                        else
                            options.ScatterC = single(options.ScatterC) .* reshape(1 ./ options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                        end
                    end
                    if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                        options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                    end
                    if options.scatter_smoothing && ~ScatterProp.smoothing
                        options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                    end
                    if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                        options.ScatterC = options.ScatterC(1:options.NSinos*options.Ndist*options.Nang);
                    end
                    if options.subtract_scatter
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                            options.SinDelayed = single(options.SinDelayed / options.TOF_bins) + single(options.ScatterC(:));
                        else
                            options.SinDelayed = double(options.SinDelayed / options.TOF_bins) + double(options.ScatterC(:));
                        end
                    end
                elseif iscell(options.ScatterC)
                    for kk = 1 : options.partitions
                        if sum(size(options.ScatterC{kk})) > 1
                            if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data
                                options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3]);
                            end
                        end
                        if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                            if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                            else
                                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(1 ./ options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                            end
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                        end
                        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                            options.ScatterC{kk} = options.ScatterC{kk}(1:options.NSinos*options.Ndist*options.Nang);
                        end
                        if options.subtract_scatter
                            if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                                options.SinDelayed = single(options.SinDelayed / options.TOF_bins) + single(options.ScatterC{kk}(:));
                            else
                                options.SinDelayed = double(options.SinDelayed / options.TOF_bins) + double(options.ScatterC{kk}(:));
                            end
                        end
                    end
                end
            elseif ~options.scatter_correction && options.TOF_bins > 1
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                    options.SinDelayed = single(options.SinDelayed / options.TOF_bins);
                else
                    options.SinDelayed = double(options.SinDelayed / options.TOF_bins);
                end
            end
        end
    else
        if options.randoms_correction && iscell(options.SinDelayed)
            for kk = 1 : options.partitions
                % SinM{kk} = SinM{kk} + SinDelayed{kk};
                % SinDelayed{kk} = 2 * SinDelayed{kk};
                if options.variance_reduction && ~RandProp.variance_reduction
                    options.SinDelayed{kk} = Randoms_variance_reduction(single(options.SinDelayed{kk}), options);
                end
                if options.randoms_smoothing && ~RandProp.smoothing
                    options.SinDelayed{kk} = randoms_smoothing(single(options.SinDelayed{kk}), options);
                end
                options.SinDelayed{kk} = options.SinDelayed{kk}(:);
                if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                    options.SinDelayed{kk} = options.SinDelayed{kk}(1:options.NSinos*options.Ndist*options.Nang);
                end
                if options.scatter_correction
                    if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                        if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data
                            options.ScatterC = permute(options.ScatterC,[2 1 3]);
                        end
                        if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                            if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                                options.ScatterC = single(options.ScatterC) .* reshape(options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                            else
                                options.ScatterC = single(options.ScatterC) .* reshape(1 ./ options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                            end
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                        end
                    elseif iscell(options.ScatterC)
                        if sum(size(options.ScatterC{1})) > 1
                            if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data
                                options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3]);
                            end
                        end
                        if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                            if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                            else
                                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(1 ./ options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                            end
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                        end
                    end
                    if options.subtract_scatter
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                            options.SinDelayed{kk} = single(options.SinDelayed{kk} / options.TOF_bins) + single(options.ScatterC{kk}(:));
                        else
                            options.SinDelayed{kk} = double(options.SinDelayed{kk} / options.TOF_bins) + double(options.ScatterC{kk}(:));
                        end
                    end
                elseif ~options.scatter_correction && options.TOF_bins > 1
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                        options.SinDelayed{kk} = single(options.SinDelayed{kk}) / single(options.TOF_bins);
                    else
                        options.SinDelayed{kk} = double(options.SinDelayed{kk}) / options.TOF_bins;
                    end
                end
            end
        elseif options.randoms_correction
            if options.variance_reduction && ~RandProp.variance_reduction
                options.SinDelayed = Randoms_variance_reduction(single(options.SinDelayed), options);
            end
            if options.randoms_smoothing && ~RandProp.smoothing
                options.SinDelayed = randoms_smoothing(single(options.SinDelayed), options);
            end
            options.SinDelayed = options.SinDelayed(:);
            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                options.SinDelayed = options.SinDelayed(1:options.NSinos*options.Ndist*options.Nang);
            end
            % Sino = Sino + SinDelayed;
            % SinDelayed = 2 * SinDelayed;
            if options.scatter_correction && isfield(options,'ScatterC')
                if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                    if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data
                        options.ScatterC = permute(options.ScatterC,[2 1 3]);
                    end
                    % if options.variance_reduction
                    %    options.ScatterC = Randoms_variance_reduction(double(options.ScatterC), options);
                    % end
                    if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                        if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                            options.ScatterC = single(options.ScatterC) .* reshape(options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                        else
                            options.ScatterC = single(options.ScatterC) .* reshape(1 ./ options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                        end
                    end
                    if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                        options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                    end
                    if options.scatter_smoothing && ~ScatterProp.smoothing
                        options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                    end
                    if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                        options.ScatterC = options.ScatterC(1:options.NSinos*options.Ndist*options.Nang);
                    end
                    if options.subtract_scatter
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                            options.SinDelayed = single(options.SinDelayed / options.TOF_bins) + single(options.ScatterC(:));
                        else
                            options.SinDelayed = double(options.SinDelayed / options.TOF_bins) + double(options.ScatterC(:));
                        end
                    end
                elseif iscell(options.ScatterC)
                    for kk = 1 : options.partitions
                        if sum(size(options.ScatterC{kk})) > 1
                            if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data
                                options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3]);
                            end
                        end
                        if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                            if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                            else
                                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(1 ./ options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                            end
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                        end
                        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                            options.ScatterC{kk} = options.ScatterC{kk}(1:options.NSinos*options.Ndist*options.Nang);
                        end
                        if options.subtract_scatter
                            if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                                options.SinDelayed = single(options.SinDelayed / options.TOF_bins) + single(options.ScatterC{kk}(:));
                            else
                                options.SinDelayed = double(options.SinDelayed / options.TOF_bins) + double(options.ScatterC{kk}(:));
                            end
                        end
                    end
                end
            elseif ~options.scatter_correction && options.TOF_bins > 1
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                    options.SinDelayed = single(options.SinDelayed) / single(options.TOF_bins);
                else
                    options.SinDelayed = double(options.SinDelayed) / options.TOF_bins;
                end
            end
        elseif options.scatter_correction
            if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data
                    options.ScatterC = permute(options.ScatterC,[2 1 3]);
                end
                if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                    if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                        options.ScatterC = single(options.ScatterC) .* reshape(options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                    else
                        options.ScatterC = single(options.ScatterC) .* reshape(1 ./ options.normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                    end
                end
                if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                    options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                end
                if options.scatter_smoothing && ~ScatterProp.smoothing
                    options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                end
                if options.subtract_scatter
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                        options.SinDelayed = single(options.ScatterC(:));
                    else
                        options.SinDelayed = double(options.ScatterC(:));
                    end
                end
            elseif iscell(options.ScatterC)
                options.SinDelayed = cell(length(options.ScatterC),1);
                for kk = 1 : options.partitions
                    if sum(size(options.ScatterC{1})) > 1
                        if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data
                            options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3]);
                        end
                    end
                    if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                        if isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
                            options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                        else
                            options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(1 ./ options.normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                        end
                    end
                    if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                        options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                    end
                    if options.scatter_smoothing && ~ScatterProp.smoothing
                        options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                    end
                    if options.subtract_scatter
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                            options.SinDelayed{kk} = single(options.ScatterC{kk}(:));
                        else
                            options.SinDelayed{kk} = double(options.ScatterC{kk}(:));
                        end
                    end
                end
            end
        end
    end
    if options.sampling > 1 && options.corrections_during_reconstruction
        options.SinDelayed = reshape(options.SinDelayed, options.Ndist, options.Nang, options.NSinos);
        options.SinDelayed = interpolateSinog(options.SinDelayed, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
        options.SinDelayed = options.SinDelayed(:);
    end
    if ~options.subtract_scatter && options.scatter_correction && options.sampling > 1 && options.corrections_during_reconstruction
        options.ScatterC = reshape(options.ScatterC, options.Ndist, options.Nang, options.NSinos);
        options.ScatterC = interpolateSinog(options.ScatterC, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
        options.ScatterC = options.ScatterC(:);
    end
elseif (options.randoms_correction || options.scatter_correction) && ~options.reconstruct_trues && ~options.reconstruct_scatter
    r_exist = isfield(options,'SinDelayed');
    if r_exist && isfield(options,'SinDelayed') && ~iscell(options.SinDelayed) && numel(options.SinDelayed) == 1 && options.randoms_correction
        r_exist = false;
    elseif r_exist && isfield(options,'SinDelayed') && iscell(options.SinDelayed) && numel(options.SinDelayed{1}) == 1 && options.randoms_correction
        r_exist = false;
    end
    options.scatter = false;
    if exist('RandProp','var') == 0 || isempty(RandProp)
        RandProp.smoothing = false;
        RandProp.variance_reduction = false;
    end
    if exist('ScatterProp','var') == 0 || isempty(ScatterProp)
        ScatterProp.smoothing = false;
        ScatterProp.variance_reduction = false;
        ScatterProp.normalization = false;
    end
    if r_exist && options.randoms_correction
        if (~options.variance_reduction && RandProp.variance_reduction) || (~options.randoms_smoothing && RandProp.smoothing)
            if options.partitions == 1
                options.SinDelayed = loadStructFromFile(sinoFile,'raw_SinDelayed');
            else
                options.SinDelayed = loadStructFromFile(sinoFile, 'raw_SinDelayed');
            end
            RandProp.variance_reduction = false;
            RandProp.smoothing = false;
        end
    end
    if options.randoms_correction && ~r_exist
        try
            options.SinDelayed = loadStructFromFile(sinoFile, 'SinDelayed');
        catch
            options = loadDelayedData(options);
        end
    end
    if options.randoms_correction
        [options,~, ~] = applyScatterRandoms(options,RandProp, options.SinDelayed, 0);
    end
    randoms_correction = false;
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
        options.SinDelayed = {single(0)};
    else
        options.SinDelayed = {0};
    end
    if options.scatter_correction
        [options,~, ~] = applyScatterRandoms(options,ScatterProp, options.ScatterC, 1);
        options = rmfield(options,'ScatterC');
        options.scatter_correction = false;
    end
    if options.verbose
        disp('Precorrections completed')
    end
else
    randoms_correction = false;
    options.scatter_correction = false;
    options.scatter = false;
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
        options.SinDelayed = {single(0)};
    else
        options.SinDelayed = {0};
    end
end

% Load (or compute) options.normalization correction coefficients
if (options.normalization_correction && options.corrections_during_reconstruction) && ~options.use_user_normalization
    normalization_correction = true;
    if ~isfield(options,'normalization')
        if options.use_raw_data
            norm_file = [folder options.machine_name '_normalization_listmode.mat'];
            if exist(norm_file, 'file') == 2
                options.normalization = loadStructFromFile(norm_file,'normalization');
            else
                error('Normalization correction selected, but no normalization data found')
            end
        else
            norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(options.Nang) '_span' num2str(options.span) '.mat'];
            if exist(norm_file, 'file') == 2
                options.normalization = loadStructFromFile(norm_file,'normalization');
            else
                error('Normalization correction selected, but no normalization data found')
            end
        end
        options.normalization(options.normalization <= 0) = 1;
        options.normalization = 1./options.normalization(:);
        options.normalization = single(options.normalization);
        if ~options.use_raw_data && options.NSinos ~= options.TotSinos
            options.normalization = options.normalization(1 : options.Ndist * options.Nang *options.NSinos);
        end
        if options.sampling > 1 && options.corrections_during_reconstruction
            options.normalization = reshape(options.normalization, options.Ndist, options.Nang, options.NSinos);
            options.normalization = interpolateSinog(options.normalization, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
            options.normalization = options.normalization(:);
        end
    end
elseif options.normalization_correction && options.use_user_normalization && options.corrections_during_reconstruction
    normalization_correction = true;
    if ~isfield(options,'normalization')
        [file, fpath] = uigetfile({'*.nrm;*.mat'},'Select normalization datafile');
        if isequal(file, 0)
            error('No file was selected')
        end
        if any(strfind(file, '.nrm'))
            if options.use_raw_data
                error('Inveon normalization data cannot be used with raw list-mode data')
            end
            fid = fopen(fullfile(fpath, file));
            options.normalization = fread(fid, inf, 'single=>single',0,'l');
            fclose(fid);
            if numel(options.normalization) ~= options.Ndist * options.Nang * options.TotSinos && ~options.use_raw_data
                error('Size mismatch between the current data and the normalization data file')
            end
            options.InveonNorm = true;
            %             options.normalization = 1 ./ options.normalization;
        else
            data = load(fullfile(fpath, file));
            variables = fieldnames(data);
            if any(strcmp(variables,'normalization'))
                options.normalization = data.(variables{strcmp(variables,'normalization')});
            else
                options.normalization = data.(variables{1});
            end
            clear data
            if numel(options.normalization) ~= options.Ndist * options.Nang * options.TotSinos && ~options.use_raw_data
                error('Size mismatch between the current data and the normalization data file')
            elseif numel(options.normalization) ~= sum(1:options.detectors) && options.use_raw_data
                error('Size mismatch between the current data and the normalization data file')
            end
        end
        options.normalization = single(options.normalization);
        options.normalization = 1./options.normalization(:);
        if ~options.use_raw_data && options.NSinos ~= options.TotSinos
            options.normalization = options.normalization(1 : options.Ndist * options.Nang *options.NSinos);
        end
        if options.sampling > 1 && options.corrections_during_reconstruction
            options.normalization = reshape(options.normalization, options.Ndist, options.Nang, options.NSinos);
            options.normalization = interpolateSinog(options.normalization, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
            options.normalization = options.normalization(:);
        end
    end
elseif options.normalization_correction && ~options.corrections_during_reconstruction
    normalization_correction = false;
    if ~isfield(options,'normalization')
        if options.use_user_normalization
            [file, fpath] = uigetfile({'*.nrm;*.mat'},'Select normalization datafile');
            if isequal(file, 0)
                error('No file was selected')
            end
            if any(strfind(file, '.nrm'))
                fid = fopen(fullfile(fpath, file));
                options.normalization = fread(fid, inf, 'single=>single',0,'l');
                fclose(fid);
                options.InveonNorm = true;
            else
                data = load(fullfile(fpath, file));
                variables = fieldnames(data);
                if any(strcmp(variables,'normalization'))
                    options.normalization = data.(variables{strcmp(variables,'normalization')});
                else
                    options.normalization = data.(variables{1});
                end
                clear data
                if iscell(options.SinM)
                    if numel(options.normalization) ~= numel(options.SinM{1})
                        error('Size mismatch between the current data and the normalization data file')
                    end
                else
                    if numel(options.normalization) ~= numel(options.SinM)
                        error('Size mismatch between the current data and the normalization data file')
                    end
                end
            end
            options.normalization = options.normalization(:);
            %         if any(strfind(file, '.nrm'))
            %             normalization = 1 ./ normalization;
            %         end
        else
            if options.use_raw_data
                norm_file = [folder options.machine_name '_normalization_listmode.mat'];
            else
                norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(options.Nang) '_span' num2str(options.span) '.mat'];
            end
            if exist(norm_file, 'file') == 2
                options.normalization = loadStructFromFile(norm_file,'normalization');
            else
                error('Normalization correction selected, but no normalization data found')
                %             normalization_coefficients(options);
                %             normalization = loadStructFromFile(norm_file,'normalization');
            end
            options.normalization = options.normalization(:);
        end
        if ~options.use_raw_data && options.NSinos ~= options.TotSinos
            options.normalization = options.normalization(1 : options.Ndist * options.Nang *options.NSinos);
        end
    end
    NN = options.Ndist * options.Nang * options.NSinos;
    if iscell(options.SinM)
        for kk = 1 : options.partitions
            options.SinM{kk} = options.SinM{kk}(:);
            if options.TOF_bins > 1
                for uu = 1 : options.TOF_bins
                    options.SinM{kk}(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM{kk}(NN * (uu - 1) + 1: NN * uu))) .* single(options.normalization);
                end
            else
                options.SinM{kk} = single(full(options.SinM{kk})) .* single(options.normalization);
            end
        end
    else
        if options.TOF_bins > 1
            for uu = 1 : options.TOF_bins
                options.SinM(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM(NN * (uu - 1) + 1: NN * uu))) .* single(options.normalization);
            end
        else
            options.SinM = options.SinM(:);
            options.SinM = single(full(options.SinM)) .* single(options.normalization);
        end
    end
    options.normalization = single(0);
else
    normalization_correction = false;
    options.normalization = single(0);
end

if ~isfield(options,'global_correction_factor') || isempty(options.global_correction_factor) || options.global_correction_factor <= 0
    options.global_correction_factor = 1;
end
if options.use_raw_data && ~options.corrections_during_reconstruction && ~options.listmode
    if iscell(options.SinM)
        for kk = 1 : options.partitions
            options.SinM{kk} = options.SinM{kk} * options.global_correction_factor;
        end
    else
        options.SinM = options.SinM * options.global_correction_factor;
    end
end
if options.implementation == 2 || options.implementation == 3
    options.global_correction_factor = single(options.global_correction_factor);
end
if ~isfield(options, 'ScatterC')
    if options.implementation == 2 || options.implementation == 3
        options.ScatterC = {single(0)};
    else
        options.ScatterC = 0;
    end
end
if ~isfield(options, 'scatter')
    options.scatter = false;
end
if (options.implementation == 2 || options.implementation == 3) && isfield(options,'deblur_iterations')
    options.deblur_iterations = uint32(options.deblur_iterations);
end
