function options = printReconstructionSettings(options)
%% Print the selected reconstruction settings
% Prints the selected data type, implementation (and the device), the
% reconstruction algorithm(s), prior(s), projector, corrections, subset
% settings and the image size if options.verbose > 0. Does not check the
% validity of the selections (see OMEGA_error_check and
% checkAlgorithmSelection).
%
% INPUT:
%   options = The options-struct (with the default values set, see
%   setMissingValues)
%
% OUTPUT:
%   options = The options-struct (returned unchanged).
%
% See also OMEGA_error_check, checkAlgorithmSelection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
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
if options.verbose <= 0
    return
end
MAP = checkAlgorithmsPriors(options, 2);
PRIOR = checkAlgorithmsPriors(options, 1);
if numel(options.partitions) > 1
    partitions = numel(options.partitions);
else
    partitions = options.partitions;
end
listModeCandidate = isListModeCandidate(options);
% Print various options that were selected if verbosity has been enabled
if options.verbose > 0
    if options.use_ASCII && options.use_machine == 0
        dispi = 'Using ASCII data';
    elseif options.use_LMF && options.use_machine == 0
        dispi = 'Using LMF data';
    elseif options.use_root && options.use_machine == 0
        dispi = 'Using ROOT data';
    elseif options.use_binary && options.use_machine == 0
        dispi = 'Using BINARY data';
    elseif options.use_machine == 1
        dispi = 'Using data obtained from list-mode file';
    elseif options.use_machine == 2
        dispi = 'Using scanner created sinogram data';
    elseif options.use_machine == 3
        dispi = 'Using 32-bit list-mode data';
    else
        dispi = 'Using input data';
    end
    if options.TOF_bins_used > 1 && options.TOF_FWHM > 0
        dispi = [dispi ' with TOF (' num2str(options.TOF_bins_used) ' bins).'];
    else
        dispi = [dispi '.'];
    end
    if ~strcmp(dispi,'.')
        disp(dispi);
    end
    if options.only_sinos
        disp('Loading only data.')
    end
    if ~(options.compute_normalization || options.only_sinos)
        dispaus = (['Using implementation ' num2str(options.implementation)]);
        if options.implementation == 1 || options.implementation == 4
            dispaus = [dispaus, '.'];
        elseif options.implementation == 2
            if ~options.use_CUDA && ~options.use_CPU
                inffo = ArrayFire_OpenCL_device_info();
                k2 = strfind(inffo, 'MB');
                if options.use_device == 0
                    k1 = strfind(inffo, '[0]');
                    dispaus = [dispaus, ' with ', inffo(k1 + 4:k2(1)+1)];
                else
                    k1 = strfind(inffo, ['-' num2str(options.use_device) '-']);
                    if isempty(k1)
                        k1 = strfind(inffo, ['[' num2str(options.use_device) ']']);
                    end
                    dispaus = [dispaus, ' with ', inffo(k1 + 4:k2(options.use_device + 1)+1)];
                end
            end
        elseif options.implementation == 3
            inffo = OpenCL_device_info();
            k6 = strfind(inffo, 'No. platforms: ');
            k2 = strfind(inffo, 'MB');
            k1 = strfind(inffo, 'Platform') + 12;
            k3 = strfind(inffo, 'Device');
            k4 = strfind(inffo, 'No. devices:');
            n_platforms = str2double(inffo(k6+15:k1(1)-15));
            n_devices = [];
            hh = 1;
            for kk = 1 : n_platforms
                n_devices = [n_devices ; str2double(inffo(k4(kk)+14:k3(hh)-2))];
                hh = hh + n_devices(kk);
            end
            k5 = strfind(inffo, 'with ');
            MB = [];
            hh = 1;
            for kk = 1 : n_platforms
                if kk == options.use_device + 1
                    for ll = 1 : n_devices(kk)
                        MB = [MB; str2double(inffo(k5(hh + ll - 1)+5:k2(hh + ll - 1)-2))];
                    end
                end
                hh = hh + n_devices(kk);
            end
            kumulatiivinen = cumsum(n_devices);
            if options.cpu_to_gpu_factor == 0
                [~, loc] = max(MB);
                ll = kumulatiivinen(options.use_device + 1) - n_devices(options.use_device + 1) + loc;
                if n_devices(options.use_device + 1) > 1
                    dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                        sprintf('\n'), inffo(k3(ll):k3(ll + 1) - 2)];
                else
                    if options.use_device + 1 == n_platforms
                        dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                            sprintf('\n'), inffo(k3(kumulatiivinen(end)):end - 1)];
                    else
                        dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                            sprintf('\n'), inffo(k3(kumulatiivinen(options.use_device + 1)):k1(options.use_device + 2) - 15)];
                    end
                end
            else
                loc = MB > 2048;
                if n_devices(options.use_device + 1) > 1
                    dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2)];
                    for hh = 1 : n_devices(options.use_device + 1)
                        if loc(hh) == 0
                            continue;
                        end
                        ll = kumulatiivinen(options.use_device + 1) - n_devices(options.use_device + 1) + hh;
                        dispaus = [dispaus, sprintf('\n'), inffo(k3(ll):k3(ll + 1) - 2)];
                    end
                else
                    if options.use_device + 1 == n_platforms
                        dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                            sprintf('\n'), inffo(k3(kumulatiivinen(end)):end - 1)];
                    else
                        dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                            sprintf('\n'), inffo(k3(kumulatiivinen(options.use_device + 1)):k1(options.use_device + 2) - 15)];
                    end
                end
            end
        end
        disp(dispaus);
        reko = {};
        varNonMAP = [recNames(6)];
        for kk = 1 : numel(varNonMAP)
            if options.OSEM && options.subsets == 1
                reko = [reko;{'MLEM'}];
                break
            elseif options.(varNonMAP{kk})
                ch = strrep(varNonMAP{kk},'_','-');
                reko = [reko;ch];
                break
            end
        end

        varMAP = recNames(2);
        for kk = 1 : numel(varMAP)
            if options.(varMAP{kk}) && strcmp(varMAP{kk},'OSL_OSEM') && PRIOR && options.subsets == 1
                reko = [reko;{'MLEM-OSL'}];
            elseif options.(varMAP{kk}) && strcmp(varMAP{kk},'OSL_COSEM')
                if options.OSL_COSEM == 1 && PRIOR
                    reko = [reko;{'ACOSEM-OSL'}];
                elseif options.OSL_COSEM == 2 && PRIOR
                    reko = [reko;{'COSEM-OSL'}];
                end
            elseif options.(varMAP{kk}) && PRIOR
                ch = strrep(varMAP{kk},'_','-');
                reko = [reko;ch];
            elseif options.(varMAP{kk}) && ~PRIOR
                ch = strrep(varMAP{kk},'_','-');
                reko = [reko;ch];
            end
        end
        for kk = 1 : length(reko)
            if kk == 1
                dispi = reko{kk};
            else
                dispi = [dispi, ', ' reko{kk}];
            end
            dispi = [dispi, ' reconstruction method selected.'];
        end
        if ~isempty(reko)
            disp(dispi)
        end
        priori = {};

        varPrior = recNames(1);
        varPriorName = recNames(10);
        for kk = 1 : numel(varPrior)
            if options.(varPrior{kk}) && MAP
                if strcmp(varPrior{kk},'TV')
                    if options.TV_use_anatomical
                        if options.TVtype == 1
                            priori = [priori;{'Anatomically weighted TV'}];
                        elseif options.TVtype == 2
                            priori = [priori;{'Joint TV'}];
                        elseif options.TVtype == 3
                            priori = [priori;{'Weighted joint TV'}];
                        end
                    else
                        if options.TVtype == 1 || options.TVtype == 2
                            priori = [priori;{'TV'}];
                        elseif options.TVtype == 3
                            priori = [priori;{'Weighted TV'}];
                        elseif options.TVtype == 4
                            priori = [priori;{'SATV'}];
                        end
                    end
                elseif strcmp(varPrior{kk},'weighted_mean')
                    if options.mean_type == 1 || options.mean_type == 4
                        priori = [priori;{'Weighted (arithmetic) mean'}];
                    elseif options.mean_type == 2 || options.mean_type == 5
                        priori = [priori;{'Weighted (harmonic) mean'}];
                    elseif options.mean_type == 3 || options.mean_type == 6
                        priori = [priori;{'Weighted (geometric) mean'}];
                    end
                elseif strcmp(varPrior{kk},'AD')
                    priori = [priori;{'AD-MRP'}];
                else
                    ch = strrep(varPriorName{kk},'_','-');
                    priori = [priori;ch];
                end
            end
        end
        for kk = 1 : length(priori)
            if kk == 1
                dispi2 = priori{kk};
            else
                dispi2 = [dispi2, ', ' priori{kk}];
            end
            if kk == length(priori) && length(priori) > 1
                dispi2 = [dispi2, ' priors selected.'];
            elseif kk == length(priori) && isscalar(priori)
                if ~options.NLM
                    dispi2 = [dispi2, ' prior selected.'];
                else
                    if options.NLTV
                        dispi2 = [dispi2, ' prior selected with NLTV.'];
                    elseif options.NLRD 
                        dispi2 = [dispi2, ' prior selected with NLRD.'];
                    elseif options.NLLange
                        dispi2 = [dispi2, ' prior selected with NL Lange.'];
                    elseif options.NLGGMRF
                        dispi2 = [dispi2, ' prior selected with NLGGMRF.'];
                    elseif options.NLGM
                        dispi2 = [dispi2, ' prior selected with NL Geman-McClure.'];
                    elseif options.NLM_MRP
                        dispi2 = [dispi2, ' prior selected with filtering mode.'];
                    else
                        dispi2 = [dispi2, ' prior selected.'];
                    end
                    if options.NLAdaptive
                        dispi2 = [dispi2, ' Using adaptive weighting.'];
                    end
                    if options.NLMaxWeight
                        dispi2 = [dispi2, ' Using the maximum weight for the reference voxel.'];
                    end
                end
            end
        end
        if ~isempty(priori)
            disp(dispi2)
            if options.NLM || options.MRP || options.quad || options.Huber || (options.RDP && options.RDPIncludeCorners) || options.GGMRF
                disp(['Using neighborhood/search window size of ' num2str(options.Ndx * 2 + 1) 'x' num2str(options.Ndy * 2 + 1) 'x' num2str(options.Ndz * 2 + 1) '.'])
            end
            if options.NLM
                disp(['Using patch size of ' num2str(options.Nlx * 2 + 1) 'x' num2str(options.Nly * 2 + 1) 'x' num2str(options.Nlz * 2 + 1) '.'])
            end
        end
        if options.projector_type == 1 && ~options.precompute_lor
            if options.implementation == 1
                disp('Improved Siddon''s algorithm selected with 1 ray.')
            else
                if options.n_rays_transaxial > 1
                    ray = 'rays';
                else
                    ray = 'ray';
                end
                if options.n_rays_axial > 1
                    aray = 'rays';
                else
                    aray = 'ray';
                end
                disp(['Improved Siddon''s algorithm selected with ' num2str(options.n_rays_transaxial) ' transaxial ' ray ' and ' ...
                    num2str(options.n_rays_axial) ' axial ' aray '.'])
            end
        elseif options.projector_type == 1 || options.projector_type == 11
            disp('Improved Siddon''s algorithm selected with 1 ray.')
        elseif options.projector_type == 2 || options.projector_type == 22
            dispi = 'Orthogonal distance-based ray tracer selected';
            if options.tube_width_z > 0
                dispi = [dispi, ' in 3D mode.'];
            else
                dispi = [dispi, ' in 2.5D mode.'];
            end
            disp(dispi)
        elseif options.projector_type == 3 || options.projector_type == 33
            disp('Volume of intersection based ray tracer selected.');
        elseif options.projector_type == 4 || options.projector_type == 44
            disp('Interpolation-based projector selected.')
        elseif options.projector_type == 5 || options.projector_type == 55
            disp('Branchless distance-driven based projector selected.')
        elseif options.projector_type == 6 || options.projector_type == 66
            disp('Rotation-based projector selected (SPECT).')
        elseif options.projector_type > 10
            fpType = floor(options.projector_type / 10);
            bpType = options.projector_type - fpType * 10;
            dispi = [];
            if fpType == 1
                dispi = [dispi, 'Improved Siddon''s algorithm selected for forward projection, '];
            elseif fpType == 2
                dispi = [dispi, 'Orthogonal distance-based ray-tracer selected for forward projection, '];
            elseif fpType == 3
                dispi = [dispi, 'Volume of intersection based ray tracer selected for forward projection, '];
            elseif fpType == 4
                dispi = [dispi, 'Interpolation-based projector selected for forward projection, '];
            elseif fpType == 5
                dispi = [dispi, 'Branchless distance-driven projector selected for forward projection, '];
            end
            if bpType == 1
                dispi = [dispi, 'improved Siddon''s algorithm for backprojection.'];
            elseif bpType == 2
                dispi = [dispi, 'orthogonal distance-based ray-tracer for backprojection.'];
            elseif bpType == 3
                dispi = [dispi, 'volume of intersection based ray tracer for backprojection.'];
            elseif bpType == 4
                dispi = [dispi, 'interpolation-based projector for backprojection.'];
            elseif bpType == 5
                dispi = [dispi, 'branchless distance-driven projector for backprojection.'];
            end
            disp(dispi)
        % elseif options.projector_type == 21
        %     disp('Improved Siddon''s algorithm selected for forward projection, orthogonal for backprojection.')
        % elseif options.projector_type == 12
        %     disp('Orthogonal selected for forward projection, improved Siddon''s algorithm for backprojection.')
        % elseif options.projector_type == 31
        %     disp('Improved Siddon''s algorithm selected for forward projection, Volume of intersection based ray tracer for backprojection.')
        % elseif options.projector_type == 13
        %     disp('Volume of intersection based ray tracer selected for forward projection, improved Siddon''s algorithm for backprojection.')
        % elseif options.projector_type == 41
        %     disp('Interpolation-based projector selected for forward projection, improved Siddon for backprojection.')
        % elseif options.projector_type == 14
        %     disp('Improved Siddon projector selected for forward projection, interpolation-based projector for backprojection.')
        % elseif options.projector_type == 15
        %     disp('Improved Siddon projector selected for forward projection, branchless distance-driven projector for backprojection.')
        % elseif options.projector_type == 45
        %     disp('Interpolation-based projector selected for forward projection, branchless distance-driven projector for backprojection.')
        % elseif options.projector_type == 54
        %     disp('Branchless distance-driven projector selected for forward projection, interpolation-based projector for backprojection.')
        % elseif options.projector_type == 51
        %     disp('Branchless distance-driven projector selected for forward projection, improved Siddon for backprojection.')
        end
        if options.use_psf
            if options.deblurring
                disp('PSF ON with deblurring phase.');
            else
                disp('PSF ON.');
            end

        end
        if options.attenuation_correction && ~options.CT
            if options.CT_attenuation
                disp('Attenuation correction (image domain) ON.')
            else
                disp('Attenuation correction (measurement domain) ON.')
            end
        end
        if options.randoms_correction && ~options.CT
            dispi = 'Randoms correction ON';
            if options.variance_reduction
                dispi = [dispi, ' with variance reduction'];
                if options.randoms_smoothing
                    dispi = [dispi, ' and smoothing'];
                end
            elseif options.randoms_smoothing
                dispi = [dispi, ' with smoothing'];
            end
            dispi = [dispi, '.'];
            disp(dispi)
        end
        % The non-subtracted scatter correction is disabled with list-mode data (see projectorClass)
        if options.scatter_correction && ~((options.useIndexBasedReconstruction || options.listmode || listModeCandidate) && ~options.CT && ~options.SPECT && ~options.subtract_scatter)
            dispi = 'Scatter correction ON';
            if options.scatter_smoothing
                dispi = [dispi, ' with smoothing'];
            end
            dispi = [dispi, '.'];
            disp(dispi)
        end
        if options.fill_sinogram_gaps && ~options.CT
            disp('Sinogram gap filling ON.')
        end
        if any(options.precondTypeImage)
            if options.precondTypeImage(1)
                disp('Using image-based preconditioning with diagonal preconditioner.')
            elseif options.precondTypeImage(2)
                disp('Using image-based preconditioning with EM preconditioner.')
            elseif options.precondTypeImage(3)
                disp('Using image-based preconditioning with IEM preconditioner.')
            end
            if options.precondTypeImage(4)
                disp('Using image-based preconditioning with momentum preconditioner.')
            end
            if options.precondTypeImage(5)
                disp('Using image-based preconditioning with normalized gradient preconditioner.')
            end
            if options.precondTypeImage(6)
                dispP = ('Using image-based preconditioning with filtering preconditioner');
                if strcmp(options.filterWindow,'hamming')
                    dispP = [dispP, ' with Hamming window.'];
                elseif strcmp(options.filterWindow,'hann')
                    dispP = [dispP, ' with Hann window.'];
                elseif strcmp(options.filterWindow,'blackman')
                    dispP = [dispP, ' with Blackman window.'];
                elseif strcmp(options.filterWindow,'nuttal')
                    dispP = [dispP, ' with Nuttal window.'];
                elseif strcmp(options.filterWindow,'gaussian')
                    dispP = [dispP, ' with Gaussian window.'];
                elseif strcmp(options.filterWindow,'shepp-logan')
                    dispP = [dispP, ' with Shepp-Logan window.'];
                elseif strcmp(options.filterWindow,'cosine')
                    dispP = [dispP, ' with cosine window.'];
                elseif strcmp(options.filterWindow,'parzen')
                    dispP = [dispP, ' with Parzen window.'];
                else
                    dispP = [dispP, ' (no windowing).'];
                end
                disp(dispP)
            end
        end
        if any(options.precondTypeMeas)
            if options.precondTypeMeas(1)
                disp('Using measurement-based preconditioning with diagonal preconditioner.')
            end
            if options.precondTypeMeas(2)
                dispP = ('Using measurement-based preconditioning with filtering preconditioner');
                if strcmp(options.filterWindow,'hamming')
                    dispP = [dispP, ' with Hamming window.'];
                elseif strcmp(options.filterWindow,'hann')
                    dispP = [dispP, ' with Hann window.'];
                elseif strcmp(options.filterWindow,'blackman')
                    dispP = [dispP, ' with Blackman window.'];
                elseif strcmp(options.filterWindow,'nuttal')
                    dispP = [dispP, ' with Nuttal window.'];
                elseif strcmp(options.filterWindow,'gaussian')
                    dispP = [dispP, ' with Gaussian window.'];
                elseif strcmp(options.filterWindow,'shepp-logan')
                    dispP = [dispP, ' with Shepp-Logan window.'];
                elseif strcmp(options.filterWindow,'cosine')
                    dispP = [dispP, ' with cosine window.'];
                elseif strcmp(options.filterWindow,'parzen')
                    dispP = [dispP, ' with Parzen window.'];
                else
                    dispP = [dispP, ' (no windowing).'];
                end
                disp(dispP)
            end
        end
        if options.oOffsetZ ~= 0 || options.oOffsetX ~= 0 || options.oOffsetY ~= 0
            disp(['Object offset is [' num2str(options.oOffsetX) ', ' num2str(options.oOffsetY) ', ' num2str(options.oOffsetZ) '] (XYZ).'])
        end
    end
    if options.normalization_correction && ~(options.compute_normalization || options.only_sinos) && ~options.CT
        disp('Normalization correction ON.')
    elseif options.compute_normalization && ~options.CT
        disp('Computing normalization coefficients.')
    end
    if options.useMultiResolutionVolumes
        if options.transaxialEFOV
            Nx = options.NxOrig + round((options.Nx - options.NxOrig) * options.multiResolutionScale);
            Ny = options.NyOrig + round((options.Ny - options.NyOrig) * options.multiResolutionScale);
        else
            Nx = options.Nx;
            Ny = options.Ny;
        end
        if options.axialEFOV
            Nz = options.NzOrig + round((options.Nz - options.NzOrig) * options.multiResolutionScale);
        else
            Nz = options.Nz;
        end
    else
        Nx = options.Nx;
        Ny = options.Ny;
        Nz = options.Nz;
    end
    if ~(options.compute_normalization || options.only_sinos)
        if options.corrections_during_reconstruction && (options.normalization_correction || options.randoms_correction || options.scatter_correction)
            disp('Corrections applied during reconstruction (ordinary Poisson).')
        elseif ~options.corrections_during_reconstruction && (options.normalization_correction || options.randoms_correction || options.scatter_correction)
            disp('Corrections applied to the measurement data.')
        end
        if options.arc_correction && ~options.use_raw_data && ~options.CT
            disp('Arc correction ON.')
        end
        if options.use_raw_data && ~options.CT
            if isempty(partitions) || partitions == 1 
                dispi = 'Using STATIC raw data';
            else
                dispi = 'Using DYNAMIC raw data';
            end
            if options.reconstruct_trues
                dispi = strcat(dispi, ' (trues)');
            elseif options.reconstruct_scatter
                dispi = strcat(dispi, ' (scatter)');
            else
                dispi = strcat(dispi, ' (prompts)');
            end
            if options.sampling_raw > 1
                dispi = strcat(dispi, [' with ' num2str(options.sampling_raw) 'x sampling']);
            end
            if partitions > 1
                if options.sampling_raw > 1
                    dispi = strcat(dispi, [' and ' num2str(partitions) ' time steps']);
                else
                    dispi = strcat(dispi, [' with ' num2str(partitions) ' time steps']);
                end
            end
            dispi = strcat(dispi, '.');
            disp(dispi)
        else
            % Type of the data: projection images (CT), list-mode events
            % (custom coordinates), index-based list-mode events (detector
            % indices) or sinograms
            if options.CT
                dataType = 'projection data';
            elseif options.useIndexBasedReconstruction
                dataType = 'index-based list-mode data';
            elseif options.listmode || listModeCandidate
                dataType = 'list-mode data';
            else
                dataType = 'sinogram data';
            end
            if isempty(partitions) || partitions == 1
                dispi = ['Using STATIC ' dataType];
            else
                dispi = ['Using DYNAMIC ' dataType];
            end
            if options.reconstruct_trues
                dispi = strcat(dispi, ' (trues)');
            elseif options.reconstruct_scatter
                dispi = strcat(dispi, ' (scatter)');
            else
                dispi = strcat(dispi, ' (prompts)');
            end
            if options.sampling > 1
                dispi = strcat(dispi, [' with ' num2str(options.sampling) 'x sampling']);
            end
            if numel(options.partitions) > 1
                if options.sampling > 1
                    dispi = strcat(dispi, [' and ' num2str(numel(options.partitions)) ' time steps']);
                else
                    dispi = strcat(dispi, [' with ' num2str(numel(options.partitions)) ' time steps']);
                end
            elseif numel(options.partitions) == 1 && options.partitions > 1
                if options.sampling > 1
                    dispi = strcat(dispi, [' and ' num2str(options.partitions) ' time steps']);
                else
                    dispi = strcat(dispi, [' with ' num2str(options.partitions) ' time steps']);
                end
            end
            dispi = strcat(dispi, '.');
            disp(dispi)
        end
        if options.subsets > 1
            if ismember(mod(options.subsets, 100), [11, 12, 13])
                abbr = 'th';
            elseif mod(options.subsets, 10) == 1
                abbr = 'st';
            elseif mod(options.subsets, 10) == 2
                abbr = 'nd';
            elseif mod(options.subsets, 10) == 3
                abbr = 'rd';
            else
                abbr = 'th';
            end
            if options.subset_type == 1
                disp(['Every ' num2str(options.subsets) abbr ' column measurement is taken per subset.'])
            elseif options.subset_type == 2
                disp(['Every ' num2str(options.subsets) abbr ' row measurement is taken per subset.'])
            elseif options.subset_type == 3
                disp('Using random subset sampling.')
            elseif options.subset_type == 4
                disp(['Every ' num2str(options.subsets) abbr ' sinogram column is taken per subset.'])
            elseif options.subset_type == 5
                disp(['Every ' num2str(options.subsets) abbr ' sinogram row is taken per subset.'])
            elseif options.subset_type == 6
                disp(['Using angle-based subset sampling with ' num2str(options.n_angles) ' angles combined per subset.'])
            elseif options.subset_type == 7
                disp('Using golden angle-based subset sampling.')
            elseif options.subset_type == 8
                disp(['Using every ' num2str(options.subsets) abbr ' sinogram/projection image.'])
            elseif options.subset_type == 9
                disp('Using sinograms/projection images in random order.')
            elseif options.subset_type == 10
                disp('Using golde angle sampling with sinograms/projections.')
            elseif options.subset_type == 11
                disp('Using prime factor ordering of projections/sinograms into subsets.')
            elseif options.subset_type == 0
                disp(['Dividing data into ' num2str(options.subsets) ' segments.'])
            end
        end
        disp(['Using an image (matrix) size of ' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) ' with ' num2str(options.Niter) ...
            ' iterations and ' num2str(options.subsets) ' subsets.'])
    elseif options.CT
        disp(['Using an image (matrix) size of ' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) ' with ' num2str(options.Niter) ...
            ' iterations and ' num2str(options.subsets) ' subsets.'])
    end
end
