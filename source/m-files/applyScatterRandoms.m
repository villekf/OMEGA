function [options,Prop, input] = applyScatterRandoms(options,Prop, input, type)
%APPLYSCATTERRANDOMS Applies the selected corrections to randoms/scatter
%data and precorrects the measurement data if applicable
if type == 0
    ichar = 'SinDelayed';
    corr = 'randoms';
elseif type == 1
    ichar = 'ScatterC';
    corr = 'scatter';
end
if isfield(options,ichar)
    if iscell(input)
        if isfield(options,'SinM') && iscell(options.SinM)
            for kk = 1 : length(options.SinM)
                if numel(input{kk}) ~= numel(options.SinM{kk}) && numel(input{kk}) * options.TOF_bins ~= numel(options.SinM{kk})
                    error(['Size mismatch between ' corr ' correction data and measurement data'])
                end
                checkTOFRandomsFilter(input{kk}, options, Prop, type);
                if options.variance_reduction && ~Prop.variance_reduction
                    input{kk} = Randoms_variance_reduction(single(input{kk}), options);
                    Prop.variance_reduction = true;
                end
                if options.randoms_smoothing && ~Prop.smoothing
                    input{kk} = randoms_smoothing(double(input{kk}), options);
                    Prop.smoothing = true;
                end
                if (~options.corrections_during_reconstruction || ~options.ordinaryPoisson)
                    if type == 0
                        if options.TOF_bins > 1 && numel(input{kk}) ~= numel(options.SinM{kk})
                            NN = options.Ndist * options.Nang * options.NSinos;
                            options.SinM{kk} = single(full(options.SinM{kk}(:)));
                            for uu = 1 : options.TOF_bins
                                options.SinM{kk}(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM{kk}(NN * (uu - 1) + 1: NN * uu))) - single(full(input{kk}(:))) / options.TOF_bins;
                            end
                        else
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) - single(full(input{kk}(:)));
                        end
                    elseif type == 1
                        if options.subtract_scatter && options.TOF_bins > 1 && numel(input{kk}) * options.TOF_bins == numel(options.SinM{kk})
                            % Non-TOF scatter data with TOF data
                            options.SinM{kk} = subtractScatterTOF(options.SinM{kk}, input{kk}, options.TOF_bins, kk == 1);
                        elseif options.subtract_scatter
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) - single(full(input{kk}(:)));
                        else
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) .* single(full(input{kk}(:)));
                        end
                    end
                end
            end
        else
            if isfield(options,'SinM') && numel(input{1}) ~= numel(options.SinM) && numel(input{1}) * options.TOF_bins ~= numel(options.SinM)
                error(['Size mismatch between ' corr ' correction data and measurement data'])
            end
            checkTOFRandomsFilter(input{1}, options, Prop, type);
            if options.variance_reduction && ~Prop.variance_reduction
                input{1} = Randoms_variance_reduction(single(input{1}), options);
                Prop.variance_reduction = true;
            end
            if options.randoms_smoothing && ~Prop.smoothing
                input{1} = randoms_smoothing(double(input{1}), options);
                Prop.smoothing = true;
            end
            if isfield(options,'SinM') && (~options.corrections_during_reconstruction || ~options.ordinaryPoisson)
                if type == 0
                    if options.TOF_bins > 1 && numel(input{1}) ~= numel(options.SinM)
                        NN = options.Ndist * options.Nang * options.NSinos;
                        options.SinM = single(full(options.SinM(:)));
                        for uu = 1 : options.TOF_bins
                            options.SinM(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM(NN * (uu - 1) + 1: NN * uu))) - single(full(input{1}(:))) / options.TOF_bins;
                        end
                    else
                        options.SinM = single(full(options.SinM(:))) - single(full(input{1}(:)));
                    end
                elseif type == 1
                    if options.subtract_scatter && options.TOF_bins > 1 && numel(input{1}) * options.TOF_bins == numel(options.SinM)
                        % Non-TOF scatter data with TOF data
                        options.SinM = subtractScatterTOF(options.SinM, input{1}, options.TOF_bins, true);
                    elseif options.subtract_scatter
                        options.SinM = single(full(options.SinM(:))) - single(full(input{1}(:)));
                    else
                        options.SinM = single(full(options.SinM(:))) .* single(full(input{1}(:)));
                    end
                end
            end
        end
    else
        if isfield(options,'SinM') && iscell(options.SinM)
            checkTOFRandomsFilter(input, options, Prop, type);
            if options.variance_reduction && ~Prop.variance_reduction
                input = Randoms_variance_reduction(single(input), options);
                Prop.variance_reduction = true;
            end
            if options.randoms_smoothing && ~Prop.smoothing
                input = randoms_smoothing(single(input), options);
                Prop.smoothing = true;
            end
            for kk = 1 : length(options.SinM)
                if numel(input) ~= numel(options.SinM{kk}) && numel(input) * options.TOF_bins ~= numel(options.SinM{kk})
                    error(['Size mismatch between ' corr ' correction data and measurement data'])
                end
                if (~options.corrections_during_reconstruction || ~options.ordinaryPoisson)
                    if type == 0
                        if options.TOF_bins > 1 && numel(input) ~= numel(options.SinM{kk})
                            NN = options.Ndist * options.Nang * options.NSinos;
                            options.SinM{kk} = single(full(options.SinM{kk}(:)));
                            for uu = 1 : options.TOF_bins
                                options.SinM{kk}(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM{kk}(NN * (uu - 1) + 1: NN * uu))) - single(full(input(:))) / options.TOF_bins;
                            end
                        else
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) - single(full(input(:))) / length(options.SinM);
                        end
                    elseif type == 1
                        if options.subtract_scatter && options.TOF_bins > 1 && numel(input) * options.TOF_bins == numel(options.SinM{kk})
                            % Non-TOF scatter data with TOF data
                            options.SinM{kk} = subtractScatterTOF(options.SinM{kk}, single(full(input(:))) / length(options.SinM), options.TOF_bins, kk == 1);
                        elseif options.subtract_scatter
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) - single(full(input(:))) / length(options.SinM);
                        else
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) .* (single(full(input(:))) / length(options.SinM));
                        end
                    end
                end
            end
        else
            if isfield(options,'SinM') && numel(input) ~= numel(options.SinM) && numel(input) * options.TOF_bins ~= numel(options.SinM)
                error(['Size mismatch between ' corr ' correction data and measurement data'])
            end
            checkTOFRandomsFilter(input, options, Prop, type);
            if options.variance_reduction && ~Prop.variance_reduction
                input = Randoms_variance_reduction(single(input), options);
                Prop.variance_reduction = true;
            end
            if options.randoms_smoothing && ~Prop.smoothing
                input = randoms_smoothing(input, options);
                Prop.smoothing = true;
            end
            if isfield(options,'SinM') && (~options.corrections_during_reconstruction || ~options.ordinaryPoisson)
                if type == 0
                    if options.TOF_bins > 1 && numel(input) ~= numel(options.SinM)
                        NN = options.Ndist * options.Nang * options.NSinos;
                        options.SinM = single(full(options.SinM(:)));
                        for uu = 1 : options.TOF_bins
                            options.SinM(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM(NN * (uu - 1) + 1: NN * uu))) - single(full(input(:))) / options.TOF_bins;
                        end
                    else
                        options.SinM = single(full(options.SinM(:))) - single(full(input(:)));
                    end
                elseif type == 1
                    if options.subtract_scatter && options.TOF_bins > 1 && numel(input) * options.TOF_bins == numel(options.SinM)
                        % Non-TOF scatter data with TOF data
                        options.SinM = subtractScatterTOF(options.SinM, input, options.TOF_bins, true);
                    elseif options.subtract_scatter
                        options.SinM = single(full(options.SinM(:))) - single(full(input(:)));
                    else
                        options.SinM = single(full(options.SinM(:))) .* single(full(input(:)));
                    end
                end
            end
        end
    end
else
    disp('Delayed coincidences not found, randoms correction not performed')
    input = 0;
    options.randoms_correction = false;
end
end

function sinM = subtractScatterTOF(sinM, scatter, tofBins, showMessage)
% Subtracts the non-TOF scatter data from the TOF measurement data: the scatter data is divided evenly among the TOF bins
% (TOF is the slowest dimension of the measurement data). Randoms are handled the same way above.
if showMessage
    disp(['The scatter data has no TOF bins and is divided evenly among the ' num2str(tofBins) ' TOF bins'])
end
sinM = single(full(sinM(:)));
scatter = single(full(scatter(:)));
NN = numel(scatter);
for uu = 1 : tofBins
    sinM(NN * (uu - 1) + 1: NN * uu) = sinM(NN * (uu - 1) + 1: NN * uu) - scatter / tofBins;
end
end

function checkTOFRandomsFilter(randoms, options, Prop, type)
% TOF-resolved randoms (one value for each TOF bin, as many elements as the TOF measurements) are subtracted as they are in each TOF
% bin. The variance reduction and the smoothing of the randoms work on a non-TOF sinogram, so they are not supported with them
if type == 0 && options.TOF_bins > 1 && isfield(options,'SinM') && ...
        ((options.variance_reduction && ~Prop.variance_reduction) || (options.randoms_smoothing && ~Prop.smoothing))
    sinM = options.SinM;
    if iscell(sinM)
        sinM = sinM{1};
    end
    if numel(randoms) == numel(sinM)
        error('Randoms variance reduction and randoms smoothing are not supported with TOF-resolved randoms data; disable them or use a non-TOF randoms sinogram')
    end
end
end
