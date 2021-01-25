function [options,Prop, input] = applyScatterRandoms(options,Prop, input, type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if type == 0
    ichar = 'SinDelayed';
    corr = 'randoms';
elseif type == 1
    ichar = 'ScatterC';
    corr = 'scatter';
end
if isfield(options,ichar)
    if iscell(input)
        if iscell(options.SinM)
            for kk = 1 : length(options.SinM)
                if numel(input{kk}) ~= numel(options.SinM{kk}) && numel(input{kk}) * options.TOF_bins ~= numel(options.SinM{kk})
                    error(['Size mismatch between ' corr ' correction data and measurement data'])
                end
                if options.variance_reduction && ~Prop.variance_reduction
                    input{kk} = Randoms_variance_reduction(single(input{kk}), options);
                    Prop.variance_reduction = true;
                end
                if options.randoms_smoothing && ~Prop.smoothing
                    input{kk} = randoms_smoothing(double(input{kk}), options);
                    Prop.smoothing = true;
                end
                if ~options.corrections_during_reconstruction
                    if type == 0
                        if options.TOF_bins > 1
                            NN = options.Ndist * options.Nang * options.NSinos;
                            options.SinM{kk} = options.SinM{kk}(:);
                            for uu = 1 : options.TOF_bins
                                options.SinM{kk}(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM{kk}(NN * (uu - 1) + 1: NN * uu))) - single(full(input{kk}(:)) / options.TOF_bins);
                            end
                        else
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) - single(full(input{kk}));
                        end
                    elseif type == 1
                        if options.subtract_scatter
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) - single(full(input{kk}(:)));
                        else
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) .* single(full(input{kk}(:)));
                        end
                    end
                    options.SinM{kk}(options.SinM{kk} < 0) = 0;
                end
            end
        else
            if numel(input{1}) ~= numel(options.SinM) && numel(input{1}) * options.TOF_bins ~= numel(options.SinM)
                error(['Size mismatch between ' corr ' correction data and measurement data'])
            end
            if options.variance_reduction && ~Prop.variance_reduction
                input{1} = Randoms_variance_reduction(single(input{1}), options);
                Prop.variance_reduction = true;
            end
            if options.randoms_smoothing && ~Prop.smoothing
                input{1} = randoms_smoothing(double(input{1}), options);
                Prop.smoothing = true;
            end
            if ~options.corrections_during_reconstruction
                if type == 0
                    if options.TOF_bins > 1
                        NN = options.Ndist * options.Nang * options.NSinos;
                        options.SinM = options.SinM(:);
                        for uu = 1 : options.TOF_bins
                            options.SinM(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM(NN * (uu - 1) + 1: NN * uu))) - single(full(input{1}(:)) / options.TOF_bins);
                        end
                    else
                        options.SinM = single(full(options.SinM(:))) - single(full(input{1}));
                    end
                elseif type == 1
                    if options.subtract_scatter
                        options.SinM = single(full(options.SinM(:))) - single(full(input{1}(:)));
                    else
                        options.SinM = single(full(options.SinM(:))) .* single(full(input{1}(:)));
                    end
                end
                options.SinM(options.SinM < 0) = 0;
            end
        end
    else
        if iscell(options.SinM)
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
                if ~options.corrections_during_reconstruction
                    if type == 0
                        if options.TOF_bins > 1
                            NN = options.Ndist * options.Nang * options.NSinos;
                            options.SinM{kk} = options.SinM{kk}(:);
                            for uu = 1 : options.TOF_bins
                                options.SinM{kk}(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM{kk}(NN * (uu - 1) + 1: NN * uu))) - single(full(input(:)) / options.TOF_bins);
                            end
                        else
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) - single(full(input(:)) / length(options.SinM));
                        end
                    elseif type == 1
                        if options.subtract_scatter
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) - single(full(input(:)) / length(options.SinM));
                        else
                            options.SinM{kk} = single(full(options.SinM{kk}(:))) .* single(full(input(:)) / length(options.SinM));
                        end
                    end
                    options.SinM{kk}(options.SinM{kk} < 0) = 0;
                end
            end
        else
            if numel(input) ~= numel(options.SinM) && numel(input) * options.TOF_bins ~= numel(options.SinM)
                error(['Size mismatch between ' corr ' correction data and measurement data'])
            end
            if options.variance_reduction && ~Prop.variance_reduction
                input = Randoms_variance_reduction(single(input), options);
                Prop.variance_reduction = true;
            end
            if options.randoms_smoothing && ~Prop.smoothing
                input = randoms_smoothing(input, options);
                Prop.smoothing = true;
            end
            if ~options.corrections_during_reconstruction
                if type == 0
                    if options.TOF_bins > 1
                        NN = options.Ndist * options.Nang * options.NSinos;
                        options.SinM = options.SinM(:);
                        for uu = 1 : options.TOF_bins
                            options.SinM(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM(NN * (uu - 1) + 1: NN * uu))) - single(full(input(:)) / options.TOF_bins);
                        end
                    else
                        options.SinM = single(full(options.SinM(:))) - single(full(input(:)));
                    end
                elseif type == 1
                    if options.subtract_scatter
                        options.SinM = single(full(options.SinM(:))) - single(full(input(:)));
                    else
                        options.SinM = single(full(options.SinM(:))) .* single(full(input(:)));
                    end
                end
                options.SinM(options.SinM < 0) = 0;
            end
        end
    end
else
    disp('Delayed coincidences not found, randoms correction not performed')
    input = 0;
    options.randoms_correction = false;
end
end

