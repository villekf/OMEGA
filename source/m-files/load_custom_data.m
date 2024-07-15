function options = load_custom_data(options)
%LOADCUSTOMDATA Loads measurement data from the selected mat-file
% Supports either raw list-mode data or sinogram data. For both cases
% either trues or prompts can be used. The input data in case of raw data
% needs to be named either true_coincidences (trues) or coincidences
% (prompts). For the sinogram case, SinTrues and SinM, respectively.

if options.use_raw_data
    [options.file, options.fpath] = uigetfile('*.mat','Select raw list-mode data');
    if isequal(options.file, 0)
        error('No file was selected')
    end
    if options.reconstruct_trues
        load(options.file,'true_coincidences')
        options.SinM = true_coincidences;
    elseif options.reconstruct_scatter
        load(options.file,'scattered_coincidences')
        options.SinM = scattered_coincidences;
    else
        load(options.file,'coincidences')
        options.SinM = coincidences;
    end
else
    [options.file, options.fpath] = uigetfile('*.mat','Select sinogram data');
    if isequal(options.file, 0)
        error('No file was selected')
    end
    if options.reconstruct_trues
        load(options.file,'SinTrues')
        options.SinM = SinTrues;
    elseif options.reconstruct_scatter
        load(options.file,'SinScatter')
        options.SinM = SinScatter;
    else
        if options.corrections_during_reconstruction || (~options.normalization_correction && ~options.randoms_correction && ~options.scatter_correction)
            load(options.file,'raw_SinM')
            options.SinM = raw_SinM;
        else
            load(options.file,'SinM')
            options.SinM = SinM;
        end
    end
end