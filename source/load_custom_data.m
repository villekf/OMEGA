function options = load_custom_data(options)


if options.use_raw_data
    [options.file, options.fpath] = uigetfile('*.mat','Select raw list-mode data');
    if options.reconstruct_trues
        load(options.file,'true_coincidences')
        options.SinM = true_coincidences;
    else
        load(options.file,'coincidences')
        options.SinM = coincidences;
    end
else
    [options.file, options.fpath] = uigetfile('*.mat','Select sinogram data');
    if options.reconstruct_trues
        load(options.file,'SinTrues')
        options.SinM = SinTrues;
    else
        load(options.file,'SinM')
        options.SinM = SinM;
    end
end