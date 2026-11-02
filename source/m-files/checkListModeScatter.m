function checkListModeScatter(options)
%CHECKLISTMODESCATTER Validates the scatter of list-mode and index-based PET data
%
% With list-mode and index-based PET data the scatter (options.ScatterC with
% options.scatter_correction = true and options.subtract_scatter = true) is
% the additive term of each event: one finite, non-negative value for each
% event, in the order of the events (options.trIndex and options.axIndex, or
% the custom coordinates) before the subsets. The unit is the measured counts
% expected for the LOR of the event, in the TOF bin of the event with TOF
% data, over the frame (the additive term s of the model n * (A x) + s), with
% no division by the normalization, the attenuation, the number of TOF bins or
% the number of LORs in a bin. A randoms estimate is added into the same
% vector by the user.
%
% Called from projectorClass once the number of events is known. Does
% nothing if the data is not list-mode PET data or the scatter correction
% does not subtract the scatter.
%
% INPUTS:
%   options = Struct containing the reconstruction options
%
% See also checkIndexBasedInputs, loadCorrections

if ~options.listmode || options.CT || options.SPECT || ~options.scatter_correction || ~options.subtract_scatter ...
        || options.reconstruct_trues || options.reconstruct_scatter
    return
end
if numel(options.partitions) > 1 || options.partitions > 1 || (isfield(options, 'Nt') && options.Nt > 1)
    error(['Scatter correction (options.ScatterC) is not supported with dynamic list-mode and index-based data! ' ...
        'Use a single time frame or disable the scatter correction.'])
end
if ~options.corrections_during_reconstruction || ~options.ordinaryPoisson
    error(['Scatter correction of list-mode and index-based data is the additive term of each event and needs corrections during the ' ...
        'reconstruction (options.corrections_during_reconstruction = true and ordinary Poisson algorithm, i.e. not PDHG, FISTA, LSQR, CGLS, ...); ' ...
        'subtracting the scatter from the events (precorrection) is not possible!'])
end
if options.largeDim
    error('Scatter correction of list-mode and index-based data is not supported with largeDim!')
end
if options.scatter_smoothing
    error(['Scatter smoothing is a sinogram operation and cannot be used with list-mode and index-based data! ' ...
        'Smooth the scatter before inputting it as a value for each event.'])
end
if options.useIndexBasedReconstruction
    nEvents = numel(options.trIndex) / 2;
elseif iscell(options.SinM)
    nEvents = 0;
    for kk = 1 : numel(options.SinM)
        nEvents = nEvents + numel(options.SinM{kk});
    end
else
    nEvents = numel(options.SinM);
end
if options.TOF_bins > 1
    tofText = ' in the TOF bin of the event';
    divText = ', the number of TOF bins or the number of LORs in a bin';
else
    tofText = '';
    divText = ' or the number of LORs in a bin';
end
units = ['one value per event, in the order of the events (options.trIndex and options.axIndex, or the custom coordinates), ' ...
    'as the measured counts expected for the LOR of the event' tofText ' over the frame, i.e. the additive term s of the model ' ...
    'n * (A x) + s, with no division by the normalization, the attenuation' divText];
if ~isfield(options, 'ScatterC') || isempty(options.ScatterC)
    error(['Scatter correction (options.scatter_correction) of list-mode and index-based data needs options.ScatterC, but it is empty! ' ...
        'It must contain ' units '. It has to contain ' num2str(nEvents) ' elements.'])
end
if iscell(options.ScatterC)
    error(['For list-mode and index-based data options.ScatterC must be a single array of ' num2str(nEvents) ' elements: ' units '!'])
end
if numel(options.ScatterC) ~= nEvents
    error(['The scatter (options.ScatterC) has ' num2str(numel(options.ScatterC)) ' elements, but the number of events is ' num2str(nEvents) ...
        '! It must contain ' units '!'])
end
if ~isnumeric(options.ScatterC) && ~islogical(options.ScatterC)
    error(['The scatter (options.ScatterC) must be a numeric array, but its class is ' class(options.ScatterC) '!'])
end
scatterData = full(options.ScatterC(:));
if any(~isfinite(scatterData))
    error(['The scatter (options.ScatterC) contains NaN or infinite values! It must contain ' units '!'])
end
if min(scatterData) < 0
    error(['The scatter (options.ScatterC) contains negative values (minimum ' num2str(double(min(scatterData))) ')! It must contain ' units '!'])
end
