function [options, state] = offsetFDKPreprocess(options)
%OFFSETFDKPREPROCESS Offset-detector (half-fan) CBCT FDK preprocessing
%   Prepares the offset-detector CT data of an FDK reconstruction. The
%   redundancy weights are applied before the ramp filtering (weighting in
%   the backprojection, after the filtering, leaves a bright central
%   artifact) and the missing side of the detector is completed with
%   virtual detector rows by conjugate completion, see offsetFDKCore.
%   The measurement data, the detector centre (options.x) and the number of
%   detector rows (nRowsD) are changed accordingly and offsetCorrection is
%   switched off for the library call, since the weights have already been
%   applied.
%
%   Applies only when options.param.CT, options.param.FDK and
%   options.param.offsetCorrection are set. Call this after the projector
%   class object has been created and the data have been linearised, before
%   anything else uses nRowsD (parseInputData, prepass_phase). The
%   half-width of the transition region of the weights is
%   options.offsetFDKDelta (mm, on the detector).
%
% Input:
%   options = The projectorClass object
% Outputs:
%   options = The modified projectorClass object
%   state = The state needed by offsetFDKRestore to put back the original
%   values, empty if the preprocessing was not applied (the in-kernel
%   offset correction is then in effect)

state = [];
p = options.param;
if ~(p.CT && p.FDK && p.offsetCorrection) || p.listmode || p.useHelical
    return
end
N = double(p.nRowsD);
M = double(p.nColsD);
K = double(p.nProjections);
pitch = double(p.dPitchX);
if numel(p.partitions) > 1 || (~isempty(p.partitions) && p.partitions > 1)
    warning('Offset correction (FDK): the redundancy weighting with conjugate completion is not supported with dynamic data, using the in-kernel weights.')
    return
end
SinM = p.SinM;
isCell = iscell(SinM);
if isCell
    if numel(SinM) ~= 1
        warning('Offset correction (FDK): unsupported measurement data layout, using the in-kernel weights.')
        return
    end
    SinM = SinM{1};
end
x = options.x;
z = options.z;
if numel(x) ~= 6 * K || (numel(z) ~= 2 * K && numel(z) ~= 6 * K)
    warning('Offset correction (FDK): unsupported detector coordinates, using the in-kernel weights.')
    return
end
if numel(SinM) ~= N * M * K
    warning('Offset correction (FDK): the measurement data size does not match the detector, using the in-kernel weights.')
    return
end
if isfield(p, 'additionalCorrection') && p.additionalCorrection && isfield(p, 'corrVector') && numel(p.corrVector) > 1
    warning('Offset correction (FDK): corrections during the reconstruction are not supported with the redundancy weighting with conjugate completion, using the in-kernel weights.')
    return
end
deltaMm = double(p.offsetFDKDelta);
if ~isscalar(deltaMm) || ~(deltaMm > 0)
    error('options.offsetFDKDelta must be a positive scalar (mm)!')
end
x6 = reshape(double(x), 6, K);
if numel(z) == 6 * K
    u = reshape(double(z), 6, K);
    u = u(1:3, :);
else
    u = reshape(double(z), 2, K);
    u = [u; zeros(1, K)];
end
y = reshape(SinM, N, M, K);
largeDim = p.largeDim;
% The largeDim FDK linearises raw (not yet linearised) data in the library (log(flat / measurement)), so such data are linearised here for
% the preprocessing and converted back to the measurement domain below. Already linearised data (usingLinearizedData) stay as they are
if largeDim && ~p.usingLinearizedData
    y = log(single(p.flat) ./ single(y));
end
if isfield(p, 'offsetFDKPad')
    pad = p.offsetFDKPad;
else
    pad = [];
end
try
    [ypad, xNew, info] = offsetFDKCore(y, x6, u, pitch, deltaMm, pad);
catch ME
    if strncmp(ME.identifier, 'offsetFDK:', 10)
        warning('Offset correction (FDK): %s, using the in-kernel weights.', ME.message)
        return
    end
    rethrow(ME)
end
clear y
if info.raised
    warning('options.offsetFDKDelta (%g mm) is smaller than the physical overlap, using %g mm.', deltaMm, info.deltaMm)
end
if info.clamped
    warning('options.offsetFDKDelta is too large for the detector, using %g mm.', info.deltaMm)
end
% Back to the measurement domain for the library, unless the input was already linearised (the library then does not linearise it either)
if largeDim && ~p.usingLinearizedData
    ypad = single(p.flat) * exp(-ypad);
end
P = info.P;
newRows = N + P;
% Save the values that are changed
state.param = struct();
fn = {'SinM', 'nRowsD', 'Ndist', 'det_per_ring', 'det_w_pseudo', 'cryst_per_block', 'offsetCorrection', 'totMeas', 'size_x', 'dSizeXBP'};
for kk = 1 : numel(fn)
    if isfield(p, fn{kk})
        state.param.(fn{kk}) = p.(fn{kk});
    end
end
state.x = options.x;
state.nMeasSubset = options.nMeasSubset;
% Apply
if isCell
    p.SinM = {ypad};
else
    p.SinM = ypad;
end
p.nRowsD = newRows;
p.Ndist = newRows;
p.det_per_ring = newRows * K;
p.det_w_pseudo = p.det_per_ring;
p.cryst_per_block = M * newRows;
if isfield(p, 'totMeas')
    p.totMeas = M * newRows * K;
end
if isfield(p, 'size_x')
    p.size_x = cast(newRows, class(p.size_x));
end
if any(p.projector_type == [5 15 45 54 51]) && isfield(p, 'dSizeXBP')
    p.dSizeXBP = (newRows + 1) * p.dPitchY;
end
p.offsetCorrection = false;
options.param = p;
options.x = cast(xNew(:), class(x));
if p.subset_type >= 8 || p.subsets == 1
    kerroin = M * newRows;
else
    kerroin = 1;
end
for kk = 1 : p.subsets
    options.nMeasSubset(kk) = options.nMeas(kk + 1) * kerroin - options.nMeas(kk) * kerroin;
end
if p.verbose >= 1
    disp(['Offset correction (FDK): redundancy weights with conjugate completion (half-width ' num2str(info.deltaMm, '%.1f') ' mm), detector padded by ' num2str(P) ' rows'])
end
end
