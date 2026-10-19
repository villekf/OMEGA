function options = checkIndexBasedInputs(options)
%CHECKINDEXBASEDINPUTS Validates the inputs of index-based reconstruction
%
% Checks the shape, the data type and the number of events of
% options.trIndex and options.axIndex (and the number of options.TOFIndices
% when TOF is used) and casts the indices to uint16. Called from
% projectorClass once the data are available.
%
% Accepted layouts for options.trIndex and options.axIndex:
%   2xN matrix (or 2xNxNt) = the two detector indices of each event in the
%   two rows, i.e. the layout returned by load_data
%   Vector with 2N elements = the 2xN matrix in column-major order, i.e.
%   [det1_ev1, det2_ev1, det1_ev2, det2_ev2, ...]
%   Cell array = one of the above for each time step (dynamic data)
% Nx2 matrices are not accepted since they would be interpreted as 2N
% events. The values must be integers between 0 and 65535.
%
% INPUTS:
%   options = Struct containing the reconstruction options
%
% OUTPUTS:
%   options = The same struct with trIndex and axIndex as uint16
%
% See also load_data, form_subset_indices

if ~isfield(options, 'trIndex') || ~isfield(options, 'axIndex')
    error('Index-based reconstruction requires options.trIndex and options.axIndex, i.e. the transaxial and axial detector indices of the events (2 x number of events), but they are not set!')
end
if iscell(options.trIndex) ~= iscell(options.axIndex)
    error('options.trIndex and options.axIndex must either both be cell arrays (dynamic data, one cell per time step) or both be numeric arrays!')
end
if iscell(options.trIndex)
    if numel(options.trIndex) ~= numel(options.axIndex)
        error(['options.trIndex has ' num2str(numel(options.trIndex)) ' time steps, but options.axIndex has ' num2str(numel(options.axIndex)) '!'])
    end
    counts = zeros(numel(options.trIndex), 1);
    for kk = 1 : numel(options.trIndex)
        [options.trIndex{kk}, nTr] = normalizeDetectorIndices(options.trIndex{kk}, 'trIndex');
        [options.axIndex{kk}, nAx] = normalizeDetectorIndices(options.axIndex{kk}, 'axIndex');
        if nTr ~= nAx
            error(['options.trIndex and options.axIndex must contain the same number of events, but trIndex has ' num2str(nTr) ' and axIndex has ' num2str(nAx) ' (time step ' num2str(kk) ')!'])
        end
        counts(kk) = nTr;
    end
else
    [options.trIndex, nTr] = normalizeDetectorIndices(options.trIndex, 'trIndex');
    [options.axIndex, nAx] = normalizeDetectorIndices(options.axIndex, 'axIndex');
    if nTr ~= nAx
        error(['options.trIndex and options.axIndex must contain the same number of events, but trIndex has ' num2str(nTr) ' and axIndex has ' num2str(nAx) '!'])
    end
    counts = nTr;
end
if isfield(options, 'TOF') && options.TOF
    if isfield(options, 'TOFIndices')
        TOFIndices = options.TOFIndices;
    else
        TOFIndices = [];
    end
    if iscell(TOFIndices)
        tofCounts = cellfun(@numel, TOFIndices(:));
        if numel(counts) == 1
            tofCounts = sum(tofCounts);
        elseif numel(tofCounts) ~= numel(counts)
            error(['options.TOFIndices has ' num2str(numel(tofCounts)) ' time steps, but the index-based data has ' num2str(numel(counts)) '!'])
        end
    else
        tofCounts = numel(TOFIndices);
        counts = sum(counts);
    end
    for kk = 1 : numel(counts)
        if tofCounts(kk) ~= counts(kk)
            if numel(counts) > 1
                suffix = [' in time step ' num2str(kk)];
            else
                suffix = '';
            end
            error(['The number of TOF indices (options.TOFIndices, ' num2str(tofCounts(kk)) ') does not correspond to the number of events (' num2str(counts(kk)) ')' suffix '!'])
        end
    end
end
end

function [A, nEvents] = normalizeDetectorIndices(A, name)
% Validates one set of detector indices and returns it as uint16 together
% with the number of events
if ~isnumeric(A)
    error(['options.' name ' must contain integer detector indices, but its class is ' class(A) '!'])
end
if isempty(A)
    error(['options.' name ' is empty!'])
end
sz = size(A);
if ndims(A) == 2 && (sz(1) == 1 || sz(2) == 1)
    if mod(numel(A), 2) ~= 0
        error(['options.' name ' is a vector of ' num2str(numel(A)) ' elements, but it must contain two detector indices per event (an even number of elements, layout [det1_ev1, det2_ev1, det1_ev2, ...])! The expected size is 2 x N.'])
    end
    nEvents = numel(A) / 2;
elseif sz(1) == 2
    nEvents = prod(sz(2:end));
else
    error(['options.' name ' has the size ' strjoin(arrayfun(@num2str, sz, 'UniformOutput', false), 'x') ', but the expected size is 2 x N, i.e. the two detector indices of each of the N events in the two rows (alternatively a vector of length 2N with the layout [det1_ev1, det2_ev1, det1_ev2, ...]). Transpose the array if it is N x 2.'])
end
if ~isa(A, 'uint16')
    if ~isinteger(A)
        if any(~isfinite(A(:))) || any(A(:) ~= floor(A(:)))
            error(['options.' name ' must contain integer detector indices, but it contains non-integer values!'])
        end
    end
    if any(A(:) < 0)
        error(['options.' name ' contains negative detector indices (minimum ' num2str(double(min(A(:)))) ')!'])
    end
    if any(A(:) > intmax('uint16'))
        error(['options.' name ' contains detector indices that do not fit in uint16 (maximum ' num2str(double(max(A(:)))) ', limit 65535)!'])
    end
    A = uint16(A);
end
end
