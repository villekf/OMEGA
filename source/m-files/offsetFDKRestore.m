function options = offsetFDKRestore(options, state)
%OFFSETFDKRESTORE Puts back the values that offsetFDKPreprocess changed
%   Restores nRowsD (and the quantities that depend on it), the detector
%   coordinates, offsetCorrection and the measurement data of the
%   projectorClass object. Does nothing if state is empty.
%
% Inputs:
%   options = The projectorClass object
%   state = The second output of offsetFDKPreprocess
% Output:
%   options = The projectorClass object with the original values

if isempty(state)
    return
end
fn = fieldnames(state.param);
for kk = 1 : numel(fn)
    options.param.(fn{kk}) = state.param.(fn{kk});
end
options.x = state.x;
options.nMeasSubset = state.nMeasSubset;
end
