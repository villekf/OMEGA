function options = gatherHostFields(options, keepFields)
%GATHERHOSTFIELDS Moves every gpuArray field of a struct back to the host
%   Fields whose names are listed in keepFields are left on the device. Used
%   before calling the MEX-files, which read every field except the designated
%   device inputs with the regular (host) mxArray API.
%
% Examples:
%   options = gatherHostFields(options, {'x0'});
%   options = gatherHostFields(options, {});
%
% INPUTS:
%   options = The options/param struct that is passed to the projector
%   MEX-files
%   keepFields = Cell array of field name strings that should be left
%   untouched, i.e. kept on the device if they currently are a gpuArray
%
% OUTPUTS:
%   options = The input struct with every gpuArray field, except those
%   listed in keepFields, replaced by its gathered (host) counterpart
%
% See also forwardProjection, backwardProjection

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

fn = fieldnames(options);
for kk = 1 : numel(fn)
    if ismember(fn{kk}, keepFields)
        continue
    end
    apu = options.(fn{kk});
    if isa(apu, 'gpuArray')
        options.(fn{kk}) = gather(apu);
    elseif iscell(apu)
        for ll = 1 : numel(apu)
            if isa(apu{ll}, 'gpuArray')
                apu{ll} = gather(apu{ll});
            end
        end
        options.(fn{kk}) = apu;
    end
end
end
