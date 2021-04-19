function options = loadDelayedData(options)
%LOAD DELAYED COINCIDENCE DATA Loads the randoms correction data
%   Loads randoms correction data from a mat file. If the mat-file contains
%   'SinDelayed' (sinogram data) or 'delayed_coincidences' (raw data)
%   variable then that will be automatically used. Otherwise the first
%   variable is loaded.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021 Ville-Veikko Wettenhovi
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


[d_file, d_fpath] = uigetfile('*.mat','Select delayed coincidence data');
if isequal(d_file, 0)
    error('No file was selected')
end

FileName = fullfile(d_fpath, d_file);
storedStructure = load(FileName);
variables = fieldnames(storedStructure);

if isfield(storedStructure, 'raw_SinDelayed') && ~options.use_raw_data
    options.SinDelayed = storedStructure.raw_SinDelayed;
elseif isfield(storedStructure, 'SinDelayed') && ~options.use_raw_data
    options.SinDelayed = storedStructure.SinDelayed;
elseif isfield(storedStructure, 'delayed_coincidences') && options.use_raw_data
    options.SinDelayed = storedStructure.delayed_coincidences;
else
    options.SinDelayed = storedStructure.(variables{1});
end
clear scatter_file s_fpath FileName storedStructure variables
disp('Delayed coincidence data successfully loaded')
end

