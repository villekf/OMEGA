function options = loadScatterData(options)
%LOAD SCATTER DATA Loads the scatter correction data
%   Loads scatter correction data from a mat file. If the mat-file contains
%   'SinScatter' (sinogram data) or 'scattered_coincidences' (raw data)
%   variable then that will be automatically used. Otherwise the first
%   variable is loaded.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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


[scatter_file, s_fpath] = uigetfile('*.mat; *.scn','Select scatter correction data');
if isequal(scatter_file, 0)
    error('No file was selected')
end

if strcmp(scatter_file(end-3:end), '.scn')
    nimi = [s_fpath scatter_file];
    fid = fopen(nimi);
    options.ScatterC = fread(fid, inf, 'single=>single',0,'l');
    fclose(fid);
    s_length = length(options.ScatterC)/(options.Nang*options.Ndist*options.TotSinos);
    options.ScatterC = reshape(options.ScatterC,options.Ndist,options.Nang,options.TotSinos,s_length);
    if s_length > 1
        if options.partitions > 1
            temp = cell(s_length, 1);
            for kk = options.partitions : -1 : 1
                temp{kk} = options.ScatterC(:,:,:,kk);
                options.ScatterC = options.ScatterC(:,:,:,1:kk-1);
            end
        else
            temp = sum(options.ScatterC,4,'native');
            if ~options.subtract_scatter
                temp = temp / s_length;
            end
        end
        options.ScatterC = temp;
    end
else
    
    FileName = fullfile(s_fpath, scatter_file);
    storedStructure = load(FileName);
    variables = fields(storedStructure);
    
    if isfield(storedStructure, 'SinScatter') && ~options.use_raw_data
        options.ScatterC = storedStructure.SinScatter;
    elseif isfield(storedStructure, 'scattered_coincidences') && options.use_raw_data
        options.ScatterC = storedStructure.scattered_coincidences;
    else
        options.ScatterC = storedStructure.(variables{1});
        if length(variables) > 1
            warning('The scatter data file contains more than one variable! The first one is used')
        end
    end
    clear scatter_file s_fpath FileName storedStructure variables
    disp('Scatter data successfully loaded')
end

