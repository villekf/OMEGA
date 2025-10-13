% Used data available from: https://doi.org/10.5281/zenodo.17185907
% The data has to be in MATLAB/Octave path!

clear

% The bin file can either be in the current working directory or
% alternatively you can specify the directory below:
path = '';

% type 0 = OSEM with improved version of Siddon ray-based projector
% type 1 = OSEM with volume of intersection ray-based projector
% type 2 = list-mode PET, PKMA with relative difference prior regularization and with improved version of Siddon ray-based projector
type = 0;

% The selected GPU device
GPUDevice = 0;

[z, t, recPar] = PET_TOF_example(type, GPUDevice, path);