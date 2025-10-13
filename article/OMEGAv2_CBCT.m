% Used data available from: https://doi.org/10.5281/zenodo.12722386
% The data has to be in MATLAB/Octave path!

clear

% type 0 = PDHG with filtering
% type 1 = PDHG with filtering and NLRDP regularization
% type 2 = PKMA with RDP regularization
type = 0;

% The selected GPU device
GPUDevice = 0;

% If your GPU doesn't have enough memory to contain all the data set the
% below variable to true. It divides the image into segments that are
% reconstructed separately. This slows down the computations, but allows
% them to be computed even on GPUs with very limited amount of memory
% 4 GB of memory is recommended even when this is set as true, 12 GB
% otherwise
largeDim = false;
[z, t, recPar] = CBCT_example(type, GPUDevice, largeDim);