function options = tileCTNormalization(options)
%TILECTNORMALIZATION Checks the size of the CT normalization and tiles it to all projections
%   The CT normalization (options.normalization) is a multiplicative factor of each measurement. It can be input either as
%   one map of the detector (options.nRowsD * options.nColsD elements), which is used for every projection, or as one value
%   for each measurement (options.nRowsD * options.nColsD * options.nProjections elements). The output is always a column
%   vector with one value for each measurement. Any other number of elements is an error. Both the input and the output
%   are the param struct of the projectorClass object.
nDet = double(options.nRowsD) * double(options.nColsD);
nProj = double(options.nProjections);
nNorm = numel(options.normalization);
if nNorm == nDet * nProj
    options.normalization = options.normalization(:);
elseif nNorm == nDet
    options.normalization = repmat(options.normalization(:), nProj, 1);
else
    error(['Incorrect size for the CT normalization (options.normalization)! Must be the size of the detector [' num2str(options.nRowsD) ' ' ...
        num2str(options.nColsD) '] (used for every projection) or the size of all projections [' num2str(options.nRowsD) ' ' ...
        num2str(options.nColsD) ' ' num2str(options.nProjections) '], but it has ' num2str(nNorm) ' elements.'])
end
end
