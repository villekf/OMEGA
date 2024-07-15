function [pz,varargout] = reconstructions_mainCT(options)
%RECONSTRUCTIONS_MAINCT Utility function for CT reconstruction
%   This function simply adds PET specific variables that are needed for
%   error-free functionality and also converts some input variables to the
%   correct format (e.g. the projection angles are converted to radians if
%   they are input as degrees).
%
%   pz = The reconstructed image. Can be a 4D matrix if N number of
%   iterations are saved of if a dynamic reconstruction is performed.
%   image_properties = Reconstruction specific variables that are initially
%   stored in the options-struct (optional)
%   options = The modified input options-struct (optional)

options.CT = true;
options.errorChecking = true;
if nargout > 1
    if nargout == 3
        [pz,varargout{1},varargout{2}] = reconstructions_main(options);
    elseif nargout == 4
        [pz,varargout{1},varargout{2},varargout{3}] = reconstructions_main(options);
    elseif nargout == 5
        [pz,varargout{1},varargout{2},varargout{3},varargout{4}] = reconstructions_main(options);
    else
        [pz,varargout{1}] = reconstructions_main(options);
    end
else
    pz = reconstructions_main(options);
end
end

