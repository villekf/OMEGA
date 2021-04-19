function im_vectors = reshape_vectors(im_vectors, options)
%RESHAPE_VECTORS This function simply reshapes the image vectors into
%actual images (matrices)

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
if options.save_iter
    Niter = options.Niter;
else
    Niter = 0;
end
fn = fieldnames(im_vectors);
fn = fn(cellfun('isempty',strfind(fn,'apu')));
for kk = 1 : numel(fn)
    im_vectors.(fn{kk}) = reshape(im_vectors.(fn{kk}),options.Nx,options.Ny,options.Nz,Niter+1);
end