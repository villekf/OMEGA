function SinM = resizeType6Sinogram(SinM, endSinogramRows, endSinogramCols, sinogramSizeRows, sinogramSizeCols)
%RESIZETYPE6SINOGRAM Resamples a SPECT detector image onto the voxel lattice of the rotation-based projector (type 6)
%   The detector image (rows x columns x projections [x frames], pixel pitch dPitchX x dPitchY) is padded with zeros or trimmed
%   symmetrically to the size of the FOV in detector pixels, endSinogramRows = FOVa_x / dPitchX and endSinogramCols = axial_fov /
%   dPitchY, and then resampled with imresize to the sinogramSizeRows x sinogramSizeCols voxels of the image. After that the rows and
%   columns correspond to the voxel lattice of the image and are centered on the center of the detector.
%
%   The resampling keeps the VALUES of the detector pixels (imresize is an interpolation of the pixel values, it does not preserve
%   the sum of the counts), i.e. the value of a voxel-pitch pixel is the (interpolated) value of the detector pixels at the same
%   physical position.
%
%   resize needs an integer number of pixels. Pad/trim + imresize is used as it always was when endSinogramRows/Cols is exactly an
%   integer and the difference of the sizes of the detector and the FOV is even (so the results of the geometries that work are not
%   changed). When the difference is odd the symmetric padding/trimming of resize is impossible: the extra pixel is put at the end, i.e.
%   the image is shifted by half a detector pixel. Such an axis is resampled in one step like a ratio that is not an integer.
%   A ratio that is not exactly an integer (the FOV is not a whole number of detector pixels, or only up to rounding error, e.g.
%   242 voxels * 2.46 mm / 4.92 mm) is resampled in one step directly from the detector pixels onto the voxel lattice, see
%   resampleDetectorAxis: this needs no fraction of a pixel of padding, has no half-pixel shift and is the same interpolation as
%   imresize (bicubic, with antialiasing when shrinking, values and not counts preserved), with zero outside of the detector as in the
%   padding. For a ratio that is an integer up to rounding error (relative 1e-9) the ratio is rounded first.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021-2025 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sinogramSizeRows = double(sinogramSizeRows);
sinogramSizeCols = double(sinogramSizeCols);
[endSinogramRows, exactRows] = nearestIntegerRatio(double(endSinogramRows));
[endSinogramCols, exactCols] = nearestIntegerRatio(double(endSinogramCols));

if exactRows && mod(size(SinM, 1) - endSinogramRows, 2) == 0
    SinM = resize(SinM, endSinogramRows, Dimension=1, FillValue=0, Side="both"); % Pad or trim sinogram rows (symmetric)
else
    SinM = resampleDetectorAxis(SinM, 1, sinogramSizeRows, endSinogramRows / sinogramSizeRows);
end
if exactCols && mod(size(SinM, 2) - endSinogramCols, 2) == 0
    SinM = resize(SinM, endSinogramCols, Dimension=2, FillValue=0, Side="both"); % Pad or trim sinogram columns (symmetric)
else
    SinM = resampleDetectorAxis(SinM, 2, sinogramSizeCols, endSinogramCols / sinogramSizeCols);
end
% Now the sinogram and FOV XZ-plane match in physical dimensions but not in resolution.
if size(SinM, 1) ~= sinogramSizeRows || size(SinM, 2) ~= sinogramSizeCols % Resampling to the same size is the identity
    SinM = imresize(SinM, [sinogramSizeRows, sinogramSizeCols]); % Resample the sinogram
end
end

function [n, isInteger] = nearestIntegerRatio(ratio)
% isInteger: the ratio is exactly an integer (pad/trim can be used if the difference to the detector size is even, see above). A
% ratio that is an integer up to rounding error (relative 1e-9) is rounded, but is resampled
n = round(ratio);
isInteger = ratio == n;
if ~isInteger && abs(ratio - n) > 1e-9 * max(1, abs(ratio))
    n = ratio;
end
end
