function options = SPECTParameters(options)
%SPECTParameters Computes the necessary variables for projector type 6
%(SPECT)
%   Computes the PSF standard deviations for the current SPECT collimator
%   if omitted. Also computes the interaction planes with the rotated
%   projector and the projection angles.
%
% Copyright (C) 2022 Ville-Veikko Wettenhovi, Matti Kortelainen
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
DistanceToFirstRow = options.radiusPerProj-(options.Nx/2-0.5)*options.dx;
Distances = repmat(DistanceToFirstRow,1,options.Nx)+repmat([[0:options.Nx-1]*options.dx],length(DistanceToFirstRow),1);
Distances = Distances-options.colL-options.colD; %these are distances to the actual detector surface

if ~isfield(options,'gFilter')
    if ~isfield(options, 'sigmaZ')
        Rg = 2*options.colR*(options.colL+options.colD+(Distances)+options.cr_p/2)/options.colL; %Anger, "Scintillation Camera with Multichannel Collimators", J Nucl Med 5:515-531 (1964)
        Rg(Rg<0) = 0;

        FWHM = sqrt(Rg.^2+options.iR^2);
        FWHM_pixel = FWHM/options.dx;
        % expr = FWHM_pixel.^2-FWHMrot^2;
        % expr(expr<=0) = 10^-16;
        % FWHM_WithinPlane = sqrt(expr);

        %Parametrit CDR-mallinnukseen
        options.sigmaZ = FWHM_pixel./(2*sqrt(2*log(2)));
        % options.sigmaXY = FWHM_WithinPlane./(2*sqrt(2*log(2)));
    end
    y = options.Nx / 2 - 0.5:-1:-options.Nx / 2 + 0.5;
    x = options.Nz / 2 - 0.5:-1:-options.Nz / 2 + 0.5;
    xx = repmat(x', 1,options.Nz);
    yy = repmat(y, size(xx,1),1);
    if ~isfield(options, 'sigmaXY')
        s1 = repmat(permute(options.sigmaZ.^2,[4 3 2 1]), size(xx,1), size(yy,2), 1);
        options.gFilter = single(1 / (2*pi*s1).*exp(-(xx.^2 + yy.^2)./(2*s1)));
    else
        s1 = repmat(permute(options.sigmaZ,[4 3 2 1]), size(xx,1), size(yy,2), 1);
        s2 = repmat(permute(options.sigmaXY,[4 3 2 1]), size(xx,1), size(yy,2), 1);
        options.gFilter = single(1 / (2*pi*(s1).*(s2)).*exp(-(xx.^2./(2*s1.^2) + yy.^2./(2*s2.^2))));
    end
    ind = max(find(options.radiusPerProj == max(options.radiusPerProj),1,'first'));
    rowE = find(options.gFilter(round(size(options.gFilter,1)/2),:,end,ind) > 1e-8,1,'last');
    rowS = find(options.gFilter(round(size(options.gFilter,1)/2),:,end,ind) > 1e-8,1,'first');
    colE = find(options.gFilter(:,round(size(options.gFilter,2)/2),end,ind) > 1e-8,1,'last');
    colS = find(options.gFilter(:,round(size(options.gFilter,2)/2),end,ind) > 1e-8,1,'first');
    options.gFilter = options.gFilter(rowS:rowE,colS:colE,:,:);
end
[~, options.blurPlanes] = max(Distances>0,[],2);
options.blurPlanes = uint32(options.blurPlanes - 1);
options.angles = single(repelem(options.startAngle, options.nProjections / 2) + repmat((0:options.angleIncrement:options.angleIncrement * (options.nProjections / 2 - 1))', options.nHeads, 1) + options.offangle);
if max(abs(options.angles)) > 8 * pi
    options.angles = options.angles / 180 * pi;
end
if options.flip_image
    options.angles = -(options.angles);
end
end