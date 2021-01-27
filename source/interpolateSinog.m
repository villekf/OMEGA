function SinM = interpolateSinog(SinM,sampling, Ndist, partitions, sampling_interpolation_method)
%INTERPOLATESINOG Interpolates additional rows to the sinogram
%   Additional rows are determined by input value sampling. The
%   interpolation method can be controlled with
%   sampling_interpolation_method.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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

Nang = size(SinM,2);
[X,Y] = meshgrid(1:Nang, 1:sampling:Ndist*sampling+1);
[Xq,Yq] = meshgrid(1:Nang, 1:Ndist*sampling+1);
if iscell(SinM)
    SinM_uus = cell(size(SinM));
    for hh = 1 : partitions
        SinM_uus{hh} = zeros(size(SinM{hh},1)*sampling, size(SinM{hh},2),size(SinM{hh},3),size(SinM{hh},4),'single');
        for uu = 1 : size(SinM{hh},4)
            for kk = 1 : size(SinM{hh},3)
                if strcmp(sampling_interpolation_method,'spline') || strcmp(sampling_interpolation_method, 'makima')
                    SinM_uus{hh}(:,:,kk,uu) = interp2(X(1:end-1,:), Y(1:end-1,:),double(SinM{hh}(:,:,kk,uu)), Xq(1:end-1,:), Yq(1:end-1,:), sampling_interpolation_method);
                else
                    apu = interp2(X, Y,double([SinM{hh}(:,:,kk,uu);SinM{hh}(end-1,:,kk,uu)]), Xq, Yq, sampling_interpolation_method);
                    SinM_uus{hh}(:,:,kk,uu) = apu(1:end-1,:);
                end
            end
        end
        SinM_uus{hh} = SinM_uus{hh}(:);
    end
    SinM = SinM_uus;
else
    SinM_uus = zeros(size(SinM,1)*sampling, size(SinM,2),size(SinM,3),size(SinM,4),'single');
    for uu = 1 : size(SinM,4)
        for kk = 1 : size(SinM,3)
        if strcmp(sampling_interpolation_method,'spline') || strcmp(sampling_interpolation_method, 'makima')
            SinM_uus(:,:,kk,uu) = interp2(X(1:end-1,:), Y(1:end-1,:) ,double(SinM(:,:,kk,uu)), Xq(1:end-1,:), Yq(1:end-1,:), sampling_interpolation_method);
        else
            apu = interp2(X, Y,double([SinM(:,:,kk,uu);SinM(end-1,:,kk,uu)]), Xq, Yq, sampling_interpolation_method);
            SinM_uus(:,:,kk,uu) = apu(1:end-1,:);
        end
        end
    end
    SinM = single(SinM_uus(:));
end
end

