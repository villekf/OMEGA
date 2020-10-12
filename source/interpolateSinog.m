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

if iscell(SinM)
    SinM_uus = cell(size(SinM));
    for hh = 1 : partitions
        SinM_uus{hh} = zeros(size(SinM{hh},1)*sampling, size(SinM{hh},2),size(SinM{hh},3),size(SinM{hh},4),'single');
        for uu = 1 : size(SinM{hh},4)
            for kk = 1 : size(SinM{hh},3)
                for ll = 1 : size(SinM{hh},2)
                    SinM_uus{hh}(:,ll,kk,uu) = single(interp1(1:sampling:Ndist*sampling+1, double([SinM{hh}(:,ll,kk,uu);SinM{hh}(end - 1,ll,kk,uu)]), ...
                        1:Ndist*sampling, sampling_interpolation_method));
                end
            end
        end
    end
    SinM = SinM_uus;
else
    SinM_uus = zeros(size(SinM,1)*sampling, size(SinM,2),size(SinM,3),size(SinM,4),'single');
    for uu = 1 : size(SinM,4)
        for kk = 1 : size(SinM,3)
            for ll = 1 : size(SinM,2)
                SinM_uus(:,ll,kk,uu) = interp1(1:sampling:Ndist*sampling+1, double([SinM(:,ll,kk,uu);SinM(end - 1,ll,kk,uu)]), ...
                    1:Ndist*sampling, sampling_interpolation_method);
            end
        end
    end
    SinM = single(SinM_uus);
end
end

