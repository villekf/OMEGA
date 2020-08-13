function SinM = interpolateSinog(SinM,sampling, Ndist, partitions, sampling_interpolation_method)
%INTERPOLATESINOG Interpolates additional rows to the sinogram
%   Additional rows are determined by input value sampling. The
%   interpolation method can be controlled with
%   sampling_interpolation_method.

if iscell(SinM)
    SinM_uus = cell(size(SinM));
    for hh = 1 : partitions
        SinM_uus{hh} = zeros(size(SinM{hh},1)*sampling, size(SinM{hh},2),size(SinM{hh},3),size(SinM{hh},4),'single');
        for uu = 1 : size(SinM{hh},4)
            for kk = 1 : size(SinM{hh},3)
                for ll = 1 : size(SinM{hh},2)
                    SinM_uus{hh}(:,ll,kk,uu) = single(interp1(1:sampling:Ndist*sampling+1, double([SinM{hh}(:,ll,kk);SinM{hh}(end - 1,ll,kk)]), ...
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
                SinM_uus(:,ll,kk,uu) = interp1(1:sampling:Ndist*sampling+1, double([SinM(:,ll,kk);SinM(end - 1,ll,kk)]), ...
                    1:Ndist*sampling, sampling_interpolation_method);
            end
        end
    end
    SinM = single(SinM_uus);
end
end

