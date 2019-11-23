function [Sin, gaps] = gapFilling(options, Sin, xp, yp, llo, gaps)
%GAPFILLING Performs gap filling to pseudo detector sinograms
%   Dedicated function for sinogram gap filling. If pseudo detectors are
%   used (e.g. Biograph mCT), these detectors don't actually receive any
%   measurements. This causes zero measurements when no mashing and reduced
%   number of measurements with mashing. The gap filling interpolation
%   attempts to fix these gaps.
% 
%   Two types of interpolation are available. The built-in fillmissing or
%   the function inpaint_nans from File exchange: 
%   https://se.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
%
% INPUTS:
%   options = Machine and sinogram properties
%   Sin = The original sinogram
%   xp, yp = Sinogram coordinates
%   llo = Current timestep
%   gaps = Location of the gaps to be filled (used for subsequent
%   timesteps)
%
% OUTPUTS:
%   Sin = The gap-filled sinogram
%   gaps = The location of the sinogram gaps for subsequent timesteps
%
% See also fillmissing, inpaint_nans, form_sinograms


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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

mashing = options.det_w_pseudo / options.Nang / 2;
if mashing == 1
    apu = sum(Sin,3);
    gaps = apu == 0;
else
    if llo == 1
        options.Nang = options.Nang * mashing;
        pseudot = options.pseudot;
        ringsp = options.rings + sum(pseudot);
        temp = pseudot;
        if ~isempty(temp) && temp > 0
            for kk = uint32(1) : temp
                pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
            end
        elseif temp == 0
            pseudot = [];
        end
        [~, ~, i, j, accepted_lors] = sinogram_coordinates_2D(options, xp, yp);
        Sinog = cell(ringsp,ringsp);
        
        ix=cellfun('isempty',Sinog);
        Sinog(ix)={zeros(options.Ndist,options.Nang,'uint16')};
        
        P1 = options.coincidences{llo};
        
        % Create the Michelograms
        for ii=1:ringsp
            if any(ii==pseudot)
                continue
            end
            for jj=1:ringsp
                if any(jj==pseudot)
                    continue
                end
                if issparse(P1{ii,jj})
                    CC = uint16(full(P1{ii,jj}));
                else
                    CC = P1{ii,jj};
                end
                CC = CC(accepted_lors);
                %                 Sinog{ii,jj}(ind) = CC;
                Sinog{ii,jj} = uint16(accumarray([i j],CC,[options.Ndist options.Nang]));
            end
        end
        Sinog = cat(3,Sinog{:});
        Sin_orig = zeros(options.Ndist,options.Nang,options.TotSinos,'uint16');
        
        
        kkj = zeros(floor((options.ring_difference-ceil(options.span/2))/options.span) + 1, 1);
        for kk = 1 : floor((options.ring_difference-ceil(options.span/2))/options.span) + 1
            kkj(kk) = ceil(options.span/2) + options.span*(kk - 1);
        end
        offset2 = cumsum(options.segment_table);
        % Create the sinograms
        % First the detectors on the same ring
        Sin_orig(:,:,1:2:options.Nz) = Sinog(:,:,1:ringsp+1:ringsp^2);
        % Then the detectors on adjacent rings
        for jh=1:floor(options.span/2)
            apu = Sinog(:,:,jh*ringsp+1:ringsp+1:ringsp^2);
            apu2 = Sinog(:,:,jh+1:ringsp+1:(ringsp-jh)*ringsp);
            Sin_orig(:,:,jh+1:2:offset2(1)-jh) = Sin_orig(:,:,jh+1:2:offset2(1)-jh) + apu + apu2;
        end
        % Lastly the rest of the detectors with the amount of combined LORs
        % specified with the span value
        for ih=1:floor(length(options.segment_table)/2)
            for jh=1:span
                apu = Sinog(:,:,(kkj(ih)+jh-1)*ringsp+1:ringsp+1:end);
                Sin_orig(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = Sin_orig(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
                apu2 = Sinog(:,:,kkj(ih)+jh:ringsp+1:(ringsp-kkj(ih)-jh+1)*ringsp);
                Sin_orig(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = Sin_orig(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
            end
        end
        clear Sinog
        apu = sum(Sin_orig,3);
        apu(apu == 0) = NaN;
        apu(~isnan(apu)) = 0;
        apu(isnan(apu)) = 1;
        for kk = 1 : mashing - 1
            apu = apu(:,1:2:end) + apu(:,2:2:end);
        end
        gaps = apu > 0;
    end
end
if strcmp('fillmissing',options.gap_filling_method)
    Sin = single(Sin);
    for kk = 1 : size(Sin,3)
        apu = Sin(:,:,kk);
        apu(gaps) = NaN;
        Sin(:,:,kk) = fillmissing(apu, options.interpolation_method_fillmissing);
    end
elseif strcmp('inpaint_nans',options.gap_filling_method)
    Sin = single(Sin);
    for kk = 1 : size(Sin,3)
        apu = Sin(:,:,kk);
        apu(gaps) = NaN;
        Sin(:,:,kk) = inpaint_nans(apu, options.interpolation_method_inpaint);
    end
else
    warning('Unsupported gap filling method! No gap filling was performed!')
end
if options.verbose
    disp('Gap filling completed')
end
end

