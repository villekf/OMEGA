function coincidences = initial_michelogram(options, coincidences)
%% Form initial raw Michelograms
% Forms the initial raw Michelograms from the raw list-mode data. I.e. a
% Michelogram without any spanning yet applied. Used by form_sinograms in
% forming the actual sinograms. 
%
% OUTPUT:
%   coincidences = A cell matrix containing the raw Michelogram for each
%   time step. The cell matrix has a size of rings x rings, where each cell
%   element contains the the coincidences between the current cell indices.
%   E.g. coincidences{2}(1,5) contains the coincidences at time step 2
%   between the detectors located on rings 1 and 5.
%
% See also form_sinograms, load_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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


folder = fileparts(which('initial_michelogram.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');
pseudot = options.pseudot;
machine_name = options.machine_name;


temp = pseudot;
if ~isempty(temp) && sum(temp) > 0
    for kk = 1 : temp
        pseudot(kk) = (options.cryst_per_block + 1) * kk;
    end
end

% Forms a raw Michelogram that is still in the detector space
if exist([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'], 'file') == 2
    load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'swap1', 'swap2', 'swap3', 'swap4');
else
    sinogram_coordinates_2D(options);
    load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'swap1', 'swap2', 'swap3', 'swap4');
end
% koko2 = options.det_per_ring * options.rings;
% koko1 = koko2;
% koko = options.det_per_ring^2/2 + options.det_per_ring/2;
% tri = tril(true(koko2,koko2), 0);
koko = options.det_per_ring*options.rings;
for llo=1:options.partitions
    %     tic
    
    prompt = coincidences{llo};
    [K, ~, V] = find(prompt);
    L = find(tril(true(koko,koko), 0));
    L = L(K);
    
    [I,J] = ind2sub([koko koko], L);
    clear L
    prompt = sparse(I,J,V,koko,koko);
    clear I J V
    
    eka=1;
    ll = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P1 = cell(options.rings + length(pseudot),(options.rings + length(pseudot) - eka + 1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     alku = 1;
%     loppu = options.det_per_ring;
%     alku2 = options.det_per_ring + 1;
%     loppu2 = options.det_per_ring * 2;
    for ring = eka:options.rings % Ring at a time
        for luuppi = 1:(options.rings-eka+1)-(ring-1) % All the LORs with other options.rings as well
            if luuppi == 1 % Both of the detectors are on the same ring
                % If pseudo ring then place the measurements to the next
                % ring
                if ismember(ll,pseudot)
                    ll = ll + 1;
                end
                lk = ll + 1;
                % Observations equaling the ones detected on the current
                % ring
                temppiP = full(prompt((1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1)),1+options.det_per_ring*(ring-1):(options.det_per_ring+options.det_per_ring*(ring-1))));
%                 apu = tril(true(koko2,koko2), 0);
%                 apu((1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1)),1+options.det_per_ring*(ring-1):(options.det_per_ring+options.det_per_ring*(ring-1))) = false;
%                 apu = ~apu;
%                 apu = apu(tril(true(koko2,koko2), 0));
%                 temppiP = full(prompt(apu));
%                 temppiP = zeros(koko,1);
%                 for kk = 1 : options.det_per_ring
%                     temppiP(sum(1:options.det_per_ring) - sum(1:options.det_per_ring - (kk-1)) + 1 : sum(1:options.det_per_ring) - sum(1:options.det_per_ring - kk)) = prompt(alku : loppu);
%                     koko2 = koko2 - 1;
%                     alku = alku + koko2;
%                     loppu = loppu + koko2 - 1;
%                 end
%                 loppu = loppu + options.det_per_ring;
                % combine same LORs, but with different detector order
                % (i.e. combine values at [325 25100] and [25100 325])
                temppiP = full(temppiP(tril(logical(true(options.det_per_ring)),0))); % Finally take only the other side
                
                if options.partitions > 1
                    P1{ll,ll} = sparse(double(temppiP(:)));
                else
                    P1{ll,ll} = uint16(temppiP(:));
                end
                
            else % Detectors on different rings
                if ismember(lk,pseudot)
                    lk = lk + 1;
                end
%                 temppiP2 = zeros(options.det_per_ring,options.det_per_ring);
                temppiP2 = (prompt(1+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1):(options.det_per_ring+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1)),(1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1))));
%                 for kk = 1 : options.det_per_ring
%                     temppiP2(:,kk) = prompt(alku2 : loppu2);
%                     koko1 = koko1 - 1;
%                     alku2 = alku2 + koko1;
%                     loppu2 = loppu2 + koko1;
%                 end                
%                 apu = tril(true(koko2,koko2), 0);
%                 apu(1+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1):(options.det_per_ring+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1)),(1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1))) = false;
%                 apu = ~apu;
%                 apu = apu(tril(true(koko2,koko2), 0));
%                 temppiP2 = full(prompt(apu));
%                 temppiP2 = reshape(temppiP2, options.det_per_ring,options.det_per_ring);
                
                % Combine LORs
                temppiP = full(triu(temppiP2));
                temppiP2 = full(tril(temppiP2));
                
                temppiP3 = temppiP;
                
                % Swap corners
                temppiP4 = zeros(size(temppiP2),'double');
                temppiP4(swap3) = temppiP2(swap3);
                temppiP(swap1) = 0;
                temppiP = temppiP + temppiP4';
                temppiP4 = zeros(size(temppiP2),'double');
                temppiP4(swap4) = temppiP2(swap4);
                temppiP(swap2) = 0;
                temppiP = temppiP + temppiP4';
                temppiP = temppiP';
                
                temppiP4 = zeros(size(temppiP3),'double');
                temppiP4(swap1) = temppiP3(swap1);
                temppiP2(swap3) = 0;
                temppiP2 = temppiP2 + temppiP4';
                temppiP4 = zeros(size(temppiP3),'double');
                temppiP4(swap2) = temppiP3(swap2);
                temppiP2(swap4) = 0;
                temppiP2 = temppiP2 + temppiP4';
                
                % Take only the other side
                temppiP2 = temppiP2(tril(true(options.det_per_ring),0));
                temppiP = temppiP(tril(true(options.det_per_ring),0));
                
                if options.partitions > 1
                    P1{ll,lk} = sparse(temppiP(:));
                    P1{lk,ll} = sparse(temppiP2(:));
                else
                    P1{ll,lk} = uint16(temppiP(:));
                    P1{lk,ll} = uint16(temppiP2(:));
                end
                lk = lk + 1;
            end
        end
        ll = ll + 1;
    end
    clear prompt
    if options.partitions > 1
        P1(cellfun(@isempty,P1)) = {sparse(size(temppiP))};
    else
        P1(cellfun(@isempty,P1)) = {zeros(size(temppiP),'uint16')};
    end
    coincidences{llo} = P1;
    
end