function im_vectors = form_image_vectors(options, N)
%FORM_IMAGE_VECTORS This function simply forms the applicable image vectors
%for each reconstruction method
%
%
%
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

if options.save_iter
    Niter = options.Niter;
else
    Niter = 0;
end
if options.implementation ~= 2
    type = 'double';
else
    type = 'single';
end


varList = recNames(3);
ll = 0;
kk = 1;
while ll == 0 && kk <= numel(varList)
    ll = ll + options.(varList{kk});
    kk = kk +1;
end
MLEM_bool = ll > 0;

varList = recNames(4);
ll = 0;
kk = 1;
while ll == 0 && kk <= numel(varList)
    ll = ll + options.(varList{kk});
    kk = kk +1;
end
OS_bool = ll > 0;
%     MLEM_bool = options.OSL_MLEM || options.mlem;
%     OS_bool = options.osem || options.rosem || options.ramla || options.OSL_OSEM || options.BSREM || options.ROSEM_MAP || options.rbi || options.drama ...
%         || options.cosem || options.ecosem || options.acosem || options.RBI_OSL || any(options.COSEM_OSL);

if options.implementation == 4 && OS_bool
    im_vectors.OSEM_apu = options.x0(:);
end
if options.implementation == 4 && MLEM_bool
    im_vectors.MLEM_apu = options.x0(:);
end
varList = ([recNames(5);recNames(6)]);

apu = false(numel(varList),1);
for kk = 1 : numel(varList)
    if options.(varList{kk})
        apu(kk) = true;
    end
end
varList = varList(apu);
varU = upper(varList);
varapu = strcat(varU,'_apu');

for kk = 1 : numel(varList)
    if options.(varList{kk})
        im_vectors.(varU{kk}) = ones(N,Niter + 1,type);
        im_vectors.(varU{kk})(:,1) = options.x0(:);
        im_vectors.(varapu{kk}) = im_vectors.(varU{kk})(:,1);
    end
end
varMAP = recNames(2);
varPrior = recNames(1);
apu = false(numel(varMAP),1);
for kk = 1 : numel(varMAP)
    if options.(varMAP{kk})
        apu(kk) = true;
    end
end
varMAP = varMAP(apu);
apu = false(numel(varPrior),1);
for kk = 1 : numel(varPrior)
    if options.(varPrior{kk})
        apu(kk) = true;
    end
end
varPrior = varPrior(apu);
nPrior = numel(varPrior);
varList = strcat(varPrior,'_');
varList = repmat(varList,numel(varMAP),1);
if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
    apu = repmat(varMAP,1,nPrior)';
    varList = strcat(varList,apu(:));
else
    varList = strcat(varList,repelem(varMAP,nPrior,1));
end
varapu = strcat(varList,'_apu');
uu = 1;
for kk = 1 : numel(varMAP)
    for ll = 1 : nPrior
        if options.(varMAP{kk}) && options.(varPrior{ll})
            im_vectors.(varList{uu}) = ones(N,Niter + 1,type);
            im_vectors.(varList{uu})(:,1) = options.x0(:);
            if options.implementation ~= 4
                im_vectors.(varapu{uu}) = im_vectors.(varList{uu})(:,1);
            end
        end
        uu = uu + 1;
    end
end
% Special ECOSEM case guarantees that OSEM and COSEM are initialized even
% if they haven't been selected
if options.ECOSEM
    if ~options.OSEM
        im_vectors.OSEM_apu = options.x0(:);
    end
    if ~options.COSEM
        im_vectors.COSEM_apu = options.x0(:);
    end
end
