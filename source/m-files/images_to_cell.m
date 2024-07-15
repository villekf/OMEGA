function [pz, varargout] = images_to_cell(im_vectors)
%IMAGES_TO_CELL Save the images at current time-step to a cell vector

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

% fn = fieldnames(im_vectors);
% fn = fn(cellfun('isempty',strfind(fn,'apu')));
gg = 1;
pz = im_vectors.recImage;
% loc = find(rekot);

% varNonMAP = [recNames(6)];
% apu = false(numel(varNonMAP),1);
% for kk = 1 : numel(varNonMAP)
%     if options.(varNonMAP{kk})
%         apu(kk) = true;
%     end
% end
% varNonMAP = varNonMAP(apu);
% for kk = 1 : numel(varNonMAP)
%     pz = im_vectors.(varNonMAP{kk});
% %     gg = gg + 1;
%     break
% end

% varMAP = recNames(2);
% varPrior = recNames(1);
% apu = false(numel(varMAP),1);
% for kk = 1 : numel(varMAP)
%     if options.(varMAP{kk})
%         apu(kk) = true;
%     end
% end
% varMAP = varMAP(apu);
% apu = false(numel(varPrior),1);
% for kk = 1 : numel(varPrior)
%     if options.(varPrior{kk})
%         apu(kk) = true;
%     end
% end
% varPrior = varPrior(apu);
% for ll = 1 : numel(varPrior)
%     for kk = 1 : numel(varMAP)
%         pz = im_vectors.([varPrior{ll} '_' varMAP{kk}]);
% %         gg = gg + 1;
%         break
%     end
% end
if nargout > 1
    varargout{1} = gg;
end