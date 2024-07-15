function options = double_to_single(options)
%DOUBLE_TO_SINGLE Convert all necessary values in options to single
%precision
if isfield(options,'x0')
    options.x0 = single(options.x0(:));
else
    options.x0 = [];
end
% options.Ndx = uint32(options.Ndx);
% options.Ndy = uint32(options.Ndy);
% options.Ndz = uint32(options.Ndz);
% options.MLEM = logical(options.MLEM);
% options.OSEM = logical(options.OSEM);
% options.MRAMLA = logical(options.MRAMLA);
% options.RAMLA = logical(options.RAMLA);
% options.ROSEM = logical(options.ROSEM);
% options.RBI = logical(options.RBI);
% options.DRAMA = logical(options.DRAMA);
% options.COSEM = logical(options.COSEM);
% options.ECOSEM = logical(options.ECOSEM);
% options.ACOSEM = logical(options.ACOSEM);
% options.OSL_MLEM = logical(options.OSL_MLEM);
% options.OSL_OSEM = logical(options.OSL_OSEM);
% options.MBSREM = logical(options.MBSREM);
% options.BSREM = logical(options.BSREM);
% options.ROSEM_MAP = logical(options.ROSEM_MAP);
% options.OSL_RBI = logical(options.OSL_RBI);
% options.MRP = logical(options.MRP);
% options.quad = logical(options.quad);
% options.L = logical(options.L);
% options.FMH = logical(options.FMH);
% options.weighted_mean = logical(options.weighted_mean);
% options.TV = logical(options.TV);
% options.AD = logical(options.AD);
% options.APLS = logical(options.APLS);
% options.TGV = logical(options.TGV);

% varPrior = recNames(1);
% varMAP = recNames(2);

% apu = false(numel(varPrior),1);
% apu2 = false(numel(varMAP),1);
% if numel(varMAP) > 0
%     for kk = 1 : numel(varPrior)
%         for ll = 1 : numel(varMAP)
%             if options.(varPrior{kk}) && options.(varMAP{ll})
%                 apu(kk) = true;
%                 apu2(ll) = true;
%             end
%         end
%     end
% end
% varPrior = varPrior(apu);
% varMAP = varMAP(apu2);


% for kk = 1 : numel(varPrior)
%     for ll = 1 : numel(varMAP)
%         options.(['beta_' [varPrior{kk} '_' varMAP{ll}]]) = single(options.(['beta_' [varPrior{kk} '_' varMAP{ll}]]));
%     end
% end
end

