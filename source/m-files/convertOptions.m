function options = convertOptions(options)
%CONVERTOPTIONS Converts some fields in options to correspond to the new
%version(s) of OMEGA
%   Guarantees backwards compatibility

if isfield(options,'mlem')
    options.MLEM = options.mlem;
    options = rmfield(options,'mlem');
end
if isfield(options,'osem')
    options.OSEM = options.osem;
    options = rmfield(options,'osem');
end
if isfield(options,'mramla')
    options.MRAMLA = options.mramla;
    options = rmfield(options,'mramla');
end
if isfield(options,'ramla')
    options.RAMLA = options.ramla;
    options = rmfield(options,'ramla');
end
if isfield(options,'rosem')
    options.ROSEM = options.rosem;
    options = rmfield(options,'rosem');
end
if isfield(options,'rbi')
    options.RBI = options.rbi;
    options = rmfield(options,'rbi');
end
if isfield(options,'drama')
    options.DRAMA = options.drama;
    options = rmfield(options,'drama');
end
if isfield(options,'cosem')
    options.COSEM = options.cosem;
    options = rmfield(options,'cosem');
end
if isfield(options,'ecosem')
    options.ECOSEM = options.ecosem;
    options = rmfield(options,'ecosem');
end
if isfield(options,'acosem')
    options.ACOSEM = options.acosem;
    options = rmfield(options,'acosem');
end
if isfield(options,'RBI_OSL')
    options.OSL_RBI = options.RBI_OSL;
    options = rmfield(options,'RBI_OSL');
end
if isfield(options,'COSEM_OSL')
    options.OSL_COSEM = options.COSEM_OSL;
    options = rmfield(options,'COSEM_OSL');
end
if isfield(options,'lambda0_mbsrem')
    options.lambda0_MBSREM = options.lambda0_mbsrem;
    options = rmfield(options,'lambda0_mbsrem');
end
if isfield(options,'lambda0_rosem')
    options.lambda0_ROSEM = options.lambda0_rosem;
    options = rmfield(options,'lambda0_rosem');
end
if ~isfield(options,'lambda')
    options.lambda = 0;
end
if (options.RAMLA || options.BSREM) && isfield(options,'lambda0') && sum(options.lambda) == 0
    options.lambda = options.lambda0;
end
if (options.ROSEM || options.ROSEM_MAP) && isfield(options,'lambda0_ROSEM') && sum(options.lambda) == 0
    options.lambda = options.lambda0_ROSEM;
end
if (options.PKMA) && isfield(options,'lambda0_PKMA') && sum(options.lambda) == 0
    options.lambda = options.lambda0_PKMA;
end
if (options.MRAMLA || options.MBSREM) && isfield(options,'lambda0_MBSREM') && sum(options.lambda) == 0
    options.lambda = options.lambda0_MBSREM;
end
if ~isfield(options,'custom')
    options.custom = false;
end
varNeg = recNames(4);
for kk = 1 : numel(varNeg)
    if options.(varNeg{kk}) && options.enforcePositivity
        options.enforcePositivity = false;
        break;
    end
end
if options.projector_type < 4
    varProj1 = recNames(3);
    for kk = 1 : numel(varProj1)
        if options.(varProj1{kk})
            if options.projector_type == 1
                options.projector_type = 11;
            elseif options.projector_type == 2
                options.projector_type = 22;
            elseif options.projector_type == 3
                options.projector_type = 33;
            end
            break;
        end
    end
end
% if ~isfield(options,'layerOffset')
%      options.nLayers = 0;
% end
if ~isfield(options,'beta')
    options.beta = 0;
end
if ~isfield(options, 'beta_temporal')
    options.beta_temporal = 0;
end
origPrior = {'mrp';'quad';'huber';'L';'fmh';'weighted';'TV';'ad';'APLS';'TGV';'NLM';'custom'};
origMAP = {'osem';'bsrem';'mbsrem';'rosem';'rbi';'cosem'};
varPrior = recNames(1);
varMAP = recNames(2);
for kk = 1 : numel(origPrior)
    for ll = 1 : numel(origMAP)
        if isfield(options,['beta_' origPrior{kk} '_' origMAP{ll}])
            options.(['beta_' varPrior{kk} '_' varMAP{ll}]) = options.(['beta_' origPrior{kk} '_' origMAP{ll}]);
            options = rmfield(options,['beta_' origPrior{kk} '_' origMAP{ll}]);
        end
    end
end
for ll = 1 : numel(varMAP)
    if ~isfield(options,varMAP{ll})
        options.(varMAP{ll}) = false;
    end
end
for ll = 1 : numel(varPrior)
    if ~isfield(options,varPrior{ll})
        options.(varPrior{ll}) = false;
    end
end
for kk = 1 : numel(varPrior)
    for ll = 1 : numel(varMAP)
        if options.(varMAP{ll}) && options.(varPrior{kk}) && ~isfield(options,['beta_' varPrior{kk} '_' varMAP{ll}])
            options.(['beta_' varPrior{kk} '_' varMAP{ll}]) = options.beta;
        end
    end
end
end