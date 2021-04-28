function options = convertOptions(options)
%CONVERTOPTIONS Converts some fields in options to correspond to the new
%version of OMEGA
%   Guarantees backwards compatibility

if ~isfield(options,'PKMA')
    options.PKMA = false;
end
if ~isfield(options,'RDP')
    options.RDP = false;
end

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
origPrior = {'mrp';'quad';'huber';'L';'fmh';'weighted';'TV';'ad';'APLS';'TGV';'NLM';'custom'};
origMAP = {'mlem';'osem';'bsrem';'mbsrem';'rosem';'rbi';'cosem'};
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
end