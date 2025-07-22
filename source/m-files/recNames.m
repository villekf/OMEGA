function varList = recNames(varargin)
%RECNAMES Obtain the OMEGA reconstruction variable names 
%   0 returns all variable names in uppercase
%   1 returns all prior names
%   2 returns all MAP names
%   3 returns all algorithms that are not built into projector type 1/2/3
%   4 returns all algorithms that should not ever use positivity constraint
%   5 returns all algorithms requiring linearized data
%   6 returns all non-MAP names
%   7 returns all algorithms that support image-based preconditioning
%   8 returns all algorithms that support measurement-based preconditioning
%   9 returns all algorithms not supported by implementation 3
%   10 returns all full prior names, for display purposes only
%   11 returns all full algorithm names, for display purposes only
%   otherwise the variable names are returned as in the options-struct
%   Add your own algorithm/prior name to the lists below

% Add Any MAP algorithm here (algorithms that can use gradient-based
% priors)
varMAP = {'OSL_OSEM';'BSREM';'MBSREM';'ROSEM_MAP';'OSL_RBI';'OSL_COSEM';'PKMA';'SPS';'PDHG';'PDHGKL';'PDHGL1';'CV';'PDDY';'SART';'ASD_POCS';'SAGA';'BB'};
% Add any non-MAP/prior-based algorithm to this list (or algorithms with
% fixed prior)
OS = {'OSEM';'MRAMLA';'RAMLA';'ROSEM';'RBI';'DRAMA';'COSEM';'ECOSEM';'ACOSEM';'LSQR';'CGLS';'FISTA'; 'FISTAL1';'FDK'};
% Algorithms that are not built into the OpenCL kernel of projector types
% 1/2/3. This will automatically select e.g. projector_type = 11 instead of
% projector_type = 1.
varProj1 = {'LSQR';'CGLS';'FISTA';'FISTAL1';'SPS';'PDHG';'PDHGKL';'PDHGL1';'PDDY';'SART';'ASD_POCS';'SAGA';'BB';'ASD_POCS'};
% Algorithms that should not have positivity constraint
varNeg = {'LSQR';'CGLS';'FDK';'SART'};
% Algorithms that support image-based preconditioners
varPreCondIm = {'MRAMLA'; 'MBSREM'; 'FISTA';'FISTAL1'; 'PKMA';'SPS';'PDHG';'PDHGKL';'PDHGL1';'PDDY';'SAGA'};
% Algorithms that support measurement-based preconditioners
varPreCondMeas = {'MRAMLA'; 'MBSREM';'FISTA';'PKMA';'SPS'; 'PDHG';'PDHGKL';'PDHGL1';'FISTAL1';'PDDY';'SAGA'};
% Algorithms that require linearized data
varLin = {'LSQR';'CGLS'; 'FISTA'; 'FISTAL1';'PDHG';'PDHGL1';'PDDY';'FDK';'SART';'ASD_POCS';'BB'};
% Add your new prior to this list (before 'custom')
varPrior = {'MRP';'quad';'Huber';'L';'FMH';'weighted_mean';'TV';'hyperbolic';'AD';'APLS';'TGV';'NLM';'RDP';'GGMRF';'ProxTV';'ProxRDP';'ProxNLM';'custom'};
% Add your prior to this list also (this is used for display purposes only
% and can be unabbreviated and contain any characters):
varPriorFull = {'Median Root';'Quadratic';'Huber';'L-filter';'FIR Median Hybrid';'Weighted mean';'Total Variation';'Hyperbolic function';'Anisotropic Diffusion';...
    'Asymmetric Parallel Level Sets';'Total Generalized Variation';'Non-Local Means';'Relative Difference';'General Gaussian Markov Random Field';...
    'Proximal Total variation';'Proximal Relative Difference';'Proximal Non-Local Means';'Custom'};
% Add your algorithm to this list also (this is used for display purposes only
% and can be unabbreviated and contain any characters):
varAlgoFull = {'Ordered Subsets Expectation Maximization (OSEM)'; 'Modified Row-action Maximum-Likelihood Algorithm (MRAMLA)'; 'Row-action Maximum-Likelihood Algorithm (RAMLA)'; ...
    'Relaxed Ordered Subsets Expectation Maximization (ROSEM)'; 'Rescaled Block-iterative algorithm (RBI)'; 'Dynamic Row-action Maximum-Likelihood Algorithm (DRAMA)';...
    'Complete-data Ordered Subsets Expectation Maximization (COSEM)'; 'Enhanced Complete-data Ordered Subsets Expectation Maximization (ECOSEM)';...
    'Accelerated Complete-data Ordered Subsets Expectation Maximization (ACOSEM)'; 'LSQR'; 'Conjugate-gradient least-squares (CGLS)'; ...
    'Fast Iterative Shrinkage Thresholding Algorithm (FISTA)'; 'Fast Iterative Shrinkage Thresholding Algorithm with L1 regularization (FISTAL1)';...
    'FDK/FBP algorithm';'Simultaneous Algebraic Reconstruction Technique';'One-step-late OSEM (OSL_OSEM)'; ...
    'Block-sequential regularized expectation maximization (BSREM)'; 'Modified block-sequential regularized expectation maximization (MBSREM)';...
    'ROSEM MAP (ROSEM_MAP)'; 'One-step-late RBI (OSL_RBI)'; 'One-step-late COSEM (OSL_COSEM)'; 'Preconditioned Krasnoselski-Mann algorithm (PKMA)';...
    'Separable paraboloidal surrogates (SPS)'; 'Primal-dual hybrid gradient (PDHG)'; 'Primal-dual hybrid gradient with Kullback-Leibler divergence (PDHGKL)';...
    'Primal-dual hybrid gradient with L1 error (PDHGL1)';'Condat-Vu algorithm';'Primal-dual Davis-Yin algorithm';'ASD-POCS';'SAGA';'Barzilia-Borwein'};
if ~isempty(varargin) && varargin{1} == 0
    varList = upper([OS(:);varMAP(:)]);
elseif ~isempty(varargin) && varargin{1} == 1
    varList = varPrior;
elseif ~isempty(varargin) && varargin{1} == 2
    varList = varMAP;
elseif ~isempty(varargin) && varargin{1} == 3
    varList = varProj1;
elseif ~isempty(varargin) && varargin{1} == 4
    varList = varNeg;
elseif ~isempty(varargin) && varargin{1} == 5
    varList = varLin;
elseif ~isempty(varargin) && varargin{1} == 6
    varList = OS;
elseif ~isempty(varargin) && varargin{1} == 7
    varList = varPreCondIm;
elseif ~isempty(varargin) && varargin{1} == 8
    varList = varPreCondMeas;
elseif ~isempty(varargin) && varargin{1} == 9
    varList = [OS(:);varMAP(:)];
    varList = varList(~strcmp(varList,'OSEM'));
elseif ~isempty(varargin) && varargin{1} == 10
    varList = varPriorFull;
elseif ~isempty(varargin) && varargin{1} == 11
    varList = varAlgoFull;
else
    varList = [OS(:);varMAP(:)];
end
end