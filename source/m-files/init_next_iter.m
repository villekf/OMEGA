function im_vectors = init_next_iter(im_vectors, options, iter, varargin)
% Initialize the next iteration

% Output slot (column) of this iteration, same layout as in the C++ code (reconstructionAF.h):
% save_iter: column 1 holds x0 and iteration k is stored to column k + 1.
% Otherwise the k-th iteration in saveNIter (zero-based iteration numbers) is stored to column k and the final
% iteration always to the last column numel(saveNIter) + 1. Other iterations are not stored (iter_n = 0)
if options.save_iter
    iter_n = iter + 1;
elseif iter == options.Niter
    iter_n = numel(options.saveNIter) + 1;
else
    iter_n = find(options.saveNIter == iter - 1, 1);
    if isempty(iter_n)
        iter_n = 0;
    end
end
if nargin > 3 && ~isempty(varargin{1})
    tt = varargin{1};
else
    tt = 1;
end

if options.BSREM || options.ROSEM_MAP
    % The prior gradient is recomputed only every regEveryIter-th iteration and the last computed gradient is used in
    % between (the first and last iteration are always computed, as in iterStep.h)
    regEvery = 1;
    if isfield(options, 'regEveryIter') && ~isempty(options.regEveryIter)
        regEvery = double(options.regEveryIter);
    end
    if regEvery <= 1 || iter == 1 || iter == options.Niter || mod(iter - 1, regEvery) == 0 || ~isfield(im_vectors, 'dU')
        if iscell(im_vectors.recApu)
            dU = applyPrior(im_vectors.recApu{1}, options, 1, iter);
        else
            dU = applyPrior(im_vectors.recApu, options, 1, iter);
        end
        im_vectors.dU = dU;
    else
        dU = im_vectors.dU;
    end
    im_vectors = MAPiter(im_vectors, options.lambda(iter), options.beta, dU, options.epps);
end

% Positivity clamp of the 1-subset image before it is saved, as in the C++ code
if options.enforcePositivity && options.subsets == 1 && ~options.CPType
    if iscell(im_vectors.recApu)
        for ii = 1 : options.nMultiVolumes + 1
            im_vectors.recApu{ii}(im_vectors.recApu{ii} < options.epps) = options.epps;
        end
    else
        im_vectors.recApu(im_vectors.recApu < options.epps) = options.epps;
    end
end

% With LSQR and CGLS recApu holds the Krylov vector except after the last iteration, so the solution iterate
% (fLSQR/fCGLS) is saved instead, as in the C++ code (iterStep.h)
saveSrc = im_vectors.recApu;
if options.LSQR && isfield(im_vectors, 'fLSQR') && ~isempty(im_vectors.fLSQR)
    saveSrc = im_vectors.fLSQR;
elseif options.CGLS && isfield(im_vectors, 'fCGLS') && ~isempty(im_vectors.fCGLS)
    saveSrc = im_vectors.fCGLS;
end
if iter_n > 0
    if iscell(saveSrc)
        for ii = 1:numel(saveSrc)
            im_vectors.recImage{ii}(:, iter_n, tt) = saveSrc{ii};
        end
    else
        im_vectors.recImage(:, iter_n, tt) = saveSrc;
    end
end

if options.verbose > 0
    disp(['Iteration ' num2str(iter) ' finished'])
end
