function im_vectors = init_next_iter(im_vectors, options, iter, varargin)
% Initialize the next iteration

if options.save_iter
    iter_n = iter + 1;
else
    % saveNIter contains zero-based iteration numbers. The output axis is
    % compact: requested iterations occupy their ordinal slots, followed by
    % one reserved slot for the final reconstruction.
    ind = find(options.saveNIter == iter - 1, 1, 'first');
    if ~isempty(ind)
        iter_n = ind;
    elseif iter == options.Niter
        iter_n = numel(options.saveNIter) + 1;
    else
        iter_n = 0;
    end
end
if nargin > 3 && ~isempty(varargin{1})
    tt = varargin{1};
else
    tt = 1;
end

if options.BSREM || options.ROSEM_MAP
    if iscell(im_vectors.recApu)
        dU = applyPrior(im_vectors.recApu{1}, options, 1, iter);
    else
        dU = applyPrior(im_vectors.recApu, options, 1, iter);
    end
    im_vectors = MAPiter(im_vectors, options.lambda(iter), options.beta, dU, options.epps);
end

if iter_n > 0
    if iscell(im_vectors.recApu)
        im_vectors.recImage = im_vectors.recApu;
    else
        im_vectors.recImage(:, iter_n, tt) = im_vectors.recApu;
    end
end

if options.verbose > 0
    disp(['Iteration ' num2str(iter) ' finished'])
end
