function im_vectors = init_next_iter(im_vectors, options, iter, varargin)
% Initialize the next iteration

if options.save_iter
    iter_n = iter + 1;
else
    ind = ismember(options.saveNIter, iter - 1);
    if any(ind)
        iter_n = options.saveNIter(ind) + 1;
    else
        iter_n = 1;
    end
end
if nargin > 3 && ~isempty(varargin{1})
    tt = varargin{1};
else
    tt = 1;
end

if options.BSREM || options.ROSEM_MAP
    if iscell(im_vectors.recApu)
        dU = applyPrior(im_vectors.recApu{1}, options, iter, options.beta);
    else
        dU = applyPrior(im_vectors.recApu, options, iter, options.beta);
    end
    im_vectors = MAPiter(im_vectors, options.lambda(iter), options.beta, dU, options.epps);
end

if iscell(im_vectors.recApu)
    im_vectors.recImage(:, iter_n, tt) = im_vectors.recApu{1};
else
    im_vectors.recImage(:, iter_n, tt) = im_vectors.recApu;
end

if options.verbose > 0
    disp(['Iteration ' num2str(iter) ' finished'])
end