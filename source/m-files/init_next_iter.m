function im_vectors = init_next_iter(im_vectors, options, iter, timestep)
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


if options.BSREM || options.ROSEM_MAP
    dU = applyPrior(im_vectors.recApu{timestep, 1}, options, iter, options.beta);
    im_vectors = MAPiter(im_vectors, options.lambda(iter), options.beta, dU, options.epps, timestep);
end
currentMultiVolume = 1; % TODO
im_vectors.recImage{currentMultiVolume}(:, iter_n, timestep) = im_vectors.recApu{timestep, 1};

if options.verbose > 0
    disp(['Iteration ' num2str(iter) ' finished'])
end