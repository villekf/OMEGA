function im_vectors = MAPiter(im_vectors, lam, beta, dU, epps)

if iscell(im_vectors.recApu)
    im_vectors.recApu{1} = im_vectors.recApu{1} - beta .* lam .* im_vectors.recApu{1} .* dU;
    im_vectors.recApu{1}(im_vectors.recApu{1} < epps) = epps;
else
    im_vectors.recApu = im_vectors.recApu - beta .* lam .* im_vectors.recApu .* dU;
    im_vectors.recApu(im_vectors.recApu < epps) = epps;
end