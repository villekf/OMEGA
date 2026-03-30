function im_vectors = MAPiter(im_vectors, lam, beta, dU, epps, timestep)
    im_vectors.recApu{timestep, 1} = im_vectors.recApu{timestep, 1} - beta .* lam .* im_vectors.recApu{timestep, 1} .* dU;
    im_vectors.recApu{timestep, 1}(im_vectors.recApu{timestep, 1} < epps) = epps;
end