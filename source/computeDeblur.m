function im_vectors = computeDeblur(im_vectors, options, iter, gaussK, Nx, Ny, Nz)
%COMPUTEDEBLUR Computes the PSF deblurring phase for all selected OS
%algorithms

if options.save_iter
    iter_n = iter;
else
    iter_n = 0;
end
if options.implementation == 4
    fn = fieldnames(im_vectors);
    fn = fn(cellfun('isempty',strfind(fn,'apu')));
    
    for kk = 1 : numel(fn)
        im_vectors.(fn{kk})(:, iter_n + 1) = deblur(im_vectors.(fn{kk})(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
    end
elseif options.implementation == 1
    
    fn = fieldnames(im_vectors);
    fn = fn(cellfun('isempty',strfind(fn,'apu')));
    for kk = 1 : numel(fn)
        im_vectors.(fn{kk})(:, iter_n + 1) = deblur(im_vectors.(fn{kk})(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
    end
end