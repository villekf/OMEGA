function im_vectors = computeDeblur(im_vectors, options, iter, gaussK, Nx, Ny, Nz, varargin)
%COMPUTEDEBLUR Computes the PSF deblurring phase for all selected OS
%algorithms

% Same output slot as in init_next_iter (the final iteration is stored to the last column)
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
if nargin > 7 && ~isempty(varargin{1})
    tt = varargin{1};
else
    tt = 1;
end

if iter_n > 0
    im_vectors.recImage(:, iter_n, tt) = deblur(im_vectors.recImage(:, iter_n, tt), options, gaussK, Nx, Ny, Nz);
end
% if options.implementation == 4
    % fn = fieldnames(im_vectors);
    % fn = fn(cellfun('isempty',strfind(fn,'Apu')));
    % 
    % for kk = 1 : numel(fn)
    %     im_vectors.(fn{kk})(:, iter_n + 1, tt) = deblur(im_vectors.(fn{kk})(:, iter_n + 1, tt), options, gaussK, Nx, Ny, Nz);
    % end
% elseif options.implementation == 1
%     
%     fn = fieldnames(im_vectors);
%     fn = fn(cellfun('isempty',strfind(fn,'Apu')));
%     for kk = 1 : numel(fn)
%         im_vectors.(fn{kk})(:, iter_n + 1, tt) = deblur(im_vectors.(fn{kk})(:, iter_n + 1, tt), options, gaussK, Nx, Ny, Nz);
%     end
% end