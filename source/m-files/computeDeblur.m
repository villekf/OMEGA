function im_vectors = computeDeblur(im_vectors, options, iter, gaussK, Nx, Ny, Nz, varargin)
%COMPUTEDEBLUR Computes the PSF deblurring phase for all selected OS
%algorithms

if options.save_iter
    iter_n = iter;
else
    iter_n = 0;
end
if nargin > 7 && ~isempty(varargin{1})
    tt = varargin{1};
else
    tt = 1;
end

im_vectors.recImage(:, iter_n + 1, tt) = deblur(im_vectors.recImage(:, iter_n + 1, tt), options, gaussK, Nx, Ny, Nz);
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