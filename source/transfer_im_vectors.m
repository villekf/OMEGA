function im_vectors = transfer_im_vectors(im_vectors, pz, options, iter)

% if options.save_iter
%     iter_n = iter + 1;
% else
    iter_n = 1;
% end
oo = find(~cellfun('isempty',pz));

fn = fieldnames(im_vectors);
fn2 = fn(~cellfun('isempty',strfind(fn,'apu')));
fn = fn(cellfun('isempty',strfind(fn,'apu')));
for kk = 1 : numel(fn)
    im_vectors.(fn2{kk}) = reshape(pz{oo(kk)}(:,:,:,iter_n),options.Nx*options.Ny*options.Nz,1);
    im_vectors.(fn{kk}) = reshape(pz{oo(kk)}(:,:,:,iter_n),options.Nx*options.Ny*options.Nz,1);
end