function im_vectors = initialize_im_vectors(im_vectors, iter, options)

if options.save_iter
    iter_n = iter;
else
    iter_n = 1;
end

fn = fieldnames(im_vectors);
fn2 = fn(cellfun('isempty',strfind(fn,'apu')));
fn = fn(~cellfun('isempty',strfind(fn,'apu')));
for kk = 1 : numel(fn)
    im_vectors.(fn{kk}) = options.x0;
end