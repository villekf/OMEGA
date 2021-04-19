function x0 = updateInitialValue(im_vectors, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if options.OSEM
    x0 = im_vectors.OSEM_apu;
end
