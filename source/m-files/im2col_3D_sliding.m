function out = im2col_3D_sliding(A,blocksize)
% 3D implementation of (sliding) im2col
% Original code taken from:
% https://stackoverflow.com/questions/36761794/rearrange-sliding-blocks-into-columns-for-a-3d-array-im2col-in-3d-matlab

%// Store blocksizes
nrows = blocksize(1);
ncols = blocksize(2);
nslices = blocksize(3);

%// Get sizes for later usages
[m,n,r] = size(A);

%// Start indices for each block
start_ind = reshape(bsxfun(@plus,[1:m-nrows+1]',[0:n-ncols]*m),[],1); %//'

%// Row indices
lin_row = permute(bsxfun(@plus,start_ind,[0:nrows-1])',[1 3 2]);  %//'

%// 2D linear indices
lidx_2D = reshape(bsxfun(@plus,lin_row,[0:ncols-1]*m),nrows*ncols,[]);

%// 3D linear indices
lidx_3D = bsxfun(@plus,lidx_2D,m*n*permute((0:(r-nslices)),[1 3 2]));

lidx_3D = repmat(lidx_3D,nslices,1,1);
apu = (0:m*n:m*n*nslices-1)';
if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
    apu = repeat_elem(apu,nrows*ncols,0);
else
    apu = repelem(apu,nrows*ncols);
end
lidx_3D = bsxfun(@plus, lidx_3D, apu);

%// Get linear indices based on row and col indices and get desired output
out = A(lidx_3D);
