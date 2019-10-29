function A = padding(A,sizeP)
%PADDING Pads the input array symmetrically
%   This function pads the input array A symmetrically. The amount of to
%   be padded is given by the matrix sizeP. The elements of sizeP give the
%   amount of padding performed on each side. E.g. sizeP = [2 2] pads
%   both X- and Y-directions with mirrored versions of the input data on
%   each side. Both A and sizeP can be either 2D or 3D matrices. Padding
%   is always performed symmetrically.
% [x, y , z] = size(A);
A = [flipud(A(1:sizeP(2),:,:));A;flipud(A(end-sizeP(2) + 1:end,:,:))];
A = [fliplr(A(:,1:sizeP(1),:)),A,fliplr(A(:,end-sizeP(1) + 1:end,:))];
% A = [zeros(sizeP(1),y + sizeP(2)*2,z);zeros(x,sizeP(2),z),A,zeros(x,sizeP(2),z);zeros(sizeP(1),y+sizeP(2)*2,z)];
if length(sizeP) == 3 && sizeP(3) ~= 0
    A = cat(3, flip(A(:,:,1:sizeP(3)),3), A);
    A = cat(3, A, flip(A(:,:,end - sizeP(3) + 1: end),3));
end
% A = A(:);
end