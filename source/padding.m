function A = padding(A,sizeP,varargin)
%PADDING Pads the input array symmetrically or with zeros
%   This function pads the input array A symmetrically. The amount of to
%   be padded is given by the matrix sizeP. The elements of sizeP give the
%   amount of padding performed on each side. E.g. sizeP = [2 2] pads
%   both X- and Y-directions with mirrored versions of the input data on
%   each side. Both A and sizeP can be either 2D or 3D matrices. Padding
%   is always performed symmetrically.
%
%   Zero padding can be achieved by inputting 'zeros' after the size.
if nargin > 3 && (isempty(varargin{1}) || ~strcmp(varargin{1},'zeros'))
    A = [flipud(A(1:sizeP(2),:,:));A;flipud(A(end-sizeP(2) + 1:end,:,:))];
    A = [fliplr(A(:,1:sizeP(1),:)),A,fliplr(A(:,end-sizeP(1) + 1:end,:))];
    if length(sizeP) == 3 && sizeP(3) ~= 0
        A = cat(3, flip(A(:,:,1:sizeP(3)),3), A);
        A = cat(3, A, flip(A(:,:,end - sizeP(3) + 1: end),3));
    end
else
    [x, y , z] = size(A);
    A = [zeros(sizeP(1),y + sizeP(2)*2,z);zeros(x,sizeP(2),z),A,zeros(x,sizeP(2),z);zeros(sizeP(1),y+sizeP(2)*2,z)];
    if length(sizeP) == 3 && sizeP(3) ~= 0
        A = cat(3, zeros(size(A,1),size(A,2),sizeP(3)), A);
        A = cat(3, A, zeros(size(A,1),size(A,2),sizeP(3)));
    end
end
% A = A(:);
end