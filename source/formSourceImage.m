function [C, SC, RA] = formSourceImage(options, S, jj, trues_index, scatter_index, randoms_index, C, SC, RA, ll, index)
%FORMSOURCEIMAGE Forms the original radioactive source image
% An ideal image can be formed with this
% Takes the columns containing the source locations for both singles
pixel_width_x = options.FOVa_x/options.Nx;
pixel_width_y = options.FOVa_y/options.Ny;
z_width = options.axial_fov/options.Nz;
x = single((-pixel_width_x*(options.Nx-1)/2:pixel_width_x:pixel_width_x*(options.Nx-1)/2)');
y = single((-pixel_width_y*(options.Ny-1)/2:pixel_width_y:pixel_width_y*(options.Ny-1)/2)');
z = single((-z_width*(options.Nz-1)/2:z_width:z_width*(options.Nz-1)/2)');
max_x = max(x) + pixel_width_x/2;
max_y = max(y) + pixel_width_y/2;
max_z = max(z) + z_width/2;

% Form the source image by using only coincidences that originate from the
% very same location (source coordinates are the same)
CC = ones(length(S),1,'int16');
t1 = all(S(:,1:3)==S(:,4:6),2);
t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
t = (t1+t2+t3==3);
CC = CC(t);
% Find the minimum distance from each pixel center
% The smallest distance means it is in the same
% pixel
[~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
[~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
[~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
%                 FOV = circshift(FOV, round(size(FOV,1)/6));
C{1,jj} = C{1,jj} + FOV;

% Form the source image by using only the first single
CC = ones(length(S),1,'int16');
t1 = true(size(S,1),1);
t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
t = (t1+t2==2);
CC = CC(t);
[~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
[~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
[~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
%                 FOV = circshift(FOV, round(size(FOV,1)/6));
C{2,jj} = C{2,jj} + FOV;

% Form the source image by using only the second single
CC = ones(length(S),1,'int16');
t1 = true(size(S,4),1);
t2 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
t = (t1+t2==2);
CC = CC(t);
[~,i] = min(abs(bsxfun(@minus,x.',S(t,4))),[],2);
[~,j] = min(abs(bsxfun(@minus,y.',S(t,5))),[],2);
[~,k] = min(abs(bsxfun(@minus,z.',S(t,6))),[],2);
FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
%                 FOV = circshift(FOV, round(size(FOV,1)/6));
C{3,jj} = C{3,jj} + FOV;

% Form the source image by using both singles (singles mode)
C{4,jj} = C{2,jj} + C{3,jj};

% Form the source image by using both singles and then dividing the counts
% by two
CC = ones(length(S),1,'int16');
t1 = true(size(S,1),1);
t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
t = (t1+t2+t3==3);
CC = CC(t);
[~,i] = min(abs(bsxfun(@minus,(x).',S(t,1))),[],2);
[~,j] = min(abs(bsxfun(@minus,(y).',S(t,2))),[],2);
[~,k] = min(abs(bsxfun(@minus,(z).',S(t,3))),[],2);
[~,i2] = min(abs(bsxfun(@minus,(x).',S(t,4))),[],2);
[~,j2] = min(abs(bsxfun(@minus,(y).',S(t,5))),[],2);
[~,k2] = min(abs(bsxfun(@minus,(z).',S(t,6))),[],2);
i = [i;i2];
j = [j;j2];
k = [k;k2];
CC = repmat(CC,uint32(2),uint32(1));
FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
%                 FOV = circshift(FOV, round(size(FOV,1)/6));
FOV = FOV/2;
C{5,jj} = C{5,jj} + FOV;

% Form the source image by using the average coordinates from both singles
CC = ones(length(S),1,'int16');
t1 = true(size(S,1),1);
t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
t = (t1+t2+t3==3);
CC = CC(t);
C1 = mean([S(t,1) S(t,4)],2);
C2 = mean([S(t,2) S(t,5)],2);
C3 = mean([S(t,3) S(t,6)],2);
[~,i] = min(abs(bsxfun(@minus,(x).',C1)),[],2);
[~,j] = min(abs(bsxfun(@minus,(y).',C2)),[],2);
[~,k] = min(abs(bsxfun(@minus,(z).',C3)),[],2);
FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
%                 FOV = circshift(FOV, round(size(FOV,1)/6));
C{6,jj} = C{6,jj} + FOV;

% Form the source image by using the true coincidences (requires
% obtain_trues = true)
if options.obtain_trues
    CC = ones(length(S),1,'int16');
    CC = CC(trues_index(ll:index));
    SS = S(trues_index(ll:index),:);
    t = all([abs(SS(:,1)) <= max_x, abs(SS(:,2)) <= max_y, abs(SS(:,3)) <= max_z],2);
    CC = CC(t);
    % Find the minimum distance from each pixel center
    % The smallest distance means it is in the same
    % pixel
    [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
    [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
    [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
    if isempty(i)
        FOV = zeros([options.Nx options.Ny options.Nz],'uint16');
    else
        FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
    end
    % FOV = circshift(FOV, round(size(FOV,1)/6));
    C{7,jj} = C{7,jj} + FOV;
    clear SS
end
% Form the source image for the scatter data (singles mode)
if options.store_scatter
    CC = ones(length(S),1,'int16');
    CC = CC(scatter_index(ll:index));
    SS = S(scatter_index(ll:index),:);
    t= all([abs(SS(:,1))<=max_x abs(SS(:,2))<=max_y abs(SS(:,3))<=max_z],2);
    CC = CC(t);
    [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
    [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
    [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
    if isempty(i)
        FOV = zeros([options.Nx options.Ny options.Nz],'uint16');
    else
        FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
    end
    % FOV = circshift(FOV, round(size(FOV,1)/6));
    SC = SC + FOV;
    CC = ones(length(S),1,'int16');
    CC = CC(scatter_index(ll:index));
    SS = S(scatter_index(ll:index),:);
    t= all([abs(SS(:,4))<=max_x abs(SS(:,5))<=max_y abs(SS(:,6))<=max_z],2);
    CC = CC(t);
    [~,i] = min(abs(bsxfun(@minus,x.',SS(t,4))),[],2);
    [~,j] = min(abs(bsxfun(@minus,y.',SS(t,5))),[],2);
    [~,k] = min(abs(bsxfun(@minus,z.',SS(t,6))),[],2);
    if isempty(i)
        FOV = zeros([options.Nx options.Ny options.Nz],'uint16');
    else
        FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
    end
    % FOV = circshift(FOV, round(size(FOV,1)/6));
    SC = SC + FOV;
    clear SS
end
% Form the source image for the randoms data (singles mode)
if options.store_randoms
    CC = ones(length(S),1,'int16');
    CC = CC(randoms_index(ll:index));
    SS = S(randoms_index(ll:index),:);
    t = all([abs(SS(:,1)) <= max_x abs(SS(:,2)) <= max_y abs(SS(:,3))<=max_z],2);
    CC = CC(t);
    [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
    [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
    [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
    if isempty(i)
        FOV = zeros([options.Nx options.Ny options.Nz],'uint16');
    else
        FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
    end
    % FOV = circshift(FOV, round(size(FOV,1)/6));
    RA = RA + FOV;
    CC = ones(length(S),1,'int16');
    CC = CC(randoms_index(ll:index));
    t= all([abs(SS(:,4))<=max_x abs(SS(:,5))<=max_y abs(SS(:,6))<=max_z],2);
    CC = CC(t);
    [~,i] = min(abs(bsxfun(@minus,x.',SS(t,4))),[],2);
    [~,j] = min(abs(bsxfun(@minus,y.',SS(t,5))),[],2);
    [~,k] = min(abs(bsxfun(@minus,z.',SS(t,6))),[],2);
    if isempty(i)
        FOV = zeros([options.Nx options.Ny options.Nz],'uint16');
    else
        FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[options.Nx options.Ny options.Nz])),1));
    end
    % FOV = circshift(FOV, round(size(FOV,1)/6));
    RA = RA + FOV;
    clear SS
end

