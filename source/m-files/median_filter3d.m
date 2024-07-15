function img = median_filter3d(img)

[Nx, Ny, Nz] = size(img);
Ndx = 1;
Ndy = 1;
Ndz = 1;
N = Nx*Ny*Nz;

s = [Nx + Ndx*2 Ny + Ndy*2 Nz + Ndz*2];
N_pad = min(3, Ndx + Ndy + Ndz);
[c1{1:N_pad}]=ndgrid(1:(Ndx*2+1));
c2(1:N_pad)={Ndy+1};
if Ndz > Ndx && Ndz > 1
    c1{1} = cat(3, c1{1}, zeros(size(c1{1},1), size(c1{1},2), Ndz));
    c1{2} = cat(3, c1{2}, zeros(size(c1{2},1), size(c1{2},2), Ndz));
    c1{3} = cat(3, c1{3}, zeros(size(c1{3},1), size(c1{3},2), Ndz));
    for kk = Ndz - 1 : - 1 : 0
        c1{1}(:,:,end-kk) = c1{1}(:,:,end - kk - 1);
        c1{2}(:,:,end-kk) = c1{2}(:,:,end - kk - 1);
        c1{3}(:,:,end-kk) = c1{3}(:,:,end - kk - 1) + 1;
    end
    c2(end) = {Ndz+1};
elseif Ndz < Ndx && Ndz > 1
    c1{1}(:,:,end-2*(Ndx-Ndz) + 1) = [];
    c1{2}(:,:,end-2*(Ndx-Ndz) + 1) = [];
    c1{3}(:,:,end-2*(Ndx-Ndz) + 1) = [];
    c2(end) = {Ndz+1};
end
offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
if Nz == 1
    tr_ind = sub2ind([Nx+Ndx*2 Ny+Ndy*2],mod((1:N)'-1,Nx)+(Ndx + 1),mod(floor(((1:double(N))'-1)/double(Nx)),double(Ny))+(Ndy + 1));
else
    tr_ind = sub2ind([Nx+Ndx*2 Ny+Ndy*2 Nz+Ndz*2],mod((1:N)'-1,Nx)+(Ndx + 1),mod(floor(((1:double(N))'-1)/double(Nx)),double(Ny))+(Ndy + 1),floor(((1:double(N))'-1)/double(Nx*Ny))+(Ndz+1));
end
tr_offsets = uint32(bsxfun(@plus,tr_ind,offsets(:)'));
padd = padding(img,[max(0,Ndx) max(0,Ndy) max(0,Ndz)]);
padd = padd(tr_offsets);
img = median(padd,2);
img = reshape(img, Nx, Ny, Nz);