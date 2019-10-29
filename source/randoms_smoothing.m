function randoms = randoms_smoothing(randoms,options)
%RANDOMS_SMOOTHING Performs a moving mean smoothing on randoms or scatter
%data
%   Input randoms/scatter sinogram/raw list-mode data is smoothed by using
%   a 8x8 moving mean smoothing. The parameters Ndx and Ndy can be used to
%   adjust the size of the mean window. Smoothing is done symmetrically by
%   taking the mirror images of the measurement data on all sides. Outputs
%   the smoothed sinogram/raw list-mode data.
if options.verbose
    disp('Beginning randoms/scatter smoothing')
end
Ndx = 6;
Ndy = 6;
Ndz = 0;
N_pad = min(2, Ndx + Ndy + Ndz);
c1 = cell(N_pad,1);
c2 = cell(N_pad,1);
[c1{1:N_pad}]=ndgrid(1:(Ndx*2+1));
c2(1:N_pad)={Ndy+1};
if options.use_raw_data
    ind = tril(true(options.detectors,options.detectors),0);
    padd = zeros(options.detectors, options.detectors);
    padd(ind) = full(randoms);
    padd = padd' + padd;
    padd = padding(padd,[Ndx Ndy Ndz]);
    tmp = zeros(options.det_w_pseudo+Ndx*2, options.det_w_pseudo+Ndy*2, options.rings^2);
    hh = 1;
    for kk = 1 : options.rings
        for ll = kk : options.rings
            tmp(:,:,hh) = padd(1 + (ll - 1)*(options.det_w_pseudo): (options.det_w_pseudo)*ll + Ndx*2,1 + (kk - 1)*(options.det_w_pseudo): (options.det_w_pseudo)*kk + Ndx*2);
            hh = hh + 1;
        end
    end
    padd = zeros(options.detectors, options.detectors);
    s = [options.det_w_pseudo + Ndx*2 options.det_w_pseudo + Ndy*2 1 + Ndz*2];
%     tmp = padding(tmp,[Ndx Ndy Ndz]);
    tmp = squeeze(num2cell(tmp,[1 2]));
    N = options.det_w_pseudo * options.det_w_pseudo;
    offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
    tr_ind = sub2ind([options.det_w_pseudo + Ndx*2 options.det_w_pseudo + Ndy*2],mod((1:N)'-1,options.det_w_pseudo)+(Ndx + 1),mod(floor(((1:double(N))'-1)/double(options.det_w_pseudo)),double(options.det_w_pseudo))+(Ndy + 1));
    tr_offsets = uint32(bsxfun(@plus,tr_ind,offsets(:)'));
    hh = 1;
    tic
    for kk = 1 : options.rings
        for ll = kk : options.rings
            grad = mean(double(tmp{hh}(tr_offsets)),2);
            padd(1 + (ll - 1)*options.det_w_pseudo: options.det_w_pseudo*ll,1 + (kk - 1)*options.det_w_pseudo: options.det_w_pseudo*kk) = reshape(grad, options.det_w_pseudo, options.det_w_pseudo);
            hh = hh + 1;
        end
    end
    toc
    randoms = padd(ind);
else
    s = [size(randoms,1) + Ndx*2 size(randoms,2) + Ndy*2 size(randoms,3) + Ndz*2];
    N = options.Ndist * options.Nang;
    offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
    tr_ind = sub2ind([options.Ndist + Ndx*2 options.Nang + Ndy*2],mod((1:N)'-1,options.Ndist)+(Ndx + 1),mod(floor(((1:double(N))'-1)/double(options.Ndist)),double(options.Nang))+(Ndy + 1));
    tr_offsets = uint32(bsxfun(@plus,tr_ind,offsets(:)'));
    padd = padding(randoms,[Ndx Ndy Ndz]);
    padd = squeeze(num2cell(padd,[1 2]));
    for kk = 1 : options.TotSinos
        grad = mean(padd{kk}(tr_offsets),2);
        randoms(:,:,kk) = reshape(grad, options.Ndist, options.Nang);
    end
    
end
if options.verbose
    disp('Smoothing complete')
end