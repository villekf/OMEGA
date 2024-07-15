function [im_vectors,options] = LSQRim(options, iter, im_vectors, ii)
% if (iter == 1)
%     if iscell(im_vectors.recApu)
%         im_vectors.wLSQR{ii} = im_vectors.recApu{ii};
%     else
%         im_vectors.wLSQR = im_vectors.recApu;
%     end
% end
if iscell(im_vectors.recApu)
    im_vectors.recApu{ii} = im_vectors.rhs{ii} - options.betaLSQR * im_vectors.recApu{ii};
else
    im_vectors.recApu = im_vectors.rhs - options.betaLSQR * im_vectors.recApu;
end
if (ii == options.nMultiVolumes + 1)
    if iscell(im_vectors.recApu)
        temp = cell2mat(im_vectors.recApu);
        options.alphaLSQR = norm(temp);
        for ll = 1 : options.nMultiVolumes + 1
            im_vectors.recApu{ll} = im_vectors.recApu{ll} / options.alphaLSQR;
        end
    else
        options.alphaLSQR = norm(im_vectors.recApu);
        im_vectors.recApu = im_vectors.recApu / options.alphaLSQR;
    end
    rho_ = sqrt(options.rhoLSQR * options.rhoLSQR + options.betaLSQR * options.betaLSQR);
    c = options.rhoLSQR / rho_;
    s = options.betaLSQR / rho_;
    options.thetaLSQR = s * options.alphaLSQR;
    options.rhoLSQR = -c * options.alphaLSQR;
    phi_ = c * options.phiLSQR;
    options.phiLSQR = s * options.phiLSQR;
    for ll = 1 : options.nMultiVolumes + 1
        if iscell(im_vectors.fLSQR)
            im_vectors.fLSQR{ll} = (phi_ / rho_) * im_vectors.wLSQR{ll} + im_vectors.fLSQR{ll};
            im_vectors.wLSQR{ll} = im_vectors.recApu{ll} - (options.thetaLSQR / rho_) * im_vectors.wLSQR{ll};
            if (iter == options.Niter)
                im_vectors.recApu{ll} = im_vectors.fLSQR{ll};
            end
        else
            im_vectors.fLSQR = (phi_ / rho_) * im_vectors.wLSQR + im_vectors.fLSQR;
            im_vectors.wLSQR = im_vectors.recApu - (options.thetaLSQR / rho_) * im_vectors.wLSQR;
            if (iter == options.Niter)
                im_vectors.recApu = im_vectors.fLSQR;
            end
        end
    end
end