function [im_vectors,options] = LSQRim(options, iter, im_vectors, ii, timestep)
    im_vectors.recApu{timestep, ii} = im_vectors.rhs{ii} - options.betaLSQR * im_vectors.recApu{timestep, ii};
    if (ii == options.nMultiVolumes + 1)
        temp = cell2mat(im_vectors.recApu);
        options.alphaLSQR = norm(temp);
        for ll = 1 : options.nMultiVolumes + 1
            im_vectors.recApu{timestep, ll} = im_vectors.recApu{timestep, ll} / options.alphaLSQR;
        end

        rho_ = sqrt(options.rhoLSQR * options.rhoLSQR + options.betaLSQR * options.betaLSQR);
        c = options.rhoLSQR / rho_;
        s = options.betaLSQR / rho_;
        options.thetaLSQR = s * options.alphaLSQR;
        options.rhoLSQR = -c * options.alphaLSQR;
        phi_ = c * options.phiLSQR;
        options.phiLSQR = s * options.phiLSQR;
        for ll = 1 : options.nMultiVolumes + 1
            im_vectors.fLSQR{timestep, ll} = (phi_ / rho_) * im_vectors.wLSQR{timestep, ll} + im_vectors.fLSQR{timestep, ll};
            im_vectors.wLSQR{timestep, ll} = im_vectors.recApu{timestep, ll} - (options.thetaLSQR / rho_) * im_vectors.wLSQR{timestep, ll};
            if (iter == options.Niter)
                im_vectors.recApu{timestep, ll} = im_vectors.fLSQR{timestep, ll};
            end
        end
    end
end