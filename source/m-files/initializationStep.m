function [im_vectors, mData, options] = initializationStep(options, mData, im_vectors, m_size, iter, subIter, timestep)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here{ii}
if (options.param.FISTA || options.param.FISTAL1)
    if (iter == 1 && subIter == 1)
        im_vectors.uFISTA = im_vectors.recApu;
    else
        im_vectors.recApu = im_vectors.uFISTA;
    end
end
if (iter == 1)
    if (options.param.verbose >= 3)
        disp("Starting initialization step");
    end
    if (options.param.LSQR && subIter == 1)
        if (options.param.verbose >= 3)
            disp("Initializing LSQR");
        end
        im_vectors.fLSQR = im_vectors.recApu;
        options.param.betaLSQR = norm(mData);
        mData = mData / options.param.betaLSQR;
        if options.param.implementation == 1
            A = formMatrix(options);
            im_vectors.rhs = A * mData;
        else
            im_vectors.rhs = options' * mData;
        end
        temp = cell2mat(im_vectors.rhs);
        for ll = 1 : options.param.nMultiVolumes + 1
            im_vectors.recApu{timestep, ll} = im_vectors.rhs{timestep, ll} / options.param.alphaLSQR;
        end
        options.param.alphaLSQR = norm(temp);

        im_vectors.wLSQR = im_vectors.recApu;
        options.param.phiLSQR = options.param.betaLSQR;
        options.param.rhoLSQR = options.param.alphaLSQR;
        if (options.param.verbose >= 3)
            disp("LSQR initialization complete");
        end
    elseif (options.param.CGLS && subIter == 1)
        if (options.param.verbose >= 3)
            disp("Initializing CGLS");
        end
        im_vectors.rCGLS = mData;
        im_vectors.fCGLS = im_vectors.recApu;
        if options.param.implementation == 1
            A = formMatrix(options);
            im_vectors.rhs = A * mData;
        else
            im_vectors.rhs = options' * mData;
        end
        im_vectors.recApu = im_vectors.rhs;
        options.param.gammaCGLS = 0;
        for ll = 1 : options.param.nMultiVolumes + 1
            options.param.gammaCGLS = options.param.gammaCGLS + (im_vectors.rhs{timestep, ll}' * im_vectors.rhs{timestep, ll});
        end

        if (options.param.verbose >= 3)
            disp("CGLS initialization complete");
        end
    end
    if (options.param.CPType)
        if (options.param.verbose >= 3)
            disp("Initializing PDHG algorithm");
        end
        if ~isfield(im_vectors,'pCP')
            im_vectors.pCP = cell(options.param.partitions, options.param.subsets);
        end
        if (options.param.subsets > 1)
            % if ~isfield(im_vectors,'p0CP')
            %     im_vectors.p0CP = cell(options.param.subsets, 1);
            % end
            if ~isfield(im_vectors,'fpCP')
                im_vectors.fpCP = cell(options.param.partitions, options.param.subsets);
            end
        end
        if subIter == 1
            im_vectors.uCP = im_vectors.recApu;
        end
        im_vectors.pCP{timestep, subIter} = zeros(m_size, 1, options.param.cType);
        if (options.param.subsets > 1)
            if subIter == 1
                im_vectors.p0CP = zeros(m_size, 1, options.param.cType);
            end
            im_vectors.fpCP{timestep, subIter} = zeros(m_size, 1, options.param.cType);
        end
        if (options.param.verbose >= 3)
            disp("PDHG initialization complete");
        end
    end
    if subIter == 1 && (options.param.COSEM || options.param.ECOSEM || options.param.ACOSEM || options.param.OSL_COSEM > 0)
        im_vectors.C_co = zeros(options.param.N(1), options.param.subsets, options.param.cType);
    end
end