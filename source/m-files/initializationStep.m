function [im_vectors, mData, options] = initializationStep(options, mData, im_vectors, m_size, iter, subIter, varargin)
%INITIALIZATIONSTEP Algorithm-specific initialization (initializationStep in functions.hpp)
%   The optional seventh input is the multiplicative correction vector
%   (corrVector) used in the forward projection of the initial image
%   (LSQR/CGLS initial residual b - A*x0)
if (options.param.FISTA || options.param.FISTAL1)
    if (iter == 1 && subIter == 1)
        % uFISTA holds the extrapolated point y_k, xFISTA the previous FISTA iterate x_k
        im_vectors.uFISTA = im_vectors.recApu;
        im_vectors.xFISTA = im_vectors.recApu;
    elseif (options.param.subsets == 1 || subIter == 1)
        % With subsets the extrapolated point is the starting point of the first subset of each iteration only, the
        % other subsets continue from the previous subiteration
        im_vectors.recApu = im_vectors.uFISTA;
    end
end
if (iter == 1)
    if (options.param.verbose >= 3)
        disp("Starting initialization step");
    end
    if ((options.param.LSQR || options.param.CGLS) && subIter == 1)
        % Krylov methods solve A*dx = b - A*x0 with x = x0 + dx, so with a nonzero initial image the initialization uses
        % the initial residual b - A*x0 instead of the measurements
        if iscell(im_vectors.recApu)
            x0NonZero = any(cellfun(@(x) any(x(:) ~= 0), im_vectors.recApu));
        else
            x0NonZero = any(im_vectors.recApu(:) ~= 0);
        end
        if options.param.implementation == 1
            A = formMatrix(options);
        end
        if x0NonZero
            if (options.param.verbose >= 3)
                disp("Computing the initial residual b - A*x0");
            end
            if options.param.implementation == 1
                if iscell(im_vectors.recApu)
                    x0 = im_vectors.recApu{1};
                else
                    x0 = im_vectors.recApu;
                end
                if options.param.use_psf
                    x0 = computeConvolution(x0, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                end
                fp0 = A' * x0;
            else
                if ~isempty(varargin) && ~isempty(varargin{1})
                    fp0 = forwardProject(options, im_vectors.recApu, subIter, varargin{1});
                else
                    fp0 = forwardProject(options, im_vectors.recApu, subIter);
                end
            end
            residual0 = mData(:) - fp0(:);
        else
            residual0 = mData;
        end
    end
    if (options.param.LSQR && subIter == 1)
        if (options.param.verbose >= 3)
            disp("Initializing LSQR");
        end
        im_vectors.fLSQR = im_vectors.recApu;
        mData = residual0;
        options.param.betaLSQR = norm(mData);
        mData = mData / options.param.betaLSQR;
        if options.param.implementation == 1
            im_vectors.rhs = A * mData;
        else
            im_vectors.rhs = options' * mData;
        end
        if iscell(im_vectors.rhs)
            temp = cell2mat(im_vectors.rhs);
            options.param.alphaLSQR = norm(temp);
            for ll = 1 : options.param.nMultiVolumes + 1
                im_vectors.recApu{ll} = im_vectors.rhs{ll} / options.param.alphaLSQR;
            end
        else
            options.param.alphaLSQR = norm(im_vectors.rhs);
            im_vectors.recApu = im_vectors.rhs / options.param.alphaLSQR;
        end
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
        % The residual of the initial image is the CGLS residual and the backprojected vector, mData is not modified
        im_vectors.rCGLS = residual0;
        im_vectors.fCGLS = im_vectors.recApu;
        if options.param.implementation == 1
            im_vectors.rhs = A * residual0;
        else
            im_vectors.rhs = options' * residual0;
        end
        im_vectors.recApu = im_vectors.rhs;
        options.param.gammaCGLS = 0;
        if iscell(im_vectors.rhs)
            for ll = 1 : options.param.nMultiVolumes + 1
                options.param.gammaCGLS = options.param.gammaCGLS + (im_vectors.rhs{ll}' * im_vectors.rhs{ll});
            end
        else
            options.param.gammaCGLS = im_vectors.rhs' * im_vectors.rhs;
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
            im_vectors.pCP = cell(options.param.subsets, 1);
        end
        if (options.param.subsets > 1)
            % if ~isfield(im_vectors,'p0CP')
            %     im_vectors.p0CP = cell(options.param.subsets, 1);
            % end
            if ~isfield(im_vectors,'fpCP')
                im_vectors.fpCP = cell(options.param.subsets, 1);
            end
        end
        if subIter == 1
            % With subsets uCP is the A^T p accumulator (p starts at zero); without subsets it is the primal iterate
            if options.param.subsets > 1
                if iscell(im_vectors.recApu)
                    im_vectors.uCP = cellfun(@(x) zeros(size(x), 'like', x), im_vectors.recApu, 'UniformOutput', false);
                else
                    im_vectors.uCP = zeros(size(im_vectors.recApu), 'like', im_vectors.recApu);
                end
            else
                im_vectors.uCP = im_vectors.recApu;
            end
        end
        im_vectors.pCP{subIter} = zeros(m_size, 1, options.param.cType);
        if (options.param.subsets > 1)
            if subIter == 1
                im_vectors.p0CP = zeros(m_size, 1, options.param.cType);
            end
            im_vectors.fpCP{subIter} = zeros(m_size, 1, options.param.cType);
        end
        if (options.param.verbose >= 3)
            disp("PDHG initialization complete");
        end
    end
    if subIter == 1 && (options.param.COSEM || options.param.ECOSEM || options.param.ACOSEM || options.param.OSL_COSEM > 0)
        im_vectors.C_co = zeros(options.param.N(1), options.param.subsets, options.param.cType);
    end
end