function [im_vectors, options] = computeOSEstimates(im_vectors, options, iter, osa_iter, uu, corrim)
for ii = 1 : options.param.nMultiVolumes + 1

    if (options.param.precondTypeImage(6) && options.paramilteringIterations > 0 && iter == options.paramilteringIterations)
        if (options.param.verbose >= 3)
            disp("Image-based filter iterations complete.");
        end
        if (options.param.CPType || options.param.FISTA || options.param.FISTAL1 || options.param.TGV || options.param.ProxTV)
            if (options.param.sigmaCP(ii) == 1)
                options.param.tauCP(ii) = options.param.tauCP2(ii);
            elseif (options.param.sigmaCP(ii) == options.param.tauCP(ii))
                options.param.tauCP(ii) = options.param.tauCP2(ii);
                options.param.sigmaCP(ii) = options.param.tauCP2(ii);
            else
                options.param.sigmaCP(ii) = options.param.tauCP2(ii);
            end
            if (options.param.CPType)
                options.param.sigma2CP = options.param.sigmaCP;
            end
        end
        if (options.param.MRAMLA || options.param.MBSREM || options.param.SPS || options.param.RAMLA || options.param.BSREM || options.param.ROSEM || options.param.ROSEM_MAP || options.param.PKMA)
            options.param.lambda = options.param.lambdaFiltered;
        end
        options.param.precondTypeImage(6) = false;
    end
    if (options.param.PDHG || options.param.PDHGKL || options.param.PDHGL1 || options.param.CV || options.param.PDDY)
        [im_vectors] = PDHG1(options.param, im_vectors, osa_iter + options.param.subsets * (iter - 1), ii);
    end
    if (options.param.PDDY && ii == 1)
        if iscell(im_vectors.recApu)
            PDDYApu = im_vectors.recApu{1};
            im_vectors.recApu{1} = im_vectors.recApu{1} - options.param.tauCP(1) .* im_vectors.uCP{1};
        else
            PDDYApu = im_vectors.recApu;
            im_vectors.recApu = im_vectors.recApu - options.param.tauCP(1) .* im_vectors.uCP;
        end
        if (options.param.verbose == 3)
            disp("Computing PDDY step\n");
        end
    end
end

if iscell(im_vectors.recApu)
    dU = applyPrior(im_vectors.recApu{1}, options.param, options.param.beta, osa_iter);
else
    dU = applyPrior(im_vectors.recApu, options.param, options.param.beta, osa_iter);
end
if (options.param.PDDY)
    if iscell(im_vectors.recApu)
        im_vectors.recApu{1} = PDDYApu;
    else
        im_vectors.recApu = PDDYApu;
    end
end

for ii = 1 : options.param.nMultiVolumes + 1
    if iscell(im_vectors.Sens)
        if size(im_vectors.Sens{ii},2) > 1
            Sens = im_vectors.Sens{ii}(:, osa_iter);
        else
            Sens = im_vectors.Sens{ii};
        end
    else
        if size(im_vectors.Sens,2) > 1
            Sens = im_vectors.Sens(:, osa_iter);
        else
            Sens = im_vectors.Sens;
        end
    end
    if iscell(im_vectors.recApu)
        im = im_vectors.recApu{ii};
        rhs = im_vectors.rhs{ii};
        if isfield(options.param, 'D')
            D = options.param.D{ii};
        end
    else
        im = im_vectors.recApu;
        rhs = im_vectors.rhs;
        if isfield(options.param, 'D')
            D = options.param.D;
        end
    end

    % Compute the (matrix free) algorithms
    % Ordered Subsets Expectation Maximization (OSEM)
    if (options.param.OSEM || options.param.ECOSEM)
        if (options.param.verbose >= 3)
            disp("Computing OSEM/ECOSEM");
        end
        if (options.param.ECOSEM)
            OSEMApu = EM(im, Sens, rhs + options.param.epps);
        else
            if (options.param.CT)
                im = EM(im, Sens, options.param.flat * rhs);
            else
                im = EM(im, Sens, rhs + options.param.epps);
            end
        end
    end

    % Modfied Row-action Maximum Likelihood (MRAMLA)
    if (options.param.MRAMLA)
        if (options.param.verbose >= 3)
            disp("Computing MRAMLA");
        end
        im = MBSREM(im, rhs, options.param, options.param.U, options.param.lambda, iter, osa_iter, ii);
    end

    % Row-action Maximum Likelihood (RAMLA)
    if (options.param.RAMLA)
        if (options.param.verbose >= 3)
            disp("Computing RAMLA");
        end
        im = BSREM(im, rhs, options.param.lambda, iter);
    end

    % Relaxed OSEM (ROSEM)
    if (options.param.ROSEM)
        if (options.param.verbose >= 3)
            disp("Computing ROSEM");
        end
        if (options.param.CT)
            im = ROSEM(im, Sens, options.param.flat * rhs, options.param.lambda, iter);
        else
            im = ROSEM(im, Sens, rhs, options.param.lambda, iter);
        end
    end

    % Rescaled Block Iterative EM (RBI)
    if (options.param.RBI)
        if (options.param.verbose >= 3)
            disp("Computing RBI");
        end
        im = RBI(im, Sens, rhs, D);
    end

    % Dynamic RAMLA
    if (options.param.DRAMA)
        if (options.param.verbose >= 3)
            disp("Computing DRAMA");
        end
        im = DRAMA(im, Sens, rhs, options.param.lambda, iter, osa_iter, options.param.subsets);
    end

    % Complete data OSEM
    if (options.param.COSEM || options.param.ECOSEM)
        if (options.param.verbose >= 3)
            disp("Computing COSEM/ECOSEM");
        end
        if (options.param.ECOSEM)
            [COSEMApu, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, D, options.param.h, 2, osa_iter);
        else
            [im, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, D, options.param.h, 2, osa_iter);
        end
    end

    % Enhanced COSEM
    if (options.param.ECOSEM)
        if (options.param.verbose >= 3)
            disp("Computing ECOSEM");
        end
        im = ECOSEM(im, D, OSEMApu, COSEMApu, options.param.epps);
    end

    % Accelerated COSEM
    if (options.param.ACOSEM)
        if (options.param.verbose >= 3)
            disp("Computing ACOSEM");
        end
        [im, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, D, options.param.h, 1, osa_iter);
        ACOSEM_rhs = computeACOSEMWeight(options, im, uu, osa_iter, corrim, ii);
        im = im .* ACOSEM_rhs;
    end
    if (options.param.OSL_OSEM)
        if (options.param.verbose >= 3)
            disp("Computing OSL-OSEM");
        end
        if (options.param.CT)
            im = EM(im, Sens + dU, options.param.flat * rhs);
        else
            im = EM(im, Sens + dU, rhs);
        end
    elseif (options.param.BSREM)
        if (options.param.verbose >= 3)
            disp("Computing BSREM");
        end
        im = BSREM(im, rhs, options.param.lambda, iter);
    elseif (options.param.MBSREM)
        if (options.param.verbose >= 3)
            disp("Computing MBSREM");
        end
        if isempty(dU)
            im = MBSREM(im, rhs, options.param, options.param.U, options.param.lambda, iter, osa_iter, ii);
        else
            im = MBSREM(im, rhs + dU, options.param, options.param.U, options.param.lambda, iter, osa_iter, ii);
        end
    elseif (options.param.ROSEM_MAP)
        if (options.param.verbose >= 3)
            disp("Computing ROSEM-MAP");
        end
        if (options.param.CT)
            im = ROSEM(im, Sens, options.param.flat * rhs, options.param.lambda, iter);
        else
            im = ROSEM(im, Sens, rhs, options.param.lambda, iter);
        end
    elseif (options.param.OSL_RBI)
        if (options.param.verbose >= 3)
            disp("Computing OSL-RBI");
        end
        im = RBI(im, Sens, rhs, D, options.param.beta, dU);
    elseif (options.param.OSL_COSEM > 0)
        if (options.param.verbose >= 3)
            disp(["Computing OSL-COSEM ", num2str(options.param.OSL_COSEM)]);
        end
        [im, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, D + dU, options.param.h, options.param.OSL_COSEM, osa_iter);
        if (options.param.OSL_COSEM == 1)
            ACOSEM_rhs = computeACOSEMWeight(options, im, uu, osa_iter, corrim, ii);
            im = im .* ACOSEM_rhs;
        end
    elseif (options.param.PKMA)
        if (options.param.verbose >= 3)
            disp("Computing PKMA");
        end
        if isempty(dU)
            im = PKMA(im, rhs, options.param, iter, osa_iter, ii);
        else
            im = PKMA(im, rhs + dU, options.param, iter, osa_iter, ii);
        end
    elseif (options.param.SPS)
        if (options.param.verbose >= 3)
            disp("Computing SPS");
        end
        im = SPS(im, rhs, options.param.U, options.param.lambda, iter, osa_iter, ii, options.param);
    elseif (options.param.LSQR)
        if (options.param.verbose >= 3)
            disp("Computing LSQR");
        end
        [im_vectors,options.param] = LSQRim(options.param, iter, im_vectors, ii);
    elseif (options.param.CGLS)
        if (options.param.verbose >= 3)
            disp("Computing CGLS");
        end
        [im_vectors,options.param] = CGLS(options.param, iter, im_vectors, ii);
    elseif (options.param.PDHG || options.param.PDHGKL || options.param.PDHGL1 || options.param.CV || options.param.PDDY)
        if (options.param.verbose >= 3)
            disp("Computing PDHG/PDHGKL/PDHGL1");
        end
        if isempty(dU)
            [im, options.param, im_vectors] = PDHG2(im, rhs, im_vectors, options.param, iter, osa_iter, ii);
        else
            [im, options.param, im_vectors] = PDHG2(im, rhs + dU, im_vectors, options.param, iter, osa_iter, ii);
        end
    elseif (options.param.FISTA)
        if (options.param.verbose >= 3)
            disp("Computing FISTA");
        end
        [im_vectors, im,options.param] = FISTA(im, rhs, options.param, im_vectors, osa_iter + options.param.subsets * (iter - 1), ii);
    elseif (options.param.FISTAL1)
        if (options.param.verbose >= 3)
            disp("Computing FISTAL1");
        end
        [im_vectors, im,options.param] = FISTAL1(im, rhs, options.param, im_vectors, options.param.beta, osa_iter + options.param.subsets * (iter - 1), ii);
    end
end
if ~options.param.LSQR && ~options.param.CGLS
    if iscell(im_vectors.recApu)
        im_vectors.recApu{ii} = im;
    else
        im_vectors.recApu = im;
    end
end