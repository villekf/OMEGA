function dU = applyPrior(recApu, options, beta, osa_iter)

dU = [];
if (options.MRP)
    if (options.verbose >= 3)
        disp("Computing MRP gradient");
    end
    dU = beta * MRP(recApu, options.medx, options.medy, options.medz, options.Nx(1), options.Ny(1), options.Nz(1), options.epps, ...
        options.tr_offsets, options.med_no_norm);
elseif (options.quad)
    if (options.verbose >= 3)
        disp("Computing quadratic prior gradient");
    end
    dU = beta * Quadratic_prior(recApu, options.weights_quad, options.Nx(1), options.Ny(1), options.Nz(1), options.Ndx, options.Ndy, options.Ndz);
elseif (options.Huber)
    if (options.verbose >= 3)
        disp("Computing Huber prior gradient");
    end
    dU = beta * Huber_prior(recApu, options.weights_huber, options.Nx(1), options.Ny(1), options.Nz(1), options.Ndx, options.Ndy, options.Ndz, options.huber_delta);
elseif (options.L)
    if (options.verbose >= 3)
        disp("Computing L-filter gradient");
    end
    dU = beta * L_filter(recApu, options.tr_offsets, options.a_L, options.Nx(1), options.Ny(1), options.Nz(1), options.Ndx, options.Ndy, options.Ndz, ...
        options.epps, options.med_no_norm);
elseif (options.FMH)
    if (options.verbose >= 3)
        disp("Computing FMH prior gradient");
    end
    dU = beta * FMH(recApu, options.tr_offsets, options.fmh_weights, options.Nx(1), options.Ny(1), options.Nz(1), options.N(1), options.Ndx, options.Ndy, options.Ndz, options.epps, ...
        options.med_no_norm);
elseif (options.weighted_mean)
    if (options.verbose >= 3)
        disp("Computing weighted mean prior gradient");
    end
    dU = beta * Weighted_mean(recApu, options.weighted_weights, options.Nx(1), options.Ny(1), options.Nz(1), options.Ndx, options.Ndy, options.Ndz, ...
        options.mean_type, options.epps, options.med_no_norm);
elseif (options.TV)
    if (options.verbose >= 3)
        disp("Computing TV prior gradient");
    end
    dU = beta * TVpriorFinal(recApu, options.TVdata, options.Nx(1), options.Ny(1), options.Nz(1), options.TV_use_anatomical, options, options.TVtype, ...
        options.tr_offsets);
elseif (options.AD)
    if (options.verbose >= 3)
        disp("Computing AD prior gradient");
    end
    if osa_iter > 1
        dU = beta * AD(recApu, options.FluxType, options.Nx(1), options.Ny(1), options.Nz(1), options);
    else
        dU = zeros(options.Nx(1) * options.Ny(1) * options.Nz(1), 1);
    end
elseif (options.APLS)
    if (options.verbose >= 3)
        disp("Computing APLS prior gradient");
    end
    dU = beta * TVpriorFinal(recApu, [], options.Nx(1), options.Ny(1), options.Nz(1), true, options, 5);
elseif (options.TGV)
    if (options.verbose >= 3)
        disp("Computing TGV prior");
    end
    warning('Use of TGV with implementation 1, 4 and 5 is NOT recommended!')
    dU = beta * TGV(recApu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx(1), options.Ny(1), options.Nz(1));
    elseif (options.ProxTV)
        error('Proximal TV is currently unsupported with implementations 1, 4 and 5')
        % if (options.verbose >= 3)
        %     disp("Computing proximal TV prior");
        % end
        % status = proxTV(vec.im_os[0], options, vec, proj, w_vec, dU, w_vec.betaReg);
elseif (options.NLM)
    if (options.verbose >= 3)
        disp("Computing NLM prior gradient");
    end
    dU = beta * NLM(recApu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
        options.sigma, options.epps, options.Nx(1), options.Ny(1), options.Nz(1), options);
elseif (options.RDP)
    if (options.verbose >= 3)
        disp("Computing RDP prior gradient");
    end
    dU = beta * RDP(recApu, options.weights_RDP, options.RDP_gamma, options.Nx(1), options.Ny(1), options.Nz(1), options.Ndx, options.Ndy, options.Ndz, options.tr_offsets);
elseif (options.GGMRF)
        error('GGMRF is currently unsupported with implementations 1, 4 and 5. As a workaround, use NLM with NLGGMRF and zero neighborhood and patch size.')
    %     if (options.verbose >= 3)
    %         disp("Computing GGMRF prior gradient");
    %     end
    %     status = GGMRF(vec.im_os[0], options, w_vec.GGMRF_p, w_vec.GGMRF_q, w_vec.GGMRF_c, w_vec.GGMRF_pqc, proj, dU, beta);
    % elseif (options.ProxRDP)
    %     if (options.verbose >= 3)
    %         disp("Computing proximal RDP prior");
    %     end
    %     status = proxRDP(vec.im_os[0], options, vec, proj, dU, beta, w_vec.RDP_gamma);
    % elseif (options.ProxNLM)
    %     if (options.verbose >= 3)
    %         disp("Computing proximal NLM prior");
    %     end
    %     status = proxNLM(vec.im_os[0], options, vec, proj, dU, w_vec, beta);
end
if (options.verbose >= 3 && (options.MRP || options.quad || options.Huber || options.L || options.FMH || options.TV ...
        || options.weighted_mean || options.AD || options.APLS || options.TGV || options.NLM || options.RDP ...
        || options.ProxTV || options.GGMRF))
    disp("Prior computed");
end