function options = computeWeights(options)
%COMPUTEWEIGTHS Computes the weights required by QP, L-filter, FMH, TV type
%3 and weighted mean
%   Computes the weights based on the neighborhood size and the pixel size
distX = options.FOVa_x/double(options.Nx);
distY = options.FOVa_y/double(options.Ny);
distZ = (double(options.axial_fov)/double(options.Nz));
if isempty(options.weights)
    options.weights = zeros(((options.Ndx*2+1) * (options.Ndy*2+1) * (options.Ndz*2+1)),1);
    edist = zeros((options.Ndx*2+1),1);
    cc = zeros((options.Ndy*2+1)*(options.Ndx*2+1),1);
    lt = 0;
    for jj = options.Ndz : -1 : -options.Ndz
        lt = lt + 1;
        ll = 0;
        for kk = options.Ndy : -1 : -options.Ndy
            ll = ll + 1;
            if options.Ndz == 0 || options.Nz == 1
                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repeat_elem(kk,options.Ndy*2+1) * distY)];
                elseif exist('OCTAVE_VERSION','builtin') == 5
                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY)];
                else
                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY)'];
                end
            else
                if options.Ndz ~= options.Ndx
                    if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                        apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repeat_elem(kk,options.Ndy*2+1) * distY), [zeros(options.Ndx-options.Ndz,1),(repeat_elem(jj,options.Ndz*2+1) * distZ),...
                            zeros(options.Ndx-options.Ndz,1)]];
                    elseif exist('OCTAVE_VERSION','builtin') == 5
                        apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY), [zeros(options.Ndx-options.Ndz,1),(repelem(jj,options.Ndz*2+1) * distZ),...
                            zeros(options.Ndx-options.Ndz,1)]];
                    else
                        apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY)', [zeros(options.Ndx-options.Ndz,1),(repelem(jj,options.Ndz*2+1) * distZ),...
                            zeros(options.Ndx-options.Ndz,1)]'];
                    end
                else
                    if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                        apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repeat_elem(kk,options.Ndy*2+1) * distY), (repeat_elem(jj,options.Ndz*2+1) * distZ)];
                    elseif exist('OCTAVE_VERSION','builtin') == 5
                        apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY), (repelem(jj,options.Ndz*2+1) * distZ)];
                    else
                        apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY)', (repelem(jj,options.Ndz*2+1) * distZ)'];
                    end
                end
            end
            for ii = 1 : length(apu)
                edist(ii) = sqrt(apu(ii,:)*apu(ii,:)');
            end
            cc((options.Ndy*2+1)*(ll-1)+1:(options.Ndy*2+1)*ll) = edist;
        end
        options.weights((options.Ndx*2+1) * (options.Ndy*2+1) * (lt - 1) + 1: (options.Ndx*2+1) * (options.Ndy*2+1) * lt) = cc;
    end
    options.weights = 1./options.weights;
end
end

