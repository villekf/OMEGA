function options = computeWeights(options, varargin)
%COMPUTEWEIGTHS Computes the weights required by QP, L-filter, FMH, TV type
%3, RDP, GGMRF and weighted mean
%   Computes the weights based on the neighborhood size and the pixel size
if isempty(options.weights)
    distX = options.FOVa_x(1)/double(options.Nx(1));
    distY = options.FOVa_y(1)/double(options.Ny(1));
    distZ = (double(options.axial_fov(1))/double(options.Nz(1)));
    options.weights = zeros(((options.Ndx*2+1) * (options.Ndy*2+1) * (options.Ndz*2+1)),1);
    edist = zeros((options.Ndx*2+1),1);
    cc = zeros((options.Ndy*2+1)*(options.Ndx*2+1),1);
    lt = 0;
    if nargin > 1 && varargin{1} == true
        for jj = options.Ndx : -1 : -options.Ndx
            lt = lt + 1;
            ll = 0;
            for kk = options.Ndy : -1 : -options.Ndy
                ll = ll + 1;
                if options.Ndx == 0 || options.Nx(1) == 1
                    if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                        apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repeat_elem(kk,options.Ndy*2+1) * distY)];
                    elseif exist('OCTAVE_VERSION','builtin') == 5
                        apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repelem(kk,options.Ndy*2+1) * distY)];
                    else
                        apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repelem(kk,options.Ndy*2+1) * distY)'];
                    end
                else
                    if options.Ndx ~= options.Ndz
                        if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                            apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repeat_elem(kk,options.Ndy*2+1) * distY), [zeros(options.Ndz-options.Ndx,1),(repeat_elem(jj,options.Ndx*2+1) * distX),...
                                zeros(options.Ndx-options.Ndx,1)]];
                        elseif exist('OCTAVE_VERSION','builtin') == 5
                            apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repelem(kk,options.Ndy*2+1) * distY), [zeros(options.Ndz-options.Ndx,1),(repelem(jj,options.Ndx*2+1) * distX),...
                                zeros(options.Ndx-options.Ndx,1)]];
                        else
                            apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repelem(kk,options.Ndy*2+1) * distY)', [zeros(options.Ndz-options.Ndx,1),(repelem(jj,options.Ndx*2+1) * distX),...
                                zeros(options.Ndx-options.Ndx,1)]'];
                        end
                    else
                        if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                            apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repeat_elem(kk,options.Ndy*2+1) * distY), (repeat_elem(jj,options.Ndx*2+1) * distX)];
                        elseif exist('OCTAVE_VERSION','builtin') == 5
                            apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repelem(kk,options.Ndy*2+1) * distY), (repelem(jj,options.Ndx*2+1) * distX)];
                        else
                            apu = [((options.Ndz:-1:-options.Ndz) * distZ)', (repelem(kk,options.Ndy*2+1) * distY)', (repelem(jj,options.Ndx*2+1) * distX)'];
                        end
                    end
                end
                for ii = 1 : length(apu)
                    edist(ii) = sqrt(apu(ii,:)*apu(ii,:)');
                end
                cc((options.Ndy*2+1)*(ll-1)+1:(options.Ndy*2+1)*ll) = edist;
            end
            options.weights((options.Ndz*2+1) * (options.Ndy*2+1) * (lt - 1) + 1: (options.Ndz*2+1) * (options.Ndy*2+1) * lt) = cc;
        end
    else
        for jj = options.Ndz : -1 : -options.Ndz
            lt = lt + 1;
            ll = 0;
            for kk = options.Ndy : -1 : -options.Ndy
                ll = ll + 1;
                if options.Ndz == 0 || options.Nz(1) == 1
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
    end
    options.weights = 1./options.weights;
    % summa = sum(options.weights(~isinf(options.weights)));
    % options.weights = options.weights / summa;
end
end

