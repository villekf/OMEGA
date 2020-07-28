function options = fmhWeights(options)
%FMHWEIGHTS Compute weights for the FMH prior
%   options.weights are needed as the input data. They can be formed with
%   computeWeights.

distX = options.FOVa_x/double(options.Nx);
distY = options.FOVa_y/double(options.Ny);
distZ = (double(options.axial_fov)/double(options.Nz));
if isempty(options.fmh_weights)
    kerroin = options.fmh_center_weight^(1/4)*distX;
    % 2D case
    if options.Nz == 1 || options.Ndz == 0
        options.fmh_weights = zeros(options.Ndx*2+1, 4);
        lll = 0;
        for jjj = 1 : 4
            lll = lll + 1;
            apu = zeros(options.Ndx*2+1,1);
            hhh = 0;
            % There are 4 different combinations where the
            % means are computed in FMH
            if jjj == 1 || jjj == 3
                for iii = options.Ndx : -1 : -options.Ndx
                    hhh = hhh + 1;
                    if iii == 0
                        apu(hhh) = options.fmh_center_weight;
                    else
                        apu(hhh) = kerroin/sqrt((distX*iii)^2+(distY*iii)^2);
                    end
                end
            elseif jjj == 2
                for iii = options.Ndx : -1 : -options.Ndx
                    hhh = hhh + 1;
                    if iii == 0
                        apu(hhh) = options.fmh_center_weight;
                    else
                        apu(hhh) = kerroin/abs(distX*iii);
                    end
                end
            elseif jjj == 4
                for iii = options.Ndx : -1 : -options.Ndx
                    hhh = hhh + 1;
                    if iii == 0
                        apu(hhh) = options.fmh_center_weight;
                    else
                        apu(hhh) = kerroin/abs(distY*iii);
                    end
                end
            end
            options.fmh_weights(:, jjj) = apu;
        end
    else
        % 3D case
        options.fmh_weights = zeros(max([options.Ndx*2+1,options.Ndz*2+1]), 13);
        lll = 0;
        for kkk = 1 : -1 : 0
            for jjj = 1 : 9
                lll = lll + 1;
                % 9 cases in 3D + the 2D cases
                if kkk == 1
                    apu = zeros(options.Ndz*2+1,1);
                    hhh = 0;
                    if jjj == 1 || jjj == 3 || jjj == 7 || jjj == 9
                        for iii = options.Ndz : -1 : -options.Ndz
                            hhh = hhh + 1;
                            if iii == 0
                                apu(hhh) = options.fmh_center_weight;
                            else
                                apu(hhh) = kerroin/sqrt(sqrt((distZ*iii)^2+(distX*iii)^2)^2+(distY*iii)^2);
                            end
                        end
                    elseif jjj == 2 || jjj == 8
                        for iii = options.Ndz : -1 : -options.Ndz
                            hhh = hhh + 1;
                            if iii == 0
                                apu(hhh) = options.fmh_center_weight;
                            else
                                apu(hhh) = kerroin/sqrt((distZ*iii)^2+(distX*iii)^2);
                            end
                        end
                    elseif jjj == 4 || jjj == 6
                        for iii = options.Ndz : -1 : -options.Ndz
                            hhh = hhh + 1;
                            if iii == 0
                                apu(hhh) = options.fmh_center_weight;
                            else
                                apu(hhh) = kerroin/sqrt((distZ*iii)^2+(distY*iii)^2);
                            end
                        end
                    elseif jjj == 5
                        for iii = options.Ndz : -1 : -options.Ndz
                            hhh = hhh + 1;
                            if iii == 0
                                apu(hhh) = options.fmh_center_weight;
                            else
                                apu(hhh) = kerroin/abs(distZ*iii);
                            end
                        end
                    end
                    options.fmh_weights(:, lll) = apu;
                else
                    % Same as in 2D case
                    apu = zeros(options.Ndx*2+1,1);
                    hhh = 0;
                    if jjj == 1 || jjj == 3
                        for iii = options.Ndx : -1 : -options.Ndx
                            hhh = hhh + 1;
                            if iii == 0
                                apu(hhh) = options.fmh_center_weight;
                            else
                                apu(hhh) = kerroin/sqrt((distX*iii)^2+(distY*iii)^2);
                            end
                        end
                    elseif jjj == 2
                        for iii = options.Ndx : -1 : -options.Ndx
                            hhh = hhh + 1;
                            if iii == 0
                                apu(hhh) = options.fmh_center_weight;
                            else
                                apu(hhh) = kerroin/abs(distX*iii);
                            end
                        end
                    elseif jjj == 4
                        for iii = options.Ndx : -1 : -options.Ndx
                            hhh = hhh + 1;
                            if iii == 0
                                apu(hhh) = options.fmh_center_weight;
                            else
                                apu(hhh) = kerroin/abs(distY*iii);
                            end
                        end
                    else
                        break
                    end
                    options.fmh_weights(:, lll) = apu;
                end
            end
        end
    end
    options.fmh_weights = options.fmh_weights./sum(options.fmh_weights,1);
end
if options.implementation == 2
    options.fmh_weights = single(options.fmh_weights);
end