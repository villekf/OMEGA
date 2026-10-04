function options = computeWeights(options, varargin)
%COMPUTEWEIGTHS Computes the weights required by QP, L-filter, FMH, TV type
%3, RDP, GGMRF and weighted mean
%   Computes the weights based on the neighborhood size and the pixel size
if isempty(options.weights)
    distX = options.FOVa_x(1)/double(options.Nx(1));
    distY = options.FOVa_y(1)/double(options.Ny(1));
    distZ = (double(options.axial_fov(1))/double(options.Nz(1)));
    % Offset vectors, each running from +N down to -N (matches the element
    % order the original loop-based implementation produced).
    xr = (options.Ndx : -1 : -options.Ndx)' * distX;
    yr = (options.Ndy : -1 : -options.Ndy)' * distY;
    zr = (options.Ndz : -1 : -options.Ndz)' * distZ;
    if nargin > 1 && varargin{1} == true
        % GGMRF-style ordering: z varies fastest, then y, then x slowest.
        [Zg, Yg, Xg] = ndgrid(zr, yr, xr);
        if options.Ndx == 0 || options.Nx(1) == 1
            dist = sqrt(Zg.^2 + Yg.^2);
        else
            dist = sqrt(Zg.^2 + Yg.^2 + Xg.^2);
        end
    else
        % Default ordering: x varies fastest, then y, then z slowest.
        [Xg, Yg, Zg] = ndgrid(xr, yr, zr);
        if options.Ndz == 0 || options.Nz(1) == 1
            dist = sqrt(Xg.^2 + Yg.^2);
        else
            dist = sqrt(Xg.^2 + Yg.^2 + Zg.^2);
        end
    end
    options.weights = 1./dist(:);
    % summa = sum(options.weights(~isinf(options.weights)));
    % options.weights = options.weights / summa;
end
end
