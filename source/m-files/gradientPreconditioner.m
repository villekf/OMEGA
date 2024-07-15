function [options] = gradientPreconditioner(options, input, ii, VAL)
f = zeros(options.Nx(ii), options.Ny(ii), options.Nz(ii));
g = zeros(options.Nx(ii), options.Ny(ii), options.Nz(ii));
h = zeros(options.Nx(ii), options.Ny(ii), options.Nz(ii));
[f, g, h] = computeGradient(input, options, f, g, h, options.derivType);
f = max(VAL, sqrt(f .* f + g .* g + h .* h) ./ mean(input(:)));
f = mean(f(:)) ./ f;
options.gradF{ii} = min(options.gradV2, max(options.gradV1, f));