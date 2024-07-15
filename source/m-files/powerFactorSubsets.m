function indices = powerFactorSubsets(nProj)
%POWERFACTORSUBSETS Computes subset indices based on the power factor of
%the input value
%   Detailed explanation goes here
% nProj = 500;
n = factor(nProj);
nn = flipud(cumprod(flipud(n(2:end)')));
N = numel(n);

p = zeros(N,nProj);

for ll = 1 : N
    p1 = 0:n(ll) - 1;
    if ll == 1
        p1 = repmat(p1,1,nProj/numel(p1));
    else
        p1 = repelem(p1, prod(n(1:ll-1)));
        p1 = repmat(p1,1,nProj/numel(p1));
    end
    p(ll,:) = p1;
end

indices = zeros(nProj,1);

for r = 1 : nProj
    tt = p(1:numel(nn),r) .* nn;
    indices(r) = sum(tt) + p(end,r) + 1;
end