function grad = TGV(im,maxits,alpha,beta, Nx, Ny, Nz)
%TGV Total Generalized Variation prior (TGV)
%
% Example:
%   grad = TGV(im, maxits, alpha, beta, Nx, Ny, Nz)
% INPUTS:
%   im = The current estimate
%   maxits = Number of TGV iterations
%   alpha = First regularization value
%   beta = Second regularization value
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%
% OUTPUTS:
%   grad = The gradient of TGV
% 
% See also MRP, FMH, L_filter, Quadratic_prior, TVpriorfinal, Weighted_mean

% Based on code by
% Christian Clason (christian.clason@uni-graz.at)
% Kristian Bredies (kristian.bredies@uni-graz.at)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im = reshape(im, Nx, Ny, Nz);

[n,m,c] = size(im);

proj   = @(x1,x2,x3) max(1,sqrt(x1.^2+x2.^2+x3.^2)/beta);  % projection on C_beta
proj2D   = @(x1,x2) max(1,sqrt(x1.^2+x2.^2)/beta);  % projection on C_beta
proj2  = @(x1,x2,x3,x4) max(1,sqrt(x1.^2+x2.^2+x3.^2+3*x4.^2)/alpha); % on C^2_alpha
proj22D  = @(x1,x2,x4) max(1,sqrt(x1.^2+x2.^2+2*x4.^2)/alpha); % on C^2_alpha
dx     = @(x) [diff(x,1,1); x(1,:,:)-x(n,:,:)];       % discrete x-derivative
dy     = @(x) [diff(x,1,2), x(:,1,:)-x(:,m,:)];       % discrete y-derivative
dz     = @(x) cat(3,diff(x,1,3), x(:,:,1)-x(:,:,c));       % discrete z-derivative
dxt    = @(x) [x(n,:,:)-x(1,:,:); -diff(x,1,1)];      % transpose x-derivative
dyt    = @(x) [x(:,m,:)-x(:,1,:), -diff(x,1,2)];      % transpose y-derivative
dzt     = @(x) cat(3,x(:,:,c)-x(:,:,1), -diff(x,1,3));       % discrete z-derivative

sigma  = 1/16;    % dual step size
tau    = 1/8;    % primal step size

% initialize iterates p=(p1,p2), q=(q1,q2,q3); grad, v1, v2
p1 = zeros(n,m,c);
p2 = zeros(n,m,c);
p3 = zeros(n,m,c);
q1 = zeros(n,m,c);
q2 = zeros(n,m,c);
q3 = zeros(n,m,c);
q4 = zeros(n,m,c);

grad = zeros(n,m,c);
v1 = zeros(n,m,c);
v2 = zeros(n,m,c);
v3 = zeros(n,m,c);

% initialize leading points  \bar grad, \bar v = (vb1,vb2)
ub  = grad;
vb1 = v1; vb2 = v2; vb3 = v3;

for iter = 1:maxits
    
    % update dual (p): (eta1,eta2) = nabla(u+grad-u0)-v,
    ukp = ub + im;
    eta1 = dx(ukp) - vb1;
    eta2 = dy(ukp) - vb2;
    if Nz > 1
        eta3 = dz(ukp) - vb3;
    end
    y1 = p1 + sigma * eta1;
    y2 = p2 + sigma * eta2;
    if Nz > 1
        y3 = p3 + sigma * eta3;
        my = proj(y1,y2,y3);
    else
        my = proj2D(y1,y2);
    end
    p1 = y1./my;
    p2 = y2./my;
    if Nz > 1
        p3 = y3./my;
    end
    
    % update dual (q): (zeta1,zeta2,zeta3) = E(v); v: (zeta4,zeta5) = -p - div q
    zeta1 = dx(vb1);
    zeta2 = dy(vb2);
    if Nz > 1
        zeta3 = dz(vb3);
        zeta4 = (dy(vb1) + dy(vb3) + dx(vb2) + dx(vb3) + dz(vb1) + dz(vb2))/6;
    else
        zeta4 = (dy(vb1) + dx(vb2))/2;
    end
    z1 = q1 + sigma * zeta1;
    z2 = q2 + sigma * zeta2;
    z4 = q4 + sigma * zeta4;
    if Nz > 1
        z3 = q3 + sigma * zeta3;
        mz = proj2(z1,z2,z3, z4);
    else
        mz = proj22D(z1,z2,z4);
    end
    q1 = z1./mz;
    q2 = z2./mz;
    if Nz > 1
        q3 = z3./mz;
    end
    q4 = z4./mz;
    
    % update primal (grad): F'*F grad  + F(u) - div(p)
    if Nz > 1
        eta4 = dxt(p1)+dyt(p2)+dzt(p3); % primal variable: image
    else
        eta4 = dxt(p1)+dyt(p2);
    end
    uold = grad;
    grad   = grad - tau * eta4;
    
    % update primal (v): (zeta4,zeta5) = -p - div q
    if Nz > 1
        zeta5 = dxt(q1) + dyt(q4) + dzt(q4) - p1;
        zeta6 = dxt(q4) + dyt(q2) + dzt(q4) - p2;
        zeta7 = dxt(q4) + dyt(q4) + dzt(q3) - p3;
    else
        zeta5 = dxt(q1) + dyt(q4) - p1;
        zeta6 = dxt(q4) + dyt(q2) - p2;
    end
    v1old = v1;
    v2old = v2;
    if Nz > 1
        v3old = v3;
    end
    v1 = v1 - tau * zeta5;
    v2 = v2 - tau * zeta6;
    if Nz > 1
        v3 = v3 - tau * zeta7;
    end
    
    % update primal leading points
    ub = 2*grad - uold;
    vb1 = 2*v1 - v1old;
    vb2 = 2*v2 - v2old;
    if Nz > 1
        vb3 = 2*v3 - v3old;
    end
end
grad = -grad(:);
