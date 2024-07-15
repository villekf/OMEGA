function alpha = lfilter_weights(Ndx, Ndy, Ndz, dx, dy, dz, oned_weights)
%LFILTER_WEIGHTS Computes the (Laplace distributed) weights for the
%L-filter
%
% The weights are either 1D (oned_weights = true) or 2D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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
N = (Ndx*2+1)*(Ndy*2+1)*(Ndz*2+1);

% 1D weights
if oned_weights
    % Precomputed weights
    if N == 3
        alpha = [0.15168;0.69663];
        alpha = [alpha;flip(alpha(1:end-1))];
    elseif N == 9
        alpha = [-0.01899;0.02904;0.06965;0.23795;0.36469];
        alpha = [alpha;flip(alpha(1:end-1))];
    elseif N == 25
        alpha = [0.0055;0.00335;-0.00427;-0.00101;-0.00008;0.00065;0.00314;0.01064;0.02907;0.06499;0.11835;0.17195;0.19541];
        alpha = [alpha;flip(alpha(1:end-1))];
    % Compute the weights
    elseif N < 67
        b = 1/sqrt(2);
        CorrM = zeros(N,N);
        for i = 1 : ceil(N/2)
            K = factorial(N)/(factorial(i-1)*factorial(N-i));
            f1 = @(x) x.^2.*(0.5*exp(x/b)).^(i-1).*(1 - 0.5*exp(x/b)).^(N-i).*(1/(2*b)*exp(x/b));
            f2 = @(x) x.^2.*(1 - 0.5*exp(-x/b)).^(i-1).*(1 - (1 - 0.5*exp(-x/b))).^(N-i).*(1/(2*b)*exp(-x/b));
            CorrM(i,i) = K*(integral(f1,-Inf,0) + integral(f2,0,Inf));
            %         f1 = @(x) x.^2.*(0.5.*(1+sign(x).*(1-exp(-abs(x)/b)))).^(i-1).*(1-(0.5.*(1+sign(x).*(1-exp(-abs(x)/b))))).^(N-i).*(1/(2*b).*exp(-abs(x)/b));
            %         CorrM(i,i) = K*integral(f1,-Inf,Inf);
            for j = i+1 : N-i+1
                K = factorial(N)/(factorial(i-1)*factorial(j-i-1)*factorial(N-j));
                f1 = @(x,y) x.*y.*(0.5*exp(x/b)).^(i-1).*(0.5*exp(y/b) - 0.5*exp(x/b)).^(j-i-1).*(1 - 0.5*exp(y/b)).^(N-j).*(1/(2*b)*exp(x/b)).*(1/(2*b)*exp(y/b));
                f2 = @(x,y) x.*y.*(1-0.5*exp(-x/b)).^(i-1).*(1 - 0.5*exp(-y/b) - (1 - 0.5*exp(-x/b))).^(j-i-1).*(1 - (1 - 0.5*exp(-y/b))).^(N-j).*(1/(2*b)*exp(-x/b)).*(1/(2*b)*exp(-y/b));
                CorrM(i,j) = K*(integral2(f1,-Inf,0,-Inf,0) + integral2(f2,0,Inf,0,Inf));
                %             f2 = @(x,y) x.*y.*(0.5.*(1+sign(x).*(1-exp(-abs(x)/b)))).^(i-1).*((0.5.*(1+sign(y).*(1-exp(-abs(y)/b)))) - (0.5.*(1+sign(x).*(1-exp(-abs(x)/b))))).^(j-i-1)...
                %                 .*(1- (0.5.*(1+sign(y).*(1-exp(-abs(y)/b))))).^(N-j).*(1/(2*b).*exp(-abs(x)/b)).*(1/(2*b).*exp(-abs(y)/b));
                %             CorrM(i,j) = K*integral2(f2,-Inf,Inf,-Inf,Inf);
            end
        end
        CorrM = CorrM + triu(CorrM,1)';
        temp = rot90(CorrM,2);
        temp(N:N-1:end-1) = 0;
        CorrM = CorrM + temp;
        
        e = ones(N,1);
        alpha = (CorrM\e)/(e'/CorrM*e);
    else
        b = 1/sqrt(2);
        x = linspace(-4,0,ceil(N/2))';
        alpha = 0.5*exp(x/b);
        alpha = [alpha;flip(alpha(1:end-1))];
        alpha = alpha/sum(alpha);
    end
else
    % 2D weights
    dxy = double(sqrt(dx^2 + dy^2));
    dxz = double(sqrt(dx^2 + dz^2));
    dyz = double(sqrt(dy^2 + dz^2));
    dxyz = double(sqrt(sqrt(dy^2 + dx^2) + dz^2));
    dmin = double(min([dx, dy, dz]));
    b = 1/sqrt(2);
    x = linspace(-1.5,0,min(1,Ndx) + min(1,Ndy) + 1)';
    alpha1 = 0.5*exp(x/b);
    if Ndz > 0
        ll = 1;
        alpha = zeros((Ndx*2+1)*(Ndy*2+1)*(Ndz*2+1),1);
        for lz = -int32(Ndz) : int32(Ndz)
            for ly = -int32(Ndy) : int32(Ndy)
                for lx = -int32(Ndx) : int32(Ndx)
                    if ly == 0 && lz == 0 && lx ~= 0
                        alpha(ll) = alpha1(2) / (double(dx) * abs(lx));
                    elseif ly == 0 && lz == 0 && lx == 0
                        alpha(ll) = alpha1(end) / (dmin);
                    elseif lx == 0 && lz == 0
                        alpha(ll) = alpha1(2) / (double(dy) * abs(ly));
                    elseif lx == 0 && ly == 0
                        alpha(ll) = alpha1(2) / (double(dz) * abs(lz));
                    elseif ly ~= 0 && lz == 0 && lx ~= 0
                        alpha(ll) = alpha1(1) / (dxy);
                    elseif ly ~= 0 && lz ~= 0 && lx == 0
                        alpha(ll) = alpha1(1) / (dyz);
                    elseif ly == 0 && lz ~= 0 && lx ~= 0
                        alpha(ll) = alpha1(1) / (dxz);
                    elseif ly ~= 0 && lz ~= 0 && lx ~= 0
                        alpha(ll) = alpha1(1) / (dxyz);
                    end
                    ll = ll + 1;
                end
            end
        end
    else
        ll = 1;
        alpha = zeros((Ndx*2+1)*(Ndy*2+1),1);
        for ly = -int32(Ndy) : int32(Ndy)
            for lx = -int32(Ndx) : int32(Ndx)
                if ly == 0 && lx ~= 0
                    alpha(ll) = alpha1(2) / (double(dx) * abs(lx));
                elseif ly == 0 && lx == 0
                    alpha(ll) = alpha1(end) / (dmin);
                elseif lx == 0
                    alpha(ll) = alpha1(2) / (double(dy) * abs(ly));
                elseif ly ~= 0 && lx ~= 0
                    alpha(ll) = alpha1(1) / (dxy);
                end
                ll = ll + 1;
            end
        end
    end
    alpha = alpha/sum(alpha);
end