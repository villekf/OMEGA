function S = assembleS(alkuarvo,T,Ny,Nx,Nz)
% This function creates the necessary weights for the anatomically weighted
% TV prior
% alkuarvo is Ny*Nx*Nz size refimage
% T is the edge threshold parameter
% Ny is the image dimension in y-direction, etc.

S = zeros(Nx*Ny*Nz*3,3,class(alkuarvo));

% compute the weighting gamma and tensor S

f = -diff(alkuarvo,1,2);
f = cat(2,f,zeros(Nx,1,Nz));
f = f(:)';
g = -diff(alkuarvo,1,1);
g = cat(1,g,zeros(1,Ny,Nz));
g = g(:)';
% if Nz > 1
    h = -diff(alkuarvo,1,3);
    h = cat(3,h,zeros(Nx,Ny,1));
    h = h(:)';
    
    gradvec = [f;g;h];
% else
%     gradvec = [f;g];
% end
if exist('vecnorm','builtin') == 5
    gradnorm = vecnorm(gradvec,2,1)';
else
    gradnorm = zeros(size(alkuarvo,1),1);
    for kk = 1 : size(gradvec,2)
        gradnorm(kk) = norm(gradvec(:,kk));
    end
end

gamma = exp(-gradnorm.^2/(T^2));

% -- construct the matrix S -------

for ll = 1: length(gradnorm)
    
    if gradnorm(ll) > 0
        
        nu = gradvec(:,ll)/gradnorm(ll);
        B = eye(3) - (1-gamma(ll))*(nu*nu');
        
    else
        
        B = eye(3);
        
    end
    
    S(3*ll-2:3*ll,:) = B;
    
end
end