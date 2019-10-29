function [grad] = TVpriorFinal(im,data,Nx, Ny ,Nz,anatomical, options, TVtype, tr_offsets)
%TVPRIORFINAL Total Variation prior (TV) or Asymmetric Parallel Level Sets (APLS) prior
%
% Example:
%   grad = TVpriorFinal(im, data, Nx, Ny, Nz, anatomical, options, TVtype)
% INPUTS:
%   im = The current estimate
%   data = Anatomical weights for TV type 1
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   anatomical = Whether anatomical weighting is used
%   options = Weights T and C, reference image for APLS (APLS_ref_image),
%   smoothing parameter for APLS (APLSsmoothing), eta parameter for APLS
%   (eta)
%   TVtype = Type of anatomical weighting to use, 1-3 = Different TVs, 4 =
%   APLS
%
% OUTPUTS:
%   grad = The gradient of TV/APLS
%
% See also MRP, FMH, L_filter, Quadratic_prior, AD, TGV, Weighted_mean

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

if TVtype ~= 3
    
    im = reshape(im,Nx,Ny,Nz);
    f = (zeros(Nx,Ny,Nz));
    f(1:Nx-1,:,:) = -diff(im);
    f(end,:,:) = im(end,:,:) - im(1,:,:);
    g = (zeros(Nx,Ny,Nz));
    g(:,1:Ny-1,:) = -diff(im,1,2);
    g(:,end,:) = im(:,end,:) - im(:,1,:);
    h = (zeros(Nx,Ny,Nz));
    h(:,:,1:Nz-1) = -diff(im,1,3);
    h(:,:,end) = im(:,:,end) - im(:,:,1);
    
    f = f(:);
    g = g(:);
    h = h(:);
    
    if anatomical
        
        if TVtype == 1
            s1 = data.s1;
            s2 = data.s2;
            s3 = data.s3;
            s4 = data.s4;
            s5 = data.s5;
            s6 = data.s6;
            s7 = data.s7;
            s8 = data.s8;
            s9 = data.s9;
            pval = sqrt(s1.*f.^2 + s5.*g.^2 + s9.*h.^2 + s4.*f.*g + s7.*f.*h + s2.*f.*g + s8.*h.*g + s3.*f.*h + s6.*h.*g + data.beta);
            apu1 = 0.5 * (2*s1.*f + s4.*g + s7.*h + s2.*g + s3.*h)./pval; % x-direction
            apu2 = 0.5 * (2*s5.*g + s4.*f + s2.*f + s8.*h + s6.*h)./pval; % y-direction
            apu3 = 0.5 * (2*s9.*h + s8.*g + s6.*g + s7.*f + s3.*f)./pval; % z-direction
            apu4 = 0.5 * (2*s1.*f + 2*s5.*g + 2*s9.*h + + s4.*f + s2.*f + s8.*h + s6.*h + s4.*g + s7.*h + s2.*g + s3.*h + s8.*g + s6.*g + s7.*f + s3.*f)./pval;
        end
        if TVtype == 2
            fp = (zeros(Nx,Ny,Nz));
            fp(1:Nx-1,:,:) = -diff(data.reference_image);
            fp(end,:,:) = data.reference_image(end,:,:) - data.reference_image(1,:,:);
            gp = (zeros(Nx,Ny,Nz));
            gp(:,1:Ny-1,:) = -diff(data.reference_image,1,2);
            gp(:,end,:) = data.reference_image(:,end,:) - data.reference_image(:,1,:);
            hp = (zeros(Nx,Ny,Nz));
            hp(:,:,1:Nz-1) = -diff(data.reference_image,1,3);
            hp(:,:,end) = data.reference_image(:,:,end) - data.reference_image(:,:,1);
            
            fp = fp(:);
            gp = gp(:);
            hp = hp(:);
            
            pval = sqrt(f.^2 + g.^2 + h.^2 + options.T*(fp.^2 + gp.^2 + hp.^2) + data.beta);
            apu1 = (f)./pval; % x-direction
            apu2 = (g)./pval; % y-direction
            apu3 = (h)./pval; % z-direction
            apu4 = ( f + g + h)./pval;
        end
        if TVtype == 4
            fp = (zeros(Nx,Ny,Nz));
            fp(1:Nx-1,:,:) = -diff(options.APLS_ref_image);
            fp(end,:,:) = options.APLS_ref_image(end,:,:) - options.APLS_ref_image(1,:,:);
            gp = (zeros(Nx,Ny,Nz));
            gp(:,1:Ny-1,:) = -diff(options.APLS_ref_image,1,2);
            gp(:,end,:) = options.APLS_ref_image(:,end,:) - options.APLS_ref_image(:,1,:);
            hp = (zeros(Nx,Ny,Nz));
            hp(:,:,1:Nz-1) = -diff(options.APLS_ref_image,1,3);
            hp(:,:,end) = options.APLS_ref_image(:,:,end) - options.APLS_ref_image(:,:,1);
            
            fp = fp(:);
            gp = gp(:);
            hp = hp(:);
            
            epsilon = [fp,gp,hp]./[fp./sqrt(fp.^2 + options.eta^2) + options.epps, gp./sqrt(gp.^2 + options.eta^2) + options.epps, hp./sqrt(hp.^2 + options.eta^2) + options.epps];
            
            pval = (f.^2 + g.^2 + h.^2 - (sum([f,g,h].*epsilon,2)).^2 + options.APLSsmoothing);
            pval(pval < 0) = options.APLSsmoothing;
            pval = sqrt(pval);
            apu1 = (f - ((sum([f,g,h].*epsilon,2)).*((epsilon(:,1)))))./pval; % x-direction
            apu2 = (g - ((sum([f,g,h].*epsilon,2)).*((epsilon(:,2)))))./pval; % y-direction
            apu3 = (h - ((sum([f,g,h].*epsilon,2)).*((epsilon(:,3)))))./pval; % z-direction
            apu4 = (g - ((sum([f,g,h].*epsilon,2)).*((epsilon(:,2)))) + f - ((sum([f,g,h].*epsilon,2)).*((epsilon(:,1)))) + ...
                h - ((sum([f,g,h].*epsilon,2)).*((epsilon(:,3)))))./pval;
        end
    else
        pval = sqrt(f.^2 + g.^2 + h.^2 + data.beta);
        apu1 = f./pval; % x-direction
        apu2 = g./pval; % y-direction
        apu3 = h./pval; % z-direction
        apu4 = ( f + g + h)./pval;
    end
    apu1=reshape(apu1,Nx,Ny,Nz);
    apu2=reshape(apu2,Nx,Ny,Nz);
    apu3=reshape(apu3,Nx,Ny,Nz);
    apu4=reshape(apu4,Nx,Ny,Nz);
    % apu1=[apu1(end,:,:);apu1(1:end-1,:,:)];
    % apu2=[apu2(:,end,:),apu2(:,1:end-1,:)];
    % apu3=cat(3,apu3(:,:,end),apu3(:,:,1:end-1));
    apu1 = circshift(apu1,1,1);
    apu2 = circshift(apu2,1,2);
    apu3 = circshift(apu3,1,3);
    grad = apu4 - apu1 - apu2 - apu3;
    
    grad = grad(:);
    extgrad = 2*options.tau*min(im(:));
    grad = grad + extgrad;
else
    if anatomical
        padd = padding(reshape(im,Nx,Ny,Nz),[options.Ndx options.Ndy options.Ndz]);
        padd = padd(tr_offsets);
        padd(:,isinf(options.weights)) = [];
        padd2 = padding(reshape(data.reference_image,Nx,Ny,Nz),[options.Ndx options.Ndy options.Ndz]);
        padd2 = padd2(tr_offsets);
        padd2(:,isinf(options.weights)) = [];
        grad = sum(bsxfun(@times,(bsxfun(@minus,im, padd)./options.C^2 .* (1./sqrt(1 + (bsxfun(@minus,im, padd)./options.C).^2 + (bsxfun(@minus,data.reference_image, padd2)./options.T).^2))),options.weights_quad(~isinf(options.weights_quad)')),2);
    else
        padd = padding(reshape(im,Nx,Ny,Nz),[options.Ndx options.Ndy options.Ndz]);
        padd = padd(tr_offsets);
        padd(:,isinf(options.weights)) = [];
        grad = sum(bsxfun(@times,(bsxfun(@minus,im, padd)./options.C^2 .* (1./sqrt(1 + (bsxfun(@minus,im, padd)./options.C).^2))),options.weights_quad(~isinf(options.weights_quad))'),2);
    end
end