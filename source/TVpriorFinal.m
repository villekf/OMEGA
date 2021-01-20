function [grad, varargout] = TVpriorFinal(im,aData,Nx, Ny ,Nz,anatomical, options, TVtype, varargin)
%TVPRIORFINAL Total Variation prior (TV) or Asymmetric Parallel Level Sets (APLS) prior
%   Computes the gradient of the specified TV prior or APLS prior. The
%   gradients are the same in TVtype = 1,2 and APLS if no anatomical prior
%   is used. TVtype = 3 is not technically a TV prior, but is very similar
%   and is such grouped in the same function. APLS is enabled by setting
%   TVtype = 5. SATV is obtained with TVtype = 4. The MRF-function can be
%   obtained as the second output optionally (grad is the gradient of
%   this).
%
% Examples:
%   grad = TVpriorFinal(im, data, Nx, Ny, Nz, anatomical, options, TVtype)
%   grad = TVpriorFinal(im, data, Nx, Ny, Nz, anatomical, options, TVtype, tr_offsets)
% INPUTS:
%   im = The current estimate
%   aData = Anatomical weights for TV type 1
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   anatomical = Whether anatomical weighting is used
%   options = Weights T and C, reference image for APLS (APLS_ref_image),
%   smoothing parameter for APLS (APLSsmoothing), eta parameter for APLS
%   (eta), tau parameter for types 1, 2, and APLS, TVsmoothing
%   TVtype = Type of anatomical weighting to use, 1-3 = Different TVs, 4 =
%   APLS
%   tr_offsets = Needed ONLY for TVtype = 3, can be omitted otherwise.
%
% OUTPUTS:
%   grad = The gradient of TV/APLS
%
% See also MRP, FMH, L_filter, Quadratic_prior, AD, TGV, Weighted_mean

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

if TVtype ~= 3
    
    im = reshape(im,Nx,Ny,Nz);
    f = (zeros(Nx,Ny,Nz));
    f(1:Nx-1,:,:) = -diff(im);
%     f(end,:,:) = im(end,:,:) - im(1,:,:);
    f(end,:,:) = -f(Nx-1,:,:);
    g = (zeros(Nx,Ny,Nz));
    g(:,1:Ny-1,:) = -diff(im,1,2);
%     g(:,end,:) = im(:,end,:) - im(:,1,:);
    g(:,end,:) = -g(:,Ny-1,:);
    h = (zeros(Nx,Ny,Nz));
    h(:,:,1:Nz-1) = -diff(im,1,3);
%     h(:,:,end) = im(:,:,end) - im(:,:,1);
    h(:,:,end) = -h(:,:,Nz-1);
    
    f = f(:);
    g = g(:);
    h = h(:);
%     
%     if TVtype == 4
%         t = f + g + h;
%         s = sign(t);
%         grad = s - s ./ (abs(t)/options.SATVPhi + 1);
%         if nargout >= 2
%             varargout{1} = abs(t) - options.SATVPhi .* log(1 + abs(t)./options.SATVPhi);
%         end
%     else
        
        if anatomical || TVtype == 5
            
            if TVtype == 1
                s1 = aData.s1;
                s2 = aData.s2;
                s3 = aData.s3;
                s4 = aData.s4;
                s5 = aData.s5;
                s6 = aData.s6;
                s7 = aData.s7;
                s8 = aData.s8;
                s9 = aData.s9;
                pval = sqrt(s1.*f.^2 + s5.*g.^2 + s9.*h.^2 + s4.*f.*g + s7.*f.*h + s2.*f.*g + s8.*h.*g + s3.*f.*h + s6.*h.*g + options.TVsmoothing);
                apu1 = 0.5 * (2*s1.*f + s4.*g + s7.*h + s2.*g + s3.*h)./pval; % x-direction
                apu2 = 0.5 * (2*s5.*g + s4.*f + s2.*f + s8.*h + s6.*h)./pval; % y-direction
                apu3 = 0.5 * (2*s9.*h + s8.*g + s6.*g + s7.*f + s3.*f)./pval; % z-direction
                apu4 = 0.5 * (2*s1.*f + 2*s5.*g + 2*s9.*h + + s4.*f + s2.*f + s8.*h + s6.*h + s4.*g + s7.*h + s2.*g + s3.*h + s8.*g + s6.*g + s7.*f + s3.*f)./pval;
            end
            if TVtype == 2
                fp = (zeros(Nx,Ny,Nz));
                fp(1:Nx-1,:,:) = -diff(aData.reference_image);
                % fp(end,:,:) = aData.reference_image(end,:,:) - aData.reference_image(1,:,:);
                fp(end,:,:) = -fp(Nx-1,:,:);
                gp = (zeros(Nx,Ny,Nz));
                gp(:,1:Ny-1,:) = -diff(aData.reference_image,1,2);
                % gp(:,end,:) = aData.reference_image(:,end,:) - aData.reference_image(:,1,:);
                gp(:,end,:) = -gp(:,Ny-1,:);
                hp = (zeros(Nx,Ny,Nz));
                hp(:,:,1:Nz-1) = -diff(aData.reference_image,1,3);
                % hp(:,:,end) = aData.reference_image(:,:,end) - aData.reference_image(:,:,1);
                hp(:,:,end) = -hp(:,:,Nz-1);
                
                fp = fp(:);
                gp = gp(:);
                hp = hp(:);
                
                pval = sqrt(f.^2 + g.^2 + h.^2 + options.T*(fp.^2 + gp.^2 + hp.^2) + options.TVsmoothing);
                apu1 = (f)./pval; % x-direction
                apu2 = (g)./pval; % y-direction
                apu3 = (h)./pval; % z-direction
                apu4 = ( f + g + h)./pval;
            end
            if TVtype == 5
                fp = (zeros(Nx,Ny,Nz));
                fp(1:Nx-1,:,:) = -diff(options.APLS_ref_image);
%                 fp(end,:,:) = options.APLS_ref_image(end,:,:) - options.APLS_ref_image(1,:,:);
                fp(end,:,:) = -fp(Nx-1,:,:);
                gp = (zeros(Nx,Ny,Nz));
                gp(:,1:Ny-1,:) = -diff(options.APLS_ref_image,1,2);
%                 gp(:,end,:) = options.APLS_ref_image(:,end,:) - options.APLS_ref_image(:,1,:);
                gp(:,end,:) = -gp(:,Ny-1,:);
                hp = (zeros(Nx,Ny,Nz));
                hp(:,:,1:Nz-1) = -diff(options.APLS_ref_image,1,3);
%                 hp(:,:,end) = options.APLS_ref_image(:,:,end) - options.APLS_ref_image(:,:,1);
                hp(:,:,end) = -hp(:,:,Nz-1);
                
                fp = fp(:) + options.epps;
                gp = gp(:) + options.epps;
                hp = hp(:) + options.epps;
                
%                 epsilon = [fp,gp,hp]./[sqrt(fp.^2 + options.eta^2)+ sqrt(gp.^2 + options.eta^2)+ sqrt(hp.^2 + options.eta^2)];
                epsilon = bsxfun(@rdivide,[fp,gp,hp], sqrt(fp.^2 + gp.^2 + hp.^2 + options.eta^2));
                
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
            if TVtype == 4
                if options.SATVPhi == 0
                    apu1 = sign(f); % x-direction
                    apu2 = sign(g); % y-direction
                    apu3 = sign(h); % z-direction
                else
                    apu1 = sign(f) - sign(f) ./ (abs(f)/options.SATVPhi + 1); % x-direction
                    apu2 = sign(g) - sign(g) ./ (abs(g)/options.SATVPhi + 1); % y-direction
                    apu3 = sign(h) - sign(h) ./ (abs(h)/options.SATVPhi + 1); % z-direction
                end
            else
                pval = sqrt(f.^2 + g.^2 + h.^2 + options.TVsmoothing);
                apu1 = f./pval; % x-direction
                apu2 = g./pval; % y-direction
                apu3 = h./pval; % z-direction
            end
            apu4 = apu1 + apu2 + apu3;
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
        if nargout >= 2
            varargout{1} = pval;
        end
%     end
else
    if isempty(varargin)
        error('Offset-values (tr_offsets) need to be input when using TVtype = 3!')
    end
    if anatomical
        padd = padding(reshape(im,Nx,Ny,Nz),[options.Ndx options.Ndy options.Ndz]);
        padd = padd(varargin{1});
        padd(:,isinf(options.weights)) = [];
        padd2 = padding(reshape(aData.reference_image,Nx,Ny,Nz),[options.Ndx options.Ndy options.Ndz]);
        padd2 = padd2(varargin{1});
        padd2(:,isinf(options.weights)) = [];
        grad = sum(bsxfun(@times,(bsxfun(@minus,im, padd)./options.C^2 .* (1./sqrt(1 + (bsxfun(@minus,im, padd)./options.C).^2 + ...
            (bsxfun(@minus,aData.reference_image, padd2)./options.T).^2))),options.weights_quad(~isinf(options.weights_quad)')),2);
        if nargout >= 2
            varargout{1} = sum(bsxfun(@times, sqrt(1 + (bsxfun(@minus,im, padd)./options.C).^2 + (bsxfun(@minus,aData.reference_image, padd2)./options.T).^2) - 1, options.weights_quad(~isinf(options.weights_quad)')),2);
        end
    else
        padd = padding(reshape(im,Nx,Ny,Nz),[options.Ndx options.Ndy options.Ndz]);
        padd = padd(varargin{1});
        padd(:,isinf(options.weights)) = [];
        grad = sum(bsxfun(@times,(bsxfun(@minus,im, padd)./options.C^2 .* (1./sqrt(1 + (bsxfun(@minus,im, padd)./options.C).^2))),options.weights_quad([1:floor(numel(options.weights_quad)/2),ceil(numel(options.weights_quad)/2) + 1 : numel(options.weights_quad)])),2);
        if nargout >= 2
            varargout{1} = sum(bsxfun(@times, sqrt(1 + (bsxfun(@minus,im, padd)./options.C).^2) - 1, options.weights_quad([1:floor(numel(options.weights_quad)/2),ceil(numel(options.weights_quad)/2) + 1 : numel(options.weights_quad)])),2);
        end
    end
end