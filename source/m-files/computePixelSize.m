function [xx,yy,zz,dx,dy,dz,bx,by,bz] = computePixelSize(FOV, N, offset, cType)
%COMPUTEPIXELSIZE Computes the pixel size and distance from origin
%   Utility function for OMEGA

% Pixel boundaries
etaisyys = -(FOV) / 2;
dx = zeros(1,size(FOV,2));
dy = zeros(1,size(FOV,2));
dz = zeros(1,size(FOV,2));
bx = zeros(1,size(FOV,2));
by = zeros(1,size(FOV,2));
bz = zeros(1,size(FOV,2));
for kk = size(FOV,2) : - 1 : 1
    xx = double(linspace(etaisyys(1,kk) + offset(1), -etaisyys(1,kk) + offset(1), N(1,kk) + 1));
    yy = double(linspace(etaisyys(2,kk) + offset(2), -etaisyys(2,kk) + offset(2), N(2,kk) + 1));
    zz = double(linspace(etaisyys(3,kk) + offset(3), -etaisyys(3,kk) + offset(3), N(3,kk) + 1));

    % Distance of adjacent pixels
    dx(kk) = diff(xx(1:2));
    dy(kk) = diff(yy(1:2));
    dz(kk) = diff(zz(1:2));

    % Distance of image from the origin
    if kk == 1
        bx(kk) = xx(1,1);
        by(kk) = yy(1,1);
        bz(kk) = zz(1,1);
    else
        if kk > 5 || (size(FOV,2) == 5 && kk > 3)
            if mod(kk,2) == 1
                by(kk) = offset(2) + FOV(2,1) / 2;
            else
                by(kk) = offset(2) - FOV(2,1) / 2 - FOV(2,kk);
            end
            bx(kk) = xx(1,1);
            bz(kk) = zz(1,1);
        elseif (kk > 3 && kk < 6) || (size(FOV,2) == 5 && kk > 1)
            if mod(kk,2) == 1
                bx(kk) = etaisyys(1,1) + offset(1) + FOV(1,1);
            else
                bx(kk) = etaisyys(1,1) + offset(1) - FOV(1,kk);
            end
            by(kk) = yy(1,1);
            bz(kk) = zz(1,1);
        elseif kk > 1 && kk < 4
            bx(kk) = xx(1,1);
            by(kk) = yy(1,1);
            if mod(kk,2) == 1
                bz(kk) = etaisyys(3,1) + offset(3) + FOV(3,1);
            else
                bz(kk) = etaisyys(3,1) + offset(3) - FOV(3,kk);
            end
        end
    end
    xx = cast(xx, cType);
    yy = cast(yy, cType);
    zz = cast(zz, cType);
    dx = cast(dx, cType);
    dy = cast(dy, cType);
    dz = cast(dz, cType);
    bx = cast(bx, cType);
    by = cast(by, cType);
    bz = cast(bz, cType);
end
end
