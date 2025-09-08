function MUvol = HU_to_mu(CTvol, keV)
% Values from https://doi.org/10.1016/j.apradiso.2008.01.002
    energies = [75 80 140 159 167 171 245 364 511];
    intercept = 0.1 * [0.16 0.16 0.15 0.14 0.13 0.14 0.12 0.09 0.09];
    slope1 = 0.1 * 1e-4 * [1.66 1.66 1.52 1.38 1.36 1.39 1.25 0.941 0.96];
    slope2 = 0.1 * 1e-4 * [1.48 1.48 1.14 1.26 1.25 1.05 0.861 0.559 0.573];

    slope1_interp = interp1(energies, slope1, keV);
    slope2_interp = interp1(energies, slope2, keV);
    intercept_interp = interp1(energies, intercept, keV);
    %error("break")

    %slope1_interp = 0.15*1e-3;
    %slope2_interp = 5.59*0.15*1e-4;
    %intercept_interp = 0.15;
%
    %slope1 = (0.045-0.0) / (2000-(-1000));
    %intercept1 = 0.015; % water
    %MUvol = (CTvol * slope1 + intercept1);

    MUvol = zeros(size(CTvol));
    MUvol(CTvol < 0) = slope1_interp * CTvol(CTvol < 0) + intercept_interp;
    MUvol(CTvol >= 0) = slope2_interp * CTvol(CTvol >= 0) + intercept_interp;
    %MUvol = 0.1 * MUvol;

    

    %MUvol = max(MUvol, 0);
end