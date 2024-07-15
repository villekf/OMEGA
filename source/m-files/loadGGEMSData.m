function options = loadGGEMSData(options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
list = dir([options.fpath '/*.mhd']);
options.nProjections = numel(list);
apu = loadMetaImage([options.fpath '/' list(1).name]);
options.SinM = zeros(size(apu,1), size(apu,2),options.nProjections,class(apu));
options.SinM(:,:,1) = apu;
options.nRowsD = size(apu,1);
options.nColsD = size(apu,2);
nimi = list(1).name;
nimi = nimi(1:end-5);
nimi = numel(nimi);
eka = list(1).name;
eka(1:nimi) = [];
eka = str2double(eka(1));
for kk = 2 : options.nProjections
    apu = loadMetaImage([options.fpath '/' list(kk).name]);
    proj = list(kk).name;
    proj((1:nimi)) = [];
    if isnan(str2double(proj(1:4)))
        if isnan(str2double(proj(1:3)))
            if isnan(str2double(proj(1:2)))
                u = 1;
            else
                u = 2;
            end
        else
            u = 3;
        end
    else
        u = 4;
    end
    proj = str2double(proj(1:u));
    if eka == 1
        proj = (proj);
    else
        proj = (proj) + 1;
    end
    options.SinM(:,:,proj) = apu;
end
d = dir(options.fpath);
isub = [d(:).isdir];
folders = {d(isub).name}';
airFound = strcmpi(folders,'air');
if any(airFound)
    path = [options.fpath '/' cell2mat(folders(airFound))];
    list = dir([path '/*.mhd']);
    air = zeros(size(options.SinM),class(options.SinM));
    nimi = list(1).name;
    nimi = nimi(1:end-5);
    nimi = numel(nimi);
    eka = list(1).name;
    eka(1:nimi) = [];
    eka = str2double(eka(1));
    for kk = 1 : options.nProjections
        apu = loadMetaImage([path '/' list(kk).name]);
        proj = list(kk).name;
        proj((1:nimi)) = [];
        if isnan(str2double(proj(1:4)))
            if isnan(str2double(proj(1:3)))
                if isnan(str2double(proj(1:2)))
                    u = 1;
                else
                    u = 2;
                end
            else
                u = 3;
            end
        else
            u = 4;
        end
        proj = str2double(proj(1:u));
        if eka == 1
            proj = (proj);
        else
            proj = (proj) + 1;
        end
        if kk == 1
            maksimi = apu(round(size(apu,1)/2),round(size(apu,2)/2));
        end
        % if proj == 53
        %     return
        % end
        apu = double(apu) ./ double(maksimi);
        air(:,:,proj) = apu;
        options.SinM(:,:,proj) = cast(double(options.SinM(:,:,proj)) ./ double(apu), class(options.SinM));
    end
    options.flat = maksimi;
else
    options.flat = max(options.SinM(:));
end