function sys = makeOCT(charArray)
%MAKEOCT Builds Octave mex-files
%   If the default method is unable to build the files due to spaces, this
%   function should correctly build the files

ind = strfind(charArray,'-L');
for ll = 1 : length(ind) - 1
    charArray = [charArray(1:ind(1) - 1), '"-', charArray(ind(1) + 1:end)];
    ind2 = strfind(charArray,'-l');
    ind2(ind2 < min(ind)) = [];
    if isempty(ind2)
        ind2 = inf;
    end
    ind = strfind(charArray,' -L');
    indi = min(ind(1), ind2(1));
    charArray = [charArray(1:indi - 1), '" ' charArray(indi:end)];
    ind = strfind(charArray,' -L') + 1;
end
tps = strfind(charArray, 'g++');
tps = [tps, numel(charArray) + 1];
for ll = 1 : length(tps) - 1
    sys = system(charArray(tps(ll):tps(ll + 1) -1));
end
end

