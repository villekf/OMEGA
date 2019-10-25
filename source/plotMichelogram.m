function plotMichelogram(rings,span, varargin)
%PLOTMICHELOGRAM Plots a Michelogram
%   Plots a Michelogram for the specified ring count, span amount and ring
%   difference. Ring difference is optional (maximum is used if omitted).
%
% Example:
%   plotMichelogram(rings, span)
%   plotMichelogram(rings, span, ring_difference)
%

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

if nargin > 2
    ring_difference = varargin{1};
else
    ring_difference = rings - 1;
end

[x,y] = meshgrid(0.5:1:rings,0.5:1:rings);

erotus = rings - ring_difference - 1;

U = rot90(triu(true(rings,rings),rings - erotus),2);

I = find(U);

for kk = I
    x(kk) = -1;
    y(kk) = -1;
end

U = rot90(triu(true(rings,rings),rings - erotus),4);

I = find(U);

for kk = I
    x(kk) = -1;
    y(kk) = -1;
end

figure

xx = x(x>0);
yy = y(y > 0);

plot(xx,yy,'b.')
hold on

sp1 = floor(span/2);
sp2 = ceil(span/2);
vaihto = 0;
if mod(sp2,2) == 1
    vaihto = 1;
    apu = sp1;
    sp1 = sp2;
    sp2 = apu;
    clear apu
end
ll = sp2;
ii = sp1;
hh = 1;
jj = 1;
minimi = min(sp1,sp2);
for kk = 2 : rings*2 - 2
    if mod(kk,2) == 0
        if kk >= minimi && kk <= rings*2 - minimi
            x1 = [x(ll,hh);x(hh,ll)];
            y1 = [y(ll,hh);y(hh,ll)];
        elseif kk < minimi
            x1 = [x(1, kk);x(kk,1)];
            y1 = [y(1,kk);y(kk,1)];
        elseif kk > rings*2 - minimi
            x1 = [x(rings, kk - rings + 1);x(kk - rings + 1,rings)];
            y1 = [y(rings,kk - rings + 1);y(kk - rings + 1,rings)];
        end
        if kk >= minimi && kk < rings*2 - minimi
            ll = ll + 1;
            hh = hh + 1;
        end
    else
        if kk >= minimi && kk <= rings*2 - minimi
            x1 = [x(ii,jj);x(jj,ii)];
            y1 = [y(ii,jj);y(jj,ii)];
        elseif kk < minimi
            x1 = [x(1, kk);x(kk,1)];
            y1 = [y(1,kk);y(kk,1)];
        elseif kk > rings*2 - minimi
            x1 = [x(rings, kk - rings + 1);x(kk - rings + 1,rings)];
            y1 = [y(rings,kk - rings + 1);y(kk - rings + 1,rings)];
        end
        if kk >= minimi && kk < rings*2 - minimi
            ii = ii + 1;
            jj = jj + 1;
        end
    end
    if sum(x1 < 0) == 0 && sum(y1 < 0) == 0
        line(x1, y1)
    end
end

for gg = 1 : floor((ring_difference + 1)/span)
    lisa = span * (gg - 1);
    maksimi = max(sp1,sp2);
    ll = maksimi+3 + lisa;%span;
    ii = maksimi+2 ++ lisa;%span;
    hh = 2;
    jj = 2;
    pp = 1;
    ss = 1;
    num1 = 1;
    if vaihto == 1
        start = span + 1;
        ii = ii + 1;
        ll = ll - 1;
    else
        start = span;
    end
    uu = ll + 1;
    oo = ii + 1;
    if span == 3
        uu = uu - 1;
        hh = hh - 1;
    end
    over1 = false;
    over2 = false;
    for kk = start : min(rings*2, rings*2)
%         if kk == 8 && gg == 6
%             break
%         end
        if mod(kk,2) == 0
            if kk >= start && ll < rings
                x1 = [y(ll,hh);y(uu,pp)];
                apu_i = ll;
                apu_j = hh;
                apu_o = uu;
                apu_s = pp;
                while(sum(x1 < 0) > 0)
                    i = find(x1 < 0);
                    if i == 1 || length(i) == 2
                        apu_i = apu_i - 1;
                        apu_j = apu_j + 1;
                        x1 = [y(apu_i,apu_j);x1(2)];
                    end
                    if i == 2 || length(i) == 2
                        apu_o = apu_o - 1;
                        apu_s = apu_s + 1;
                        x1 = [x1(1);y(apu_o,apu_s)];
                    end
                end
                y1 = [x(ll,hh);x(uu,pp)];
                apu_i = ll;
                apu_j = hh;
                apu_o = uu;
                apu_s = pp;
                while(sum(y1 < 0) > 0)
                    i = find(y1 < 0);
                    if i == 1 || length(i) == 2
                        apu_i = apu_i - 1;
                        apu_j = apu_j + 1;
                        y1 = [x(apu_i,apu_j);y1(2)];
                    end
                    if i == 2 || length(i) == 2
                        apu_o = apu_o - 1;
                        apu_s = apu_s + 1;
                        y1 = [y1(1);x(apu_o,apu_s)];
                    end
                end
            end
            if kk >= span && ~over2
                hh = hh + 1;
                ll = ll + 1;
                if hh >= maksimi && vaihto == 0 || hh > maksimi && vaihto == 1
                    pp = pp + 1;
                    uu = uu + 1;
                    uu = min(rings, uu);
                else
                    uu = uu + 2;
                end
                if uu >= rings
                    over2 = true;
                    uu = min(rings, uu);
                    if pp == 1 && over1 || pp == 1 && kk == start
                        pp = pp + 1;
                    end
                end
            elseif kk > span && over2
                ll = ll + 1;
                hh = hh + 1;
%                 if pp > 1
                    pp = pp + 2;
%                 else
%                     pp = pp + 1;
%                 end
                uu = min(rings, uu);
            end
        else
            if kk >= start && ii < rings
                x1 = [y(ii,jj);y(oo,ss)];
                apu_i = ii;
                apu_j = jj;
                apu_o = oo;
                apu_s = ss;
                while(sum(x1 < 0) > 0)
                    i = find(x1 < 0);
                    if i == 1 || length(i) == 2
                        apu_i = apu_i - 1;
                        apu_j = apu_j + 1;
                        x1 = [y(apu_i,apu_j);x1(2)];
                    end
                    if i == 2 || length(i) == 2
                        apu_o = apu_o - 1;
                        apu_s = apu_s + 1;
                        x1 = [x1(1);y(apu_o,apu_s)];
                    end
                end
                y1 = [x(ii,jj);x(oo,ss)];
                apu_i = ii;
                apu_j = jj;
                apu_o = oo;
                apu_s = ss;
                while(sum(y1 < 0) > 0)
                    i = find(y1 < 0);
                    if i == 1 || length(i) == 2
                        apu_i = apu_i - 1;
                        apu_j = apu_j + 1;
                        y1 = [x(apu_i,apu_j);y1(2)];
                    end
                    if i == 2 || length(i) == 2
                        apu_o = apu_o - 1;
                        apu_s = apu_s + 1;
                        y1 = [y1(1);x(apu_o,apu_s)];
                    end
                end
            end
            if kk >= span && ~over1
                jj = jj + 1;
                if jj > maksimi && vaihto == 0 || jj >= maksimi && vaihto == 1
                    ii = ii + 1;
                    if kk - span <= minimi && vaihto == 0 && span > 3
                        oo = oo + 2;
                    else
                        oo = oo + 1;
                    end
                    if num1 > minimi || vaihto == 1 || span == 3
                        ss = ss + 1;
                    end
                else
                    ii = ii + 1;
                    oo = oo + 2;
                end
                if oo >= rings
                    over1 = true;
                    oo = min(rings, oo);
                    if ss == 1 && over2 || ss == 1 && kk == start
                        ss = ss + 1;
                    end
                end
                %             num2 = num2 + 1;
            elseif kk > span && over1
                ii = ii + 1;
                jj = jj + 1;
%                 if ss > 1
                    ss = ss + 2;
%                 else
%                     ss = ss + 1;
%                 end
                oo = min(rings, oo);
            end
        end
        if sum(x1 < 0) == 0 && sum(y1 < 0) == 0
            line(x1, y1)
        end
        num1 = num1 + 1;
    end
end

for gg = 1 : floor((ring_difference + 1)/span)
    lisa = span * (gg - 1);
    maksimi = max(sp1,sp2);
    ll = maksimi+3 + lisa;%span;
    ii = maksimi+2 ++ lisa;%span;
    hh = 2;
    jj = 2;
    pp = 1;
    ss = 1;
    num1 = 1;
    if vaihto == 1
        start = span + 1;
        ii = ii + 1;
        ll = ll - 1;
    else
        start = span;
    end
    uu = ll + 1;
    oo = ii + 1;
    if span == 3
        uu = uu - 1;
        hh = hh - 1;
    end
    over1 = false;
    over2 = false;
    for kk = start : min(rings*2, rings*2)
%             if kk == 23 && gg == 6
%                 return
%             end
        if mod(kk,2) == 0
            if kk >= start && ll < rings
                x1 = [x(ll,hh);x(uu,pp)];
                apu_i = ll;
                apu_j = hh;
                apu_o = uu;
                apu_s = pp;
                while(sum(x1 < 0) > 0)
                    i = find(x1 < 0);
                    if i == 1 || length(i) == 2
                        apu_i = apu_i - 1;
                        apu_j = apu_j + 1;
                        x1 = [x(apu_i,apu_j);x1(2)];
                    end
                    if i == 2 || length(i) == 2
                        apu_o = apu_o - 1;
                        apu_s = apu_s + 1;
                        x1 = [x1(1);x(apu_o,apu_s)];
                    end
                end
                y1 = [y(ll,hh);y(uu,pp)];
                apu_i = ll;
                apu_j = hh;
                apu_o = uu;
                apu_s = pp;
                while(sum(y1 < 0) > 0)
                    i = find(y1 < 0);
                    if i == 1 || length(i) == 2
                        apu_i = apu_i - 1;
                        apu_j = apu_j + 1;
                        y1 = [y(apu_i,apu_j);y1(2)];
                    end
                    if i == 2 || length(i) == 2
                        apu_o = apu_o - 1;
                        apu_s = apu_s + 1;
                        y1 = [y1(1);y(apu_o,apu_s)];
                    end
                end
            end
            if kk >= span && ~over2
                hh = hh + 1;
                ll = ll + 1;
                if hh >= maksimi && vaihto == 0 || hh > maksimi && vaihto == 1
                    pp = pp + 1;
                    uu = uu + 1;
                    uu = min(rings, uu);
                else
                    uu = uu + 2;
                end
                if uu >= rings
                    over2 = true;
                    uu = min(rings, uu);
                    if pp == 1 && over1 || pp == 1 && kk == start
                        pp = pp + 1;
                    end
                end
            elseif kk > span && over2
                ll = ll + 1;
                hh = hh + 1;
                pp = pp + 2;
                uu = min(rings, uu);
            end
        else
            if kk >= start && ii < rings
                y1 = [y(ii,jj);y(oo,ss)];
                apu_i = ii;
                apu_j = jj;
                apu_o = oo;
                apu_s = ss;
                while(sum(y1 < 0) > 0)
                    i = find(y1 < 0);
                    if i == 1 || length(i) == 2
                        apu_i = apu_i - 1;
                        apu_j = apu_j + 1;
                        y1 = [y(apu_i,apu_j);y1(2)];
                    end
                    if i == 2 || length(i) == 2
                        apu_o = apu_o - 1;
                        apu_s = apu_s + 1;
                        y1 = [y1(1);y(apu_o,apu_s)];
                    end
                end
                x1 = [x(ii,jj);x(oo,ss)];
                apu_i = ii;
                apu_j = jj;
                apu_o = oo;
                apu_s = ss;
                while(sum(x1 < 0) > 0)
                    i = find(x1 < 0);
                    if i == 1 || length(i) == 2
                        apu_i = apu_i - 1;
                        apu_j = apu_j + 1;
                        x1 = [x(apu_i,apu_j);x1(2)];
                    end
                    if i == 2 || length(i) == 2
                        apu_o = apu_o - 1;
                        apu_s = apu_s + 1;
                        x1 = [x1(1);x(apu_o,apu_s)];
                    end
                end
            end
            if kk >= span && ~over1
                jj = jj + 1;
                if jj > maksimi && vaihto == 0 || jj >= maksimi && vaihto == 1
                    ii = ii + 1;
                    if kk - span <= minimi && vaihto == 0 && span > 3
                        oo = oo + 2;
                    else
                        oo = oo + 1;
                    end
                    if num1 > minimi || vaihto == 1 || span == 3
                        ss = ss + 1;
                    end
                else
                    ii = ii + 1;
                    oo = oo + 2;
                end
                if oo >= rings
                    over1 = true;
                    oo = min(rings, oo);
                    if ss == 1 && over2 || ss == 1 && kk == start
                        ss = ss + 1;
                    end
                end
                %             num2 = num2 + 1;
            elseif kk > span && over1
                ii = ii + 1;
                jj = jj + 1;
                ss = ss + 2;
                oo = min(rings, oo);
            end
        end
        if sum(x1 < 0) == 0 && sum(y1 < 0) == 0
            line(x1, y1)
        end
        num1 = num1 + 1;
    end
end

for gg = 1 : min(round((ring_difference + 1)/span),round(rings/span))
    lisa = span * (gg - 1);
    line([maksimi - 0.5 + lisa; rings], [0;rings-maksimi + 0.5 - lisa])
    line([0;rings-maksimi + 0.5 - lisa], [maksimi - 0.5 + lisa; rings])
end

hold off
xlim([0,rings])
ylim([0,rings])