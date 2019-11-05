function [index, pituus] = goldenAngleSubsets(options)
%GOLDENANGLESUBSETS Chooses the subsets based on golden angle sampling
%   Detailed explanation goes here

[~, ~, kulma] = subset_angles(options);
kulma = reshape(kulma, options.Ndist, options.Nang);

ga = 180/((1+sqrt(5))/2);

koko = ceil(options.Nang/options.subsets);
if mod(options.Nang,options.subsets) > 0
    last = options.Nang - (options.subsets - 1) * koko;
else
    last = koko;
end

loc = diff(kulma);
loc = find(abs(loc(1,:)) > mean(abs(loc(1,:))),1);

kulma(:,loc) = max(kulma(:,loc));

keskiarvo = mean(kulma);
orig_keskiarvo = keskiarvo;

angle0 = 0;
indices = zeros(koko, options.subsets);
for kk = 1 : options.subsets
    if kk < options.subsets
        for ll = 1 : koko
            [~,closestIndex] = min(abs(bsxfun(@minus,angle0, keskiarvo')));
            closestIndex2 = find(orig_keskiarvo == keskiarvo(closestIndex));
            indices(ll,kk) = closestIndex2;
            keskiarvo(closestIndex) = [];
            angle0 = angle0 + ga;
            if angle0 >= 180
                angle0 = angle0 - 180;
            end
        end
    else
        for ll = 1 : last
            [~,closestIndex] = min(abs(bsxfun(@minus,angle0, keskiarvo')));
            closestIndex2 = find(orig_keskiarvo == keskiarvo(closestIndex));
            indices(ll,kk) = closestIndex2;
            keskiarvo(closestIndex) = [];
            angle0 = angle0 + ga;
            if angle0 >= 180
                angle0 = angle0 - 180;
            end
        end
    end
end

index1 = uint32(1:options.Ndist*options.Nang*options.NSinos)';
index1 = reshape(index1, options.Ndist, options.Nang, options.NSinos);

index = cell(options.subsets,1);
pituus = zeros(options.subsets,1,'uint32');
for i=1:options.subsets
    index3 = [];
    for ll = 1 : size(indices,1)
        if indices(ll,i) == 0
            break
        end
        index2 = squeeze(index1(:,indices(ll,i),:));
        index3 = [index3; index2(:)];
    end
    index{i} = index3(:);
    pituus(i) = uint32(length(index{i}));
end
