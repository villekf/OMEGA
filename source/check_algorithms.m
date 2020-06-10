function img = check_algorithms(pz, algorithms, color_from_algo)
%CHECK_ALGORITHMS Check if the input algorithm exists
%   This function checks the input cell matrix and whether or not the input
%   algorithm actually exists there. If not, the available algorithms are
%   displayed.

algo_char = algorithms_char();
if color_from_algo == 0
    img = pz{algorithms(1)};
else
    img = pz{algorithms(color_from_algo)};
end
for jj = 1:numel(algorithms)
    if isempty(pz{algorithms(jj)})
        warning('The current selected algorithm does not contain any estimes!')
        fprintf('The following are contained in the input array:\n')
        indeksi = ~cellfun(@isempty,pz);
        indeksi = indeksi(:,1);
        char_ar = algo_char(indeksi);
        loc = find(indeksi);
        for kk = 1 : nnz(indeksi)-1
            fprintf('Element %d: %s\n', loc(kk), char_ar{kk})
        end
        img = [];
        return
    end
end
end

