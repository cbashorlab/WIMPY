function [nums, locs] = fastar(pregions, ref, step, bw)
% fastar - Finds approximate substring matches in regions of interest
%
% Syntax:  [nums, locs] = fastar(pregions, ref, step, bw)
%
% Inputs:
%    pregions - Cell array of strings, each representing a region of interest
%    ref - Reference string to search for substrings
%    step - Length of the substring to search for
%    bw - Bandwidth for kernel density estimation
%
% Outputs:
%    nums - Array of integers, each representing the number of approximate matches found in the corresponding region
%    locs - Cell array of arrays, each containing the locations of the approximate matches in the corresponding region
%
% Example:
%    pregions = {'ATCG', 'GCTA', 'CGAT'};
%    ref = 'ATCG';
%    step = 2;
%    bw = 0.5;
%    [nums, locs] = fastar(pregions, ref, step, bw);
%

nums = zeros(size(pregions));
locs = cell(size(pregions));

for j = 1:length(pregions)
    x = cell2mat(pregions(j));
    a1 = [];
    if contains(x, 'X')
        continue
    end

    for i = 1:length(ref) - step
        a = ref(i:i+step);
        y = strfind(x, a);
        if ~isempty(y)
            a1(end+1:end+length(y)) = y;
        end
    end

    if ~isempty(a1)
        [f, x] = ksdensity(unique(a1), 'bandwidth', bw);
        a1 = x(islocalmax(f));
        locs(j) = {a1};
        nums(j) = length(a1);
    end
end
