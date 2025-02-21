function [FP_tiles, FP_tiles_nonnorm, conf_chart] = viscount(reads_correct,step,FP, thresh, confchart)
% viscount - Counts the number of tiles identified for each map and read
%
% Syntax: [FP_tiles, FP_tiles_nonnorm, conf_chart] = viscount(reads_correct, step, FP, thresh, confchart)
%
% Inputs:
%    reads_correct - Cell array of sequences to be analyzed
%    step - Integer, the step size for tiling
%    FP - Cell array of feature patterns to be matched
%    thresh - Threshold value for confidence chart
%    confchart - Character, 'T' to generate confidence chart, otherwise no chart
%
% Outputs:
%    FP_tiles - Normalized count of tiles identified for each map and read
%    FP_tiles_nonnorm - Non-normalized count of tiles identified for each map and read
%    conf_chart - Confidence chart matrix if confchart is 'T', otherwise empty

s = ' ';
disp(['Viscounting for', s, inputname(3)])

FP_lengths = zeros(size(FP));

for i = 1:length(FP)
    a = upper(cell2mat(FP(i)));
    FP_lengths(i) = length(a);
    FP(i) = cellstr(a);
end
    
%map_length = min(FP_lengths);

FP_tiles = zeros(length(reads_correct), length(FP));

for j = 1:length(FP)
    a = cell2mat(FP(j));
    %disp(j)
    for i = 1:(length(a)-step)
        %disp(i)
        b = a(i:i+step);
        FP_tiles(:, j) = FP_tiles(:,j) + contains(reads_correct, b);

    end
end

ls = repmat(FP_lengths', length(reads_correct), 1);
FP_tiles_nonnorm = FP_tiles;
FP_tiles = FP_tiles./ls;

conf_chart = zeros(length(FP), length(FP));

if confchart == 'T'
    for i = 1:length(FP)
        for j = 1:length(FP)
            conf_chart(i, j) = sum(FP_tiles(:, i) > thresh & FP_tiles(:, j) > thresh);
        end
    end
end

end

