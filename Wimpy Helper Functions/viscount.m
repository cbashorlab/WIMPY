function [FP_tiles, FP_tiles_nonnorm, conf_chart] = viscount(reads_correct,step,FP, thresh, confchart)
%UNTITLED Summary of this function goes here
%   Inputs = sequence maps for each part
%   Outputs = number of tiles identified for each map and read

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

