function [found, positions, protiles_F] = tilepin(q, ref, thresh, suppressmessage)
% WELCOME to TilePin
%   Detailed explanation goes here
s = ' ';
if suppressmessage ~= 'T'
   disp(['Tilepinning for', s, inputname(2)])
end
warning('off')

ref = upper(ref);

found = zeros(length(q), 1);
protiles_F = cell(length(q), length(ref) - 10);

for i = 1:length(ref)-10
    %disp(i)
    a = ref(i:i+10);
    protiles_F(:, i) = strfind(q, a);
end



positions = zeros(length(q), 1);

for i = 1:length(q)
    x = cell2mat(cat(1, protiles_F(i, :)));
    found(i) = length(x);
    if length(x)/length(ref) > thresh
        positions(i) = median(sort(x));
    else
        positions(i) = -1;
    end
end

warning('on')



