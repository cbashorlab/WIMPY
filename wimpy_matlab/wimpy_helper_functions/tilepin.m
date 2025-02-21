function [found, positions, protiles_F] = tilepin(q, ref, thresh, suppressmessage)
% WELCOME to TilePin
%   Detailed explanation goes here

% TILEPIN finds occurrences of substrings in a reference string.
%
%   [FOUND, POSITIONS, PROTILES_F] = TILEPIN(Q, REF, THRESH, SUPPRESSMESSAGE)
%   searches for substrings of length 11 in the reference string REF within
%   the query cell array Q. It returns the following:
%       FOUND - A binary array indicating whether each element of Q has
%               substrings found in REF.
%       POSITIONS - An array indicating the median position of the found
%                   substrings in REF for each element of Q. If the number
%                   of found substrings divided by the length of REF is
%                   less than THRESH, the position is set to -1.
%       PROTILES_F - A cell array containing the positions of found
%                    substrings for each element of Q.
%
%   Inputs:
%       Q - A cell array of query strings.
%       REF - A reference string in which to search for substrings.
%       THRESH - A threshold value to determine if the found substrings are
%                significant.
%       SUPPRESSMESSAGE - A flag to suppress the display message. If set to
%                         'T', the message is suppressed.
%
%   Example:
%       q = {'ATCG', 'GCTA', 'TACG'};
%       ref = 'ATCGTACGATCG';
%       thresh = 0.1;
%       suppressmessage = 'F';
%       [found, positions, protiles_F] = tilepin(q, ref, thresh, suppressmessage);
%
%   See also STRFIND, MEDIAN, SORT.

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



