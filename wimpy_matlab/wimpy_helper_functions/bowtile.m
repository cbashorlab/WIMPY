function [new_seq, right_seq, flip, positions, protiles_F, protiles_R] = bowtile(q,ref,thresh)
% BOWTILE - Aligns query sequences to a reference sequence and determines the best orientation.
%
%   [NEW_SEQ, RIGHT_SEQ, FLIP, POSITIONS, PROTILES_F, PROTILES_R] = BOWTILE(Q, REF, THRESH)
%   aligns the sequences in cell array Q to the reference sequence REF using a threshold THRESH.
%
%   Inputs:
%       Q       - Cell array of query sequences to be aligned.
%       REF     - Reference sequence to which the query sequences are aligned.
%       THRESH  - Threshold value to determine the best orientation of the query sequences.
%
%   Outputs:
%       NEW_SEQ     - Cell array of re-arranged query sequences based on alignment.
%       RIGHT_SEQ   - Cell array of query sequences in the correct orientation.
%       FLIP        - Array indicating the orientation of each query sequence (1 for reverse, 0 for forward, -1 for ambiguous).
%       POSITIONS   - Array of positions where the query sequences align to the reference sequence (forward and reverse).
%       PROTILES_F  - Cell array of positions of forward alignments for each query sequence.
%       PROTILES_R  - Cell array of positions of reverse alignments for each query sequence.

s = ' ';
disp(['Bowtiling to', s, inputname(2)])
warning('off')

if length(ref) > 100
    ref = upper(ref(1:100));
else
    ref = upper(ref);
end

protiles_F = cell(length(q), length(ref) - 10);
protiles_R = cell(length(q), length(ref) - 10);
flip = zeros(length(q), 1);
right_seq = cell(length(q), 1);
new_seq = right_seq;

for i = 1:length(ref)-10
%    disp(i)
    a = ref(i:i+10);
    b = seqrcomplement(a);

    protiles_F(:, i) = strfind(q, a);
    protiles_R(:, i) = strfind(q, b);
end

positions = zeros(length(q), 2); %FWD first, REV second

for i = 1:length(q)
    %disp(i)
    read = cell2mat(q(i));
    x = cell2mat(cat(1, protiles_F(i, :)));
    y = cell2mat(cat(1, protiles_R(i, :)));
    
    if ~isempty(x)
        positions(i, 1) = median(sort(x));
    end
    if ~isempty(y)
        positions(i, 2) = median(sort(y));
    end
    
    px = length(x)/(length(ref)-10);
    py = length(y)/(length(ref)-10);
    
    if py > px + thresh
        flip(i) = 1;
        right_seq(i) = cellstr(seqrcomplement(read));
%        d = length(read) - positions(i, 2);
%        if d > 0
%            a2 = read(d:end);
%            read = strcat(a2, read(1:d-1)); %Re arrange the read so it starts from the alignment point
%            new_seq(i) = cellstr(seqrcomplement(read));
          new_seq(i) = cellstr(seqrcomplement(strcat(read(positions(i, 2):end), read(1:positions(i, 2)))));
%         else
%            new_seq(i) = cellstr('X');
%         end
    elseif px > py + thresh
        right_seq(i) = cellstr(read);
        new_seq(i) = cellstr(strcat(read(positions(i, 1):end), read(1:positions(i, 1)))); %Re-arrange read
    else
        new_seq(i) = cellstr('X');
        flip(i) = -1;
        right_seq(i) = cellstr(read);
    end
    
end

warning('on')