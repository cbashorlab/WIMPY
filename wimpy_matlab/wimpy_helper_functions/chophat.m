function out = chophat(reads_correct, positions, truncationcriteria, retain)
% CHOPHAT - Trims sequences based on given positions and criteria.
% 
% out = CHOPHAT(reads_correct, positions, truncationcriteria, retain)
% trims the sequences in reads_correct based on the positions and
% truncation criteria provided. The function returns a cell array of
% trimmed sequences.
%
% Inputs:
%     reads_correct - Cell array of sequences to be trimmed.
%     positions - Array of positions indicating where to start and end
%                 trimming for each sequence. Can be a scalar or a two-element array.
%     truncationcriteria - Scalar indicating the maximum length of the
%                          trimmed sequence. If the remaining sequence
%                          length is greater than this value, it will be
%                          truncated.
%     retain - Scalar (0 or 1). If 1, retains the sequence from the
%              starting position to the end if truncationcriteria is not met.
%
% Outputs:
%     out - Cell array of trimmed sequences.
%
% Example:
%     reads_correct = {'ATCGTACGATCG', 'GCTAGCTAGCTA'};
%     positions = [3; 5];
%     truncationcriteria = 4;
%     retain = 1;
%     out = chophat(reads_correct, positions, truncationcriteria, retain);
%     % out will be {'CGTA', 'AGCT'}
    
out = cell(size(reads_correct));

%Trim reads into regions
for i = 1:length(reads_correct)
    a = cell2mat(reads_correct(i));
    b = positions(i,:);
    if isscalar(b)
        if b > 0
            if truncationcriteria > 0
                if b > 0 && length(a(b:end)) > truncationcriteria
                    out(i) = cellstr(a(b:b+truncationcriteria));
                elseif b > 0 && retain == 1
                    out(i) = cellstr(a(b:end));
                else
                    out(i) = cellstr('X');
                end
            else
                out(i) = cellstr(a(b:end));
            end
        else
            out(i) = cellstr('X');
        end
    else
        if b(1) > 0 && b(2) > b(1)
            out(i) = cellstr(a(b(1):b(2)));
        else
            out(i) = cellstr('X');
        end
    end
end

