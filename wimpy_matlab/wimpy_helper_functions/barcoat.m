function [bc, bc_l, bc_score, bc_pos] = barcoat(tregion_BFP, preset, scoring_matrix)
% barcoat - Aligns barcodes using semi-constrained structures.
%
%
% Inputs:
%    tregion_BFP - Cell array of sequences to be analyzed.
%    preset = either 'BBA', 'DDC' or 'None' (BBA and DDCs are currently
%    supported with predefined custom scoring matrices. 
%    scoring_matrix - Required if using a preset other than 'BBA' or 'DDC'.
%    scoring_matrix must be defined as a 4X4 matrix representing match
%    scores for 'A', 'C', 'G', 'T' in that order.
%
% Outputs:
%    barcodes1_correct - Cell array of corrected barcodes for BC1.
%    barcodes2_correct - Cell array of corrected barcodes for BC2.
%    bc1 - Cell array of aligned barcodes for BC1.
%    bc2 - Cell array of aligned barcodes for BC2.
%
% Example:
%    [barcodes1_correct, barcodes2_correct, bc1, bc2] = barcoat(tregion_BFP, 'align')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: Blosum62_BC1.xlsx, Blosum62_BC2.xlsx
%
% See also: swalign
if sum(preset == 'BBA') == 3
   cus = 'ATTATTATTATTATTATTA';
   m = [5 -5 -5 -5; -5 5 5 5; -5 5 5 5; -5 5 5 5]; %A, C, G, T
elseif sum(preset == 'DDC')
   cus = 'CTTCTTCTTCTTCTTCTTC';
   m = [5 -5 5 5; -5 5 -5 -5; 5 -5 5 5; 5 -5 5 5];
else
    m = scoring_matrix;
end

if sum(size(m) == [4, 4]) ~= 2
    error('scoring_matrix must be a 4X4 matrix if presets are not BBA or DDC')
end


bc = cell(size(tregion_BFP)); bc_pos = zeros(length(tregion_BFP),1);
bc_l = bc_pos; bc_score = bc_pos;

for i = 1:length(tregion_BFP)
    a = cell2mat(tregion_BFP(i));
    if contains(a, 'X')
        bc_score(i) = 0; bc(i) = cellstr('X'); bc_pos(i) = -1; bc_l(i) = 0;
    else
        [score, barcode, pos] = swalign(cus, a, 'alphabet', 'nt', 'ScoringMatrix', m);
        bc_score(i) = score;
        bc(i) = cellstr(barcode(3,:));
        bc_pos(i) = pos(2);
        bc_l(i) = length(barcode(3, barcode(3, :) ~= '-'));
    end
end

