function [barcodes1_correct, barcodes2_correct, bc1, bc2] = barcoat(tregion_BFP, method)

if sum(method == 'scar') == 4
    
elseif sum(method == 'align') == 5
    %load barcode structure and BLOSUM62 matrix
    bc1_cons = 'aTTaTTaTTaTTaTTaTTa';
    bc2_cons = 'cTTcTTcTTcTTcTTcTTc';
    
    m = readmatrix('Blosum62_BC1.xlsx'); m = m(:,2:25);
    n = readmatrix('Blosum62_BC2.xlsx'); n = n(:,2:25);
    
    bc1 = cell(size(tregion_BFP)); bc1_pos = zeros(length(tregion_BFP),1);
    bc1_l = bc1_pos; bc1_scores = bc1_pos;
    bc2 = bc1; bc2_pos = bc1_pos;
    bc2_l = bc1_l; bc2_scores = bc1_scores;
    
    for i = 1:length(tregion_BFP)
        a = tregion_BFP{i};
        [score1,barcode1,~] = swalign(bc1_cons, a, 'ScoringMatrix',m);
        bc1_scores(i) = score1;
        bc1(i) = cellstr(barcode1(3,:));
        bc1_l(i) = length(barcode1(3, barcode1(3, :) ~= '-'));
        [score2,barcode2,~] = swalign(bc2_cons, a, 'ScoringMatrix',n);
        bc2_scores(i) = score2;
        bc2(i) = cellstr(barcode2(3,:));
        bc2_l(i) = length(barcode2(3,barcode2(3, :) ~= '-'));
    end
    barcodes1_correct = bc1(bc1_scores > 85 & bc2_scores > 90);
    barcodes2_correct = bc2(bc1_scores > 85 & bc2_scores > 90);

else
    error('Method must be either "scar" or "align"')
end