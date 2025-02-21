%% Analysing A subset of data from Experiment 68 (166K-member Level 3 library)

addpath('./wimpy_helper_functions');
%Load reference sequences
Puro = upper('cgctccgcatcggcctaaggaaccggcgtggttcctggctacggtgggagtctcacctgaccatcaaggaaagggattgggaagtgctgtcgttcttcca');
GFP = upper('atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa');
A4 = upper('ggagcatcgcccttccccggccctcaggtaagaggaccaaataccgtagccgtttccaatttcagtcctttagcgccacctggtgctaactactctatcacgcttttatccaataactacctttgtaaatttcctttcaaaagttctggccgggcgcggtggcTGTAC');
A4 = A4(end-50:end);
minP = upper('gctcactctcccttacacggagtggataTAGT');
mRuby = upper('GTGAGTAAAGGCGAAGAACTTATCAAGGAAAATATGCGGATGAAAGTGGTTATGGAGGGTAGCGTGAACGGACACCAGTTCAAATGCACGGGAGAGGGCGAGGGGCGACCCTACGAGGGAGTCCAAACAATGAGGATTAAGGTTATAGAAGGTGGTCCGCTGCCATTCGCATTCGATATTTTGGCCACGTCCTTCATGTACGGCTCCCGAACCTTTATCAAATACCCTGCGGATATCCCAGACTTTTTCAAGCAATCCTTTCCGGAAGGGTTCACGTGGGAGCGAGTCACGAGATATGAGGATGGAGGCGTAGTAACAGTAACCCAAGACACATCACTTGAGGACGGTGAGCTTGTCTACAATGTGAAGGTACGCGGCGTCAATTTCCCCTCAAATGGCCCGGTGATGCAAAAGAAAACTAAAGGATGGGAGCCCAACACCGAAATGATGTACCCGGCAGATGGGGGGCTTAGGGGCTATACGGACATCGCATTGAAGGTTGATGGCGGGGGCCATCTCCATTGTAACTTTGTAACTACATATCGGTCAAAAAAGACTGTGGGGAACATTAAAATGCCGGGAGTACACGCTGTTGATCATCGCCTGGAAAGGATAGAGGAAAGCGACAATGAAACGTATGTAGTACAGCGGGAGGTCGCCGTCGCCAAATATAGTAATCTGGGCGGTGGCATGGACGAGCTTTATAAA');
mRuby = mRuby(200:300);
BFP = upper('GTCAGTAAAGGGGAAgagcttataaaggaaaatatgcacatgaagctctacatggagggcactgtagataaccaccatttcaaatgtacctctgaaggggagggcaagccatacgaaggtactcaaaccatgcgaataaaagtagttgaaggcgggcctcttccctttgcattcgacattctcgcaacctcttttctgtacggcagtaagactttcataaaccacactcaaggcattccagacttcttcaagcaatcattccccgaaggattcacctgggagcgagttactacttatgaggacggaggagtccttactgcaacccaagacacctcactgcaagatgggtgcctgatttacaatgtaaagatcagaggggtgaatttcacaagcaatgggccagttatgcaaaaaaagacccttggatgggaggccttcaccgagacactgtacccagccgatggtggactggagggcaggaatgacatggccctcaagctcgtcggaggcagtcacttgattgccaacgccaaaaccacttaccgctctaagaaacctgctaaaaacctgaagatgcccggcgtctattatgtggactatcgacttgagagaattaaggaggcaaacaacgagacttatgtcgaacagcatgaagttgccgtggctaggtactgtgatttgcccagcaaattgggtcataaacttaac');
BFP = BFP(200:300);
minP_100k = readcell('../example_ref_sequences/100k_Minimal-Promoters.xlsx'); minP_100k = minP_100k(2:end, 2);
terminators_100k = readcell('../example_ref_sequences/100k_Terminators.xlsx'); terminators_100k = terminators_100k(2:end, 2);
spacers_100k = readcell('../example_ref_sequences/100k_Terminator_spacers.xlsx'); spacers_100k = spacers_100k(2:end, 2);
promoters_100k = readcell('../example_ref_sequences/100k-Promoters.xlsx'); promoters_100k = promoters_100k(2:end, 2);
ORF_parts = readcell('../example_ref_sequences/100k_ORF-Parts_SynTF.xlsx'); AD_parts = ORF_parts(2:5, 3);
WT = upper('tcagaatcagggcatctcaaacgacatctccgc');
C1_WT = upper('ccaggagaaCGAccatttca tctgtatgCGTaatttctct atttgcatgAGGaacttttc atttgcatgAGGaattttag atctgtatgAGAaatttctc atatgtatgCGcaacttttc');
C1_6x = upper('CCAGGAGAAGCCCCATTTCA TCTGTATGGCCAATTTCTCT ATTTGCATGGCCAACTTTTC ATTTGCATGGCCAATTTTAG ATCTGTATGGCTAATTTCTC ATATGTATGGCCAACTTTTC');
ZF_parts = {WT C1_WT C1_6x};

%load in data using fastqall function
[~, l, seq] = fastqall('../example_fastq', 'fastq');
%Filter based on read length
seq = seq(l > 9500 & l < 15000); l = l(l > 9500 & l < 15000);

%Index to Puromycin using bowtile
thresh = 0.03;
[new_seq, ~, ~, ~, ~, ~] = bowtile(seq, Puro, thresh);
%Remove any reads that didn't align to Puro (marked as 'X' in new_seq)
reads_correct = new_seq(~contains(new_seq, 'X')); l_readscorrect = l(~contains(new_seq, 'X'));

%% REPORTER Identification
% Locate A4, GFP, constant region upstream of pMin, BFP and mRuby via tilepin
[~, positionsA4, ~] = tilepin(reads_correct, A4, thresh, 'F');
[~, positionsGFP, ~] = tilepin(reads_correct, GFP(end-100:end), thresh, 'F');
[~, positionsminP, ~] = tilepin(reads_correct, minP, thresh, 'F');
[~, positionsmRuby, ~] = tilepin(reads_correct, mRuby(end-100:end), thresh, 'F');
[~, positionsBFP, ~] = tilepin(reads_correct, BFP, thresh, 'F');

positions2 = floor([positionsA4, positionsGFP, positionsminP, positionsmRuby, positionsBFP]);

%Assessing Reporter diversity
%pregions and tregions using chophat
pregions = chophat(reads_correct, positions2(:, 1:2), 0, 0);
tregions = chophat(reads_correct, positions2(:, 2), 2000, 1);

%Minimal Promoter
[minP_variants_scaled, minP_variants, minPconf] = viscount(pregions, 6, minP_100k, 0.2, 'T');

[~,minP_variants_scaled(:,4)] = max(minP_variants_scaled, [], 2);
minP_variants_scaled((sum(minP_variants_scaled(:,1:3),2) < 0.2),4) = 0;

%Number of Binding sites using FASTar
BS10_1 = upper('cGGCGTAGCCGATGTCGCGc');
[variants_bs, ~] = FASTar(pregions, BS10_1, 6, 8);

%Everything between
variants_bs(variants_bs > 6.9 & variants_bs < 9.1) = 8;
variants_bs(variants_bs > 9.2 & variants_bs < 14) = 12;
variants_bs(variants_bs > 3.9 & variants_bs < 6.1) = 4;
variants_bs(variants_bs ~= 2 & variants_bs ~= 4 & variants_bs ~= 8 & variants_bs ~= 12) = 0;

variants_bs(variants_bs == 2) = 1; %2 binding sites
variants_bs(variants_bs == 4) = 2; %4 binding sites
variants_bs(variants_bs == 8) = 3; %8 binding sites
variants_bs(variants_bs == 12) = 4; %12 binding sites

%Terminators
[~, variants_term, ~] = viscount(tregions, 10, terminators_100k, 0.2, 'F');

%scale number of tiles to length of reference
term_lengths = [length(cell2mat(terminators_100k(1))), length(cell2mat(terminators_100k(2))), length(cell2mat(terminators_100k(3)))];
variants_term_scaled = zeros(size(variants_term));
variants_term_scaled(:,1) = variants_term(:,1);
variants_term_scaled(:,2) = variants_term(:,2)-variants_term(:,1).*1.5;
variants_term_scaled(:,3) = variants_term(:,3)-variants_term(:,2);
threshold = 10; variants_term_scaled(variants_term_scaled < threshold) = 0;

[~,term_indices] = max(variants_term_scaled, [], 2);
x = sum(variants_term_scaled,2); term_indices(x < 1) = 0;
%Calculate assignments: minP_variants_scaled, variants_BS, term_indices
reporter_variants = [variants_bs minP_variants_scaled(:,4) term_indices];

%% SynTF Identification
%Promoters
pregions_synTF = chophat(reads_correct, [ones(size(positions2, 1), 1) positions2(:, 4)], 0, 0);
[synTFprom_variants_scaled, ~, synTFprom_conf] = viscount(pregions_synTF,10,promoters_100k, 0.03, 'T');
[~,variants_synTF_prom] = max(synTFprom_variants_scaled'); variants_synTF_prom = variants_synTF_prom';
variants_synTF_prom(sum(synTFprom_variants_scaled,2) < 0.03) = 0;

% ORF Parts
[~, f, AD_conf] = viscount(pregions_synTF, 10, AD_parts, 0.2, 'T');
%match mRuby index to AD indices
[f5, ~, p5] = tilepin(pregions_synTF, mRuby(1:50), 1, 'F');

AD_conf(4, 4) = sum(f(:, 4) > 10 & f(:, 3) < 100);
AD_conf(3, 3) = sum(f(:, 3) > 100);
AD_conf(3, 4) = sum(f(:, 4) < 10 & f(:, 3) > 100);
AD_conf(4, 3) = AD_conf(3, 4);

%Assigning AD variants
f(f < 30) = 0; [~,variants_AD] = max(f, [], 2); variants_AD(sum(f,2) == 0) = -1;

%Identify IDRs in synTF_pregion

IDR_parts = ORF_parts(6:9, 3);
[IDR_tiles, ~, IDR_conf] = viscount(pregions_synTF, 10, IDR_parts, 0.15, 'T');
[~, variants_IDR] = max(IDR_tiles, [], 2);

%Identify Zinc Fingers
tregions_synTF = chophat(reads_correct, [positions2(:, 4) positions2(:, 1)], 0, 0);
[ZF_tiles, ~, ZF_conf] = viscount(tregions_synTF, 10, ZF_parts', 0.03, 'T');

[~,variants_ZF] = max(ZF_tiles, [], 2);
variants_ZF(sum(ZF_tiles,2) == 0) = -1;

%Looking for terminators
length_term = 100;
%assign terminators
[synTF_term, ~, term_conf] = viscount(tregions_synTF, 10, terminators_100k, 0.4, 'T');
[~,variants_synTF_term] = max(synTF_term, [], 2);

%assign spacers
[spacer_synTF, ~, ~] = viscount(tregions_synTF, 10, spacers_100k, 0.4, 'F');
spacer_synTF(spacer_synTF < 0.2) = 0;
spacer_synTF = [zeros(size(spacer_synTF, 1), 1) spacer_synTF];
[~, variants_synTF_spacer] = max(spacer_synTF, [], 2);

%Calculate assignments
synTF_variants = [variants_synTF_prom variants_AD variants_IDR variants_ZF variants_synTF_term variants_synTF_spacer];

%% Calculate orientation

%based on location of GFP and mRuby
variants_orientation = ones(size(reads_correct, 1), 1);
variants_orientation(positions2(:, 2) < positions2(:, 4)) = 2;

%% Sum all synTF and reporter indices to create full list of 100k identified circuits
all_100k_variants = [synTF_variants reporter_variants variants_orientation];
%remove all rows containing zero or -1
all_100k_variants(any(all_100k_variants < 1,2),:) = [];
%library_100k = all_100k_variants(:,1) + 4*(all_100k_variants(:,2)-1) + 16*(all_100k_variants(:,3)-1) + 80*(all_100k_variants(:,4)-1) + 240*(all_100k_variants(:,5)-1) + 2880*(all_100k_variants(:,6)-1) + 11520*(all_100k_variants(:,7)-1) +  34560*(all_100k_variants(:,8)-1);
library_100k = (all_100k_variants - [0 ones(1, 9)]).*[1 4 16 64 192 768 2304 9216 27648 82944];

%% Barcoding & Barcoded AssTable Generation
tregion_BFP = chophat(reads_correct, positions2(:, 5), 1000, 1);
%load barcode structure and BLOSUM62 matrix
bc1_cons = 'aTTaTTaTTaTTaTTaTTa';
bc2_cons = 'cTTcTTcTTcTTcTTcTTc';

m = readmatrix('Blosum62_BC1.xlsx');
m = m(:,2:25);

n = readmatrix('Blosum62_BC2.xlsx');
n = n(:,2:25);

bc1 = cell(size(tregion_BFP));
bc1_pos = zeros(length(tregion_BFP),1);
bc1_l = bc1_pos;
bc1_scores = bc1_pos;

bc2 = bc1;
bc2_pos = bc1_pos;
bc2_l = bc1_l;
bc2_scores = bc1_scores;

count = 0;

for i = 1:length(tregion_BFP)
    a = tregion_BFP{i};
    [score1,barcode1,bc1_p] = swalign(bc1_cons, a, 'ScoringMatrix',m);
    bc1_scores(i) = score1;
    bc1(i) = cellstr(barcode1(3,:));
    bc1_pos(i) = bc1_p(2);
    bc1_l(i) = length(barcode1(3,:));
    [score2,barcode2,bc2_p] = swalign(bc2_cons, a, 'ScoringMatrix',n);
    bc2_scores(i) = score2;
    bc2(i) = cellstr(barcode2(3,:));
    bc2_pos(i) = bc2_p(2);
    bc2_l(i) = length(barcode2(3,:));    
    count = count + logical(score1 > 85 & score2 > 90);
end

barcodes1_correct = bc1(bc1_scores > 85 & bc2_scores > 90);
barcodes2_correct = bc2(bc1_scores > 85 & bc2_scores > 90);
barcoded_variants = all_100k_variants(bc1_scores > 85 & bc2_scores > 90,:);

%Remove all reads with incomplete circuit IDs
barcodes1_correct(any(barcoded_variants < 1,2)) = [];
barcodes2_correct(any(barcoded_variants < 1,2)) = [];
barcoded_variants(any(barcoded_variants < 1,2),:) = [];

%Making the ass table
AB_100k_asstable = zeros(length(barcoded_variants),1);
for i = 1:length(barcoded_variants)
    AB_100k_asstable(i) = sum(barcoded_variants(i,1) + 4*(barcoded_variants(i,2)-1) + 16*(barcoded_variants(i,3)-1) + 64*(barcoded_variants(i,4)-1) + 192*(barcoded_variants(i,5)-1) + 2304*(barcoded_variants(i,6)-1) + 9216*(barcoded_variants(i,7)-1) +  27486*(barcoded_variants(i,8)-1));
end

%% Barcodes

barcode_regions = cell(size(reads_correct));
for i = 1:length(reads_correct)
    
    disp(i)
        a = cell2mat(reads_correct(i));
        b = floor(positions2(i,:));
        
        if b(5) + 1000 < length(a) && b(5) > 0
            barcode_regions(i) = cellstr(a(b(5):b(5) + 1000));
        elseif b(5) + 1000 > length(a) && b(5) < length(a) && b(5) > 0
            barcode_regions(i) = cellstr(a(b(5):end));
        else
            barcode_regions(i) = cellstr('X');
        end
end

upscar_bc1 = strfind(barcode_regions, 'GAAACG');
downscar_bc1 = strfind(barcode_regions, 'ACAGT');

BC1 = cell(size(barcode_regions));
BC1_lengths = zeros(size(BC1));

for i = 1:length(BC1)
    disp(i)
    y = cell2mat(barcode_regions(i));
    a = cell2mat(upscar_bc1(i));
    if ~isempty(a)
        a = a(1);
    end
    b = cell2mat(downscar_bc1(i));
    
    if ~isempty(a)
        if length(b) > 1
            for j = 1:length(b)
                x = b(j);
                if x < a + 50 && x > a
                    b = x;
                    break
                end
            end
        else
            BC1(i) = cellstr('X');
            BC1_lengths(i) = -1000;
        end

        if length(a) == 1 && length(b) == 1
            BC1(i) = cellstr(y((a + 6):(b - 1)));
            BC1_lengths(i) = b - a - 6;
        else
            BC1(i) = cellstr('X');
            BC1_lengths(i) = -1000;
        end
    else
        BC1(i) = cellstr('X');
        BC1_lengths(i) = -1000;
    end
    
end
% 
% barcodes_correct = BC1(BC1_lengths == 17);
% [bc1_unique,~,b] = unique(barcodes_correct);
% bc_count = zeros(size(bc1_unique));
% 
% for i = 1:length(bc1_unique)    
%     bc_count(i) = sum(b == i);
% end

upscar_bc2 = strfind(barcode_regions, 'AGGTAC');
downscar_bc2 = strfind(barcode_regions, 'CGATA');

BC2 = cell(size(barcode_regions));
BC2_lengths = zeros(size(BC2));

for i = 1:length(BC2)
    disp(i)
    y = cell2mat(barcode_regions(i));
    a = cell2mat(upscar_bc2(i));
    if ~isempty(a)
        a = a(end);
    end
    b = cell2mat(downscar_bc2(i));
    
    if ~isempty(a)
        if length(b) > 1
            for j = 1:length(b)
                x = b(j);
                if x < a + 50 && x > a
                    b = x;
                    break
                end
            end
        else
            BC2(i) = cellstr('X');
            BC2_lengths(i) = -1000;
        end

        if length(a) == 1 && length(b) == 1
            BC2(i) = cellstr(y((a + 6):(b - 1)));
            BC2_lengths(i) = b - a - 6;
        else
            BC2(i) = cellstr('X');
            BC2_lengths(i) = -1000;
        end
    else
        BC2(i) = cellstr('X');
        BC2_lengths(i) = -1000;
    end
    
end


    
