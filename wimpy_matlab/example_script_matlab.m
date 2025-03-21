%% Analysing A subset of data from O'Connell, Rai et al
%This is a plain text version of example_live_script_matlab.mlx. For more
%context behind the functions used in this script, and detailed
%explanations on how the WIMPY pipeline works in matlab, please refer to
%the live script file.

addpath('./wimpy_helper_functions');
%Load reference sequences
ref_seqs = struct2cell(fastaread('../example_ref_sequences/ref_sequences.fasta'));
Puro = upper(cell2mat(ref_seqs(2, contains(ref_seqs(1, :), 'Puro'))));
GFP = upper(cell2mat(ref_seqs(2, contains(ref_seqs(1, :), 'GFP'))));
A4 = upper(cell2mat(ref_seqs(2, contains(ref_seqs(1, :), 'A4'))));
A4 = A4(end-50:end);
minP = upper(cell2mat(ref_seqs(2, contains(ref_seqs(1, :), 'minP'))));
mRuby = upper(cell2mat(ref_seqs(2, contains(ref_seqs(1, :), 'mRuby'))));
mRuby = mRuby(200:300);
BFP = upper(cell2mat(ref_seqs(2, contains(ref_seqs(1, :), 'BFP'))));
BFP = BFP(200:300);
BS10_1 = upper(cell2mat(ref_seqs(2, contains(ref_seqs(1, :), 'BS10_1'))));

%Load variable sequences
ZF_parts = ref_seqs(2, 7:9)';
minP_100k = readcell('../example_ref_sequences/100k_Minimal-Promoters.xlsx'); minP_100k = minP_100k(2:end, 2);
terminators_100k = readcell('../example_ref_sequences/100k_Terminators.xlsx'); terminators_100k = terminators_100k(2:end, 2);
spacers_100k = readcell('../example_ref_sequences/100k_Terminator_spacers.xlsx'); spacers_100k = spacers_100k(2:end, 2);
promoters_100k = readcell('../example_ref_sequences/100k-Promoters.xlsx'); promoters_100k = promoters_100k(2:end, 2);
ORF_parts = readcell('../example_ref_sequences/100k_ORF-Parts_SynTF.xlsx'); 
AD_parts = ORF_parts(2:5, 3); IDR_parts = ORF_parts(6:9, 3);

%% Consolidate reads using fastqall 
[~, l, seq] = fastqall('../example_fastq', 'fastq');
%Filter based on read length
seq = seq(l > 9500 & l < 15000); l = l(l > 9500 & l < 15000);

%% Index to Puromycin using bowtile
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
[minP_variants_scaled, ~, minPconf] = viscount(pregions, 6, minP_100k, 0.2, 'T');
[~,minP_variants] = max(minP_variants_scaled, [], 2);
minP_variants((sum(minP_variants_scaled(:,1:3),2) < 0.2)) = 0;

%Number of Binding sites using fastar
[variants_bs, ~] = fastar(pregions, BS10_1, 6, 8);

%Everything between
variants_bs(variants_bs > 6.9 & variants_bs < 9.1) = 8;
variants_bs(variants_bs > 9.2 & variants_bs < 14) = 12;
variants_bs(variants_bs > 3.9 & variants_bs < 6.1) = 4;
variants_bs(variants_bs ~= 2 & variants_bs ~= 4 & variants_bs ~= 8 & variants_bs ~= 12) = 0;

variants_bs(variants_bs == 2) = 1; %2 binding sites
variants_bs(variants_bs == 4) = 2; %4 binding sites
variants_bs(variants_bs == 8) = 3; %8 binding sites
variants_bs(variants_bs == 12) = 4; %12 binding sites

%Terminator spacing
[variants_term, ~, ~] = viscount(tregions, 10, spacers_100k(2), 0.2, 'F');
term_indices = ones(size(variants_term));
term_indices(variants_term > 0.2 & variants_term < 0.52) = 2;
term_indices(variants_term > 0.52) = 3;

%Calculate assignments: minP_variants_scaled, variants_BS, term_indices
reporter_variants = [variants_bs minP_variants term_indices];

%% SynTF Identification
%Create pregions and tregions using chophat
pregions_synTF = chophat(reads_correct, [ones(size(positions2, 1), 1) positions2(:, 4)], 0, 0);
tregions_synTF = chophat(reads_correct, [positions2(:, 4) positions2(:, 1)], 0, 0);

%Assign Promoters
[synTFprom_variants_scaled, ~, synTFprom_conf] = viscount(pregions_synTF,10,promoters_100k, 0.03, 'T');
[~,variants_synTF_prom] = max(synTFprom_variants_scaled, [], 2); variants_synTF_prom = variants_synTF_prom';
variants_synTF_prom(sum(synTFprom_variants_scaled,2) < 0.03) = 0;

%Assign ADs
[~, f, AD_conf] = viscount(pregions_synTF, 10, AD_parts, 0.2, 'T');
f(f < 30) = 0; 
[~,variants_AD] = max(f, [], 2); variants_AD(sum(f,2) == 0) = -1;

AD_conf(4, 4) = sum(f(:, 4) > 10 & f(:, 3) < 100);
AD_conf(3, 3) = sum(f(:, 3) > 100);
AD_conf(3, 4) = sum(f(:, 4) < 10 & f(:, 3) > 100);
AD_conf(4, 3) = AD_conf(3, 4);

%Assign IDRs
[IDR_tiles, ~, IDR_conf] = viscount(pregions_synTF, 10, IDR_parts, 0.15, 'T');
[~, variants_IDR] = max(IDR_tiles, [], 2);

%Assign Zinc Fingers
[ZF_tiles, ~, ZF_conf] = viscount(tregions_synTF, 10, ZF_parts, 0.03, 'T');
[~,variants_ZF] = max(ZF_tiles, [], 2);
variants_ZF(sum(ZF_tiles,2) == 0) = -1;

%Assign terminators
[synTF_term, ~, term_conf] = viscount(tregions_synTF, 10, terminators_100k, 0.4, 'T');
[~,variants_synTF_term] = max(synTF_term, [], 2);

%assign spacers
[spacer_synTF, ~, ~] = viscount(tregions_synTF, 10, spacers_100k(2), 0.4, 'F');
variants_synTF_spacer = ones(size(spacer_synTF));
variants_synTF_spacer(spacer_synTF > 0.2 & spacer_synTF < 0.52) = 2;
variants_synTF_spacer(spacer_synTF > 0.52) = 3;

%Calculate assignments
synTF_variants = [variants_synTF_prom variants_AD variants_IDR variants_ZF variants_synTF_term variants_synTF_spacer];
%% Assign orientation
%based on location of GFP and mRuby
variants_orientation = ones(size(reads_correct, 1), 1);
variants_orientation(positions2(:, 2) < positions2(:, 4)) = 2;
%% Sum all synTF and reporter indices to create full list of 100k identified circuits
all_100k_variants = [synTF_variants reporter_variants variants_orientation];
%remove all rows containing zero or -1
all_100k_variants(any(all_100k_variants < 1,2),:) = [];
%library_100k = all_100k_variants(:,1) + 4*(all_100k_variants(:,2)-1) + 16*(all_100k_variants(:,3)-1) + 80*(all_100k_variants(:,4)-1) + 240*(all_100k_variants(:,5)-1) + 2880*(all_100k_variants(:,6)-1) + 11520*(all_100k_variants(:,7)-1) +  34560*(all_100k_variants(:,8)-1);
library_100k = sum((all_100k_variants - [0 ones(1, 9)]).*[1 4 16 64 192 768 2304 9216 27648 82944], 2);
library_100k(library_100k < 1) = 0;
%% Barcode Assignments
tregion_BFP = chophat(reads_correct, positions2(:, 5), 1000, 1);
%Identify bc1 and 2
[bc1, bc1_l, bc1_score, bc1_pos] = barcoat(tregion_BFP, 'BBA', 0);
[bc2, bc2_l, bc2_score, bc2_pos] = barcoat(tregion_BFP, 'DDC', 0);
%% Data QC and filtering
is_valid_assignment = bc1_score > 85 & bc2_score > 90;
is_fully_assigned = library_100k > 0;

bc1_correct = bc1(is_valid_assignment & is_fully_assigned);
bc2_correct = bc2(is_valid_assignment & is_fully_assigned);
barcoded_variants = all_100k_variants(is_valid_assignment & is_fully_assigned,:);
library_variant_indices = library_100k(is_valid_assignment & is_fully_assigned);
