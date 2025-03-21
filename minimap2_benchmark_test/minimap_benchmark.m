
% Number of binding sites

data = readcell('alignment_numBS.csv');

read_name = data(2:end, 1);
percalign = cell2mat(data(2:end, 13));
nbs = cell2mat(data(2:end, 6));

[un, ~, b] = unique(read_name);
un_counts = accumarray(b, 1);

nbs_assignments = zeros(size(un));
nbs_max = nbs_assignments;
nbs_conf = zeros(4, 4);

for i = 1:length(un)
    x = un(i);
    a = [nbs(contains(read_name, x)) percalign(contains(read_name, x))];
    amax = [max([a(a(:, 1) == 2, 2); 0]) max([a(a(:, 1) == 4, 2); 0]) max([a(a(:, 1) == 8, 2); 0]) max([a(a(:, 1) == 12, 2); 0])];
    y = find(amax > 0.64);
    if isscalar(y)
        nbs_conf(y, y) = nbs_conf(y, y) + 1;
    elseif length(y) == 2
        nbs_conf(y(1), y(2)) = nbs_conf(y(1), y(2)) + 1;
        nbs_conf(y(2), y(1)) = nbs_conf(y(2), y(1)) + 1;
    end
        
    nbs_assignments(i) = find(amax == max(amax));
    nbs_max(i) = max(amax);

end

minimapout_bs = table(un, nbs_assignments, nbs_max, 'VariableNames', {'Read Name', 'Assignment', 'Perc. Align'});
writetable(minimapout_bs, 'Minimap_bs_assignments.csv');


%% Minimap synTF Promoters

data = readcell('alignment_prom.csv');

read_name = data(2:end, 1);
percalign = cell2mat(data(2:end, 13));
prom = data(2:end, 6);

[un, ~, b] = unique(read_name);
un_counts = accumarray(b, 1);

prom_assignments = zeros(size(un));
prom_max = nbs_assignments;
prom_conf = zeros(4, 4);

for i = 1:length(un)
    x = un(i);
    a = prom(contains(read_name, x)); ap = percalign(contains(read_name, x));
    amax = [max([ap(contains(a, 'hPGK')); 0]) max([ap(contains(a, 'CMV')); 0]) max([ap(contains(a, 'RSV')); 0]) max([ap(contains(a, 'hEf1a1')); 0])];
    y = find(amax > 0.64);
    if isscalar(y)
        prom_conf(y, y) = prom_conf(y, y) + 1;
    elseif length(y) == 2
        prom_conf(y(1), y(2)) = prom_conf(y(1), y(2)) + 1;
        prom_conf(y(2), y(1)) = prom_conf(y(2), y(1)) + 1;
    end
    
    if sum(amax(1:3) > 0.65) > 0
        prom_assignments(i) = find(amax(1:3) == max(amax(1:3)));
        prom_max(i) = max(amax(1:3));
    elseif amax(4) > 0.65
        prom_assignments(i) = 4;
        prom_max(i) = amax(4);
    end

end

for i = 1:3
    prom_conf(i, i) = prom_conf(i, i) + prom_conf(i, 4);
    prom_conf(i, 4) = 0;
    prom_conf(4, i) = 0;
end

minimapout_prom = table(un, prom_assignments, prom_max, 'VariableNames', {'Read Name', 'Assignment', 'Perc. Align'});
writetable(minimapout_prom, 'Minimap_prom_assignments.csv');

%% WIMPY SynTF Promoters

%Promoters
addpath('../wimpy_matlab/wimpy_helper_functions/')

Puro = upper('cgctccgcatcggcctaaggaaccggcgtggttcctggctacggtgggagtctcacctgaccatcaaggaaagggattgggaagtgctgtcgttcttcca');
mRuby = upper('GTGAGTAAAGGCGAAGAACTTATCAAGGAAAATATGCGGATGAAAGTGGTTATGGAGGGTAGCGTGAACGGACACCAGTTCAAATGCACGGGAGAGGGCGAGGGGCGACCCTACGAGGGAGTCCAAACAATGAGGATTAAGGTTATAGAAGGTGGTCCGCTGCCATTCGCATTCGATATTTTGGCCACGTCCTTCATGTACGGCTCCCGAACCTTTATCAAATACCCTGCGGATATCCCAGACTTTTTCAAGCAATCCTTTCCGGAAGGGTTCACGTGGGAGCGAGTCACGAGATATGAGGATGGAGGCGTAGTAACAGTAACCCAAGACACATCACTTGAGGACGGTGAGCTTGTCTACAATGTGAAGGTACGCGGCGTCAATTTCCCCTCAAATGGCCCGGTGATGCAAAAGAAAACTAAAGGATGGGAGCCCAACACCGAAATGATGTACCCGGCAGATGGGGGGCTTAGGGGCTATACGGACATCGCATTGAAGGTTGATGGCGGGGGCCATCTCCATTGTAACTTTGTAACTACATATCGGTCAAAAAAGACTGTGGGGAACATTAAAATGCCGGGAGTACACGCTGTTGATCATCGCCTGGAAAGGATAGAGGAAAGCGACAATGAAACGTATGTAGTACAGCGGGAGGTCGCCGTCGCCAAATATAGTAATCTGGGCGGTGGCATGGACGAGCTTTATAAA');
mRuby = mRuby(200:300);
promoters_100k = readcell('./info/100k-Promoters.xlsx'); promoters_100k = promoters_100k(2:end, 2);
promoters_100k = promoters_100k([2, 3, 4, 1]);

[hs, l, seq] = fastqall('./fastq', 'fastq_runid');

thresh = 0.03;
[new_seq, ~, ~, ~, ~, ~] = bowtile(seq,Puro,thresh);
reads_correct = new_seq(~contains(new_seq, 'X')); l_readscorrect = l(~contains(new_seq, 'X'));
hs_correct = hs(~contains(new_seq, 'X'));

[~, positionsmRuby, ~] = tilepin(reads_correct, mRuby(end-100:end), thresh, 'F');
positionsmRuby = floor(positionsmRuby);

pregions_synTF = chophat(reads_correct, [ones(size(positionsmRuby, 1), 1) positionsmRuby], 0, 0);
[synTFprom_variants_scaled, ~, synTFprom_conf] = viscount(pregions_synTF,10,promoters_100k, 0.03, 'T');
synTFprom_variants_scaled = [zeros(size(synTFprom_variants_scaled, 1), 1) synTFprom_variants_scaled];

[vp,variants_synTF_prom] = max(synTFprom_variants_scaled, [], 2);
variants_synTF_prom = variants_synTF_prom - 1;

wimpyout_prom = table(hs_correct, variants_synTF_prom, vp, 'VariableNames', {'Read Name', 'Assignment', 'Perc. Align'});
writetable(wimpyout_prom, 'Wimpy_prom_assignments.csv');

%% WIMPY Binding sites

GFP = upper('atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa');
A4 = upper('ggagcatcgcccttccccggccctcaggtaagaggaccaaataccgtagccgtttccaatttcagtcctttagcgccacctggtgctaactactctatcacgcttttatccaataactacctttgtaaatttcctttcaaaagttctggccgggcgcggtggcTGTAC');
A4 = A4(end-50:end);

[~, positionsA4, ~] = tilepin(reads_correct, A4, thresh, 'F');
[~, positionsGFP, ~] = tilepin(reads_correct, GFP(end-100:end), thresh, 'F');
positions2 = floor([positionsA4, positionsGFP]);

pregions = chophat(reads_correct, positions2(:, 1:2), 0, 0);
BS10_1 = upper('cGGCGTAGCCGATGTCGCGc');
[bs_num, ~] = FASTar(pregions, BS10_1, 6, 8);

bs_num(bs_num > 6.9 & bs_num < 9.1) = 8;
bs_num(bs_num > 9.2 & bs_num < 14) = 12;
bs_num(bs_num > 3.9 & bs_num < 6.1) = 4;

bs_num(bs_num ~= 2 & bs_num ~= 4 & bs_num ~= 8 & bs_num ~= 12) = 0;

bs_num_wimpy = zeros(size(bs_num));
bs_num_wimpy(bs_num == 2) = 1;
bs_num_wimpy(bs_num == 4) = 2;
bs_num_wimpy(bs_num == 8) = 3;
bs_num_wimpy(bs_num == 12) = 4;

wimpyout_bs = table(hs_correct, bs_num_wimpy, ones(size(bs_num_wimpy)), 'VariableNames', {'Read Name', 'Assignment', 'Perc. Align'});
writetable(wimpyout_bs, 'Wimpy_bs_assignments.csv');

%% Ground truths - num_bs

minP = upper('gctcactctcccttacacggagtggataTAGT');
[~, positionsminP, ~] = tilepin(reads_correct, minP, thresh, 'F');

positions2(:, 3) = floor(positionsminP);

variants_bs_gt = zeros(length(reads_correct), 1); a = positions2(:, 3) - positions2(:, 1);
variants_bs_gt(a > 60 & a < 160) = 1; %2 binding sites
variants_bs_gt(a > 160 & a < 255) = 2; %4 binding sites
variants_bs_gt(a > 255 & a < 415) = 3; %8 binding sites
variants_bs_gt(a > 415 & a < 555) = 4; %12 binding sites

%% Ground truths - SynTF promoters

synTFprom_variants_scaled = [zeros(size(synTFprom_variants_scaled, 1), 1) synTFprom_variants_scaled];
[vp,variants_synTF_prom] = max(synTFprom_variants_scaled, [], 2);

variants_synTF_gt = variants_synTF_prom -1;
variants_synTF_gt(vp < 0.5) = 0;

%% Compare results

%Num_bs first

x = readcell("Wimpy_bs_assignments.csv");
reads_wimpy = x(2:end, 1);
numbs_wimpy = cell2mat(x(2:end, 2));

for i = 1:length(reads_wimpy)
    a = cell2mat(reads_wimpy(i));
    b = strfind(a, 'runid');
    reads_wimpy(i) = cellstr(a(1:b-2));
end


x = readcell("Minimap_bs_assignments.csv");
reads_minimap = x(2:end, 1);
numbs_minimap = cell2mat(x(2:end, 2));

reads_unique = unique([reads_wimpy; reads_minimap]);

wimpy_stats = zeros(1, 3); minimap_stats = wimpy_stats; gt_corr = 0;

for i = 1:length(reads_unique)
    a = numbs_wimpy(contains(reads_wimpy, cell2mat(reads_unique(i))));
    if isempty(a); a = 0; end
    b = numbs_minimap(contains(reads_minimap, cell2mat(reads_unique(i))));
    if isempty(b); b = 0; end
    c = variants_bs_gt(contains(hs_correct, cell2mat(reads_unique(i))));
    if isempty(c); c = 0; end

    if c ~= 0
        gt_corr = gt_corr + 1;
        if a == c
            wimpy_stats(1) = wimpy_stats(1) + 1;
        elseif a == 0
            wimpy_stats(2) = wimpy_stats(2) + 1;
        else
            wimpy_stats(3) = wimpy_stats(3) + 1;
        end
        if b == c
            minimap_stats(1) = minimap_stats(1) + 1;
        elseif b == 0
            minimap_stats(2) = minimap_stats(2) + 1;
        else
            minimap_stats(3) = minimap_stats(3) + 1;
        end
    end
end

mm_perc = 100*minimap_stats/gt_corr;
wp_perc = 100*wimpy_stats/gt_corr;


% Proms next

x = readcell("Wimpy_prom_assignments.csv");
reads_wimpy = x(2:end, 1);
proms_wimpy = cell2mat(x(2:end, 2));

for i = 1:length(reads_wimpy)
    a = cell2mat(reads_wimpy(i));
    b = strfind(a, 'runid');
    reads_wimpy(i) = cellstr(a(1:b-2));
end

x = readcell("Minimap_prom_assignments.csv");
reads_minimap = x(2:end, 1);
proms_minimap = cell2mat(x(2:end, 2));

reads_unique = unique([reads_wimpy; reads_minimap]);
wimpy_stats_prom = zeros(1, 3); minimap_stats_prom = wimpy_stats_prom; gt_corr_prom = 0;

for i = 1:length(reads_unique)
    a = proms_wimpy(contains(reads_wimpy, cell2mat(reads_unique(i))));
    if isempty(a); a = 0; end
    b = proms_minimap(contains(reads_minimap, cell2mat(reads_unique(i))));
    if isempty(b); b = 0; end
    c = variants_synTF_gt(contains(hs_correct, cell2mat(reads_unique(i))));
    if isempty(c); c = 0; end

    if c ~= 0
        gt_corr_prom = gt_corr_prom + 1;
        if a == c
            wimpy_stats_prom(1) = wimpy_stats_prom(1) + 1;
        elseif a == 0
            wimpy_stats_prom(2) = wimpy_stats_prom(2) + 1;
        else
            wimpy_stats_prom(3) = wimpy_stats_prom(3) + 1;
        end
        if b == c
            minimap_stats_prom(1) = minimap_stats_prom(1) + 1;
        elseif b == 0
            minimap_stats_prom(2) = minimap_stats_prom(2) + 1;
        else
            minimap_stats_prom(3) = minimap_stats_prom(3) + 1;
        end
    end
end

mm_perc_prom = 100*minimap_stats_prom/gt_corr_prom;
wp_perc_prom = 100*wimpy_stats_prom/gt_corr_prom;

