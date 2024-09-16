function [corrseq, bc_ident]  = bc_detect(seq, lengths, bc_sequence, thresh)

a = callbarcode;

if isnumeric(bc_sequence)
    bin_bc = a(2*bc_sequence - 1, :);
    bin_bc(2, :) = a(2*bc_sequence, :);
else
    bin_bc = bc_sequence;
    bin_bc(2, :) = seqrcomplement(bc_sequence);
end


corrseq = cell(size(seq));
bc_start = zeros(size(seq));
bc_ident = bc_start;

if isempty(lengths)
    for i = 1:length(seq)
            %disp(i)
            a = cell2mat(seq(i));
            lengths(i) = length(a);
    end
end


for i = 1:length(seq)
    if mod(i, 1000) == 0
        disp(i)
    end
    a = cell2mat(seq(i));
    if lengths(i) > 150
        b = a(1:150);
        [~, y, z] = swalign(b, bin_bc(1, :), 'Alphabet', 'NT');
        [~, y1, z1] = swalign(b, bin_bc(2, :), 'Alphabet', 'NT');
        s1 = sum(y(2, :) == '|')/size(bin_bc, 2);
        s2 = sum(y1(2, :) == '|')/size(bin_bc, 2);
        if s1 >= thresh
            corrseq(i) = cellstr(a);
            bc_start(i) = z(1);
            bc_ident(i) = 1;
        elseif s2 >= thresh
                corrseq(i) = cellstr(seqrcomplement(a));
                bc_start(i) = z1(1);
                bc_ident(i) = 2;
        else
            corrseq(i) = cellstr('No BC Found');
            bc_ident(i) = 0;
            bc_start(i) = 0;
        end
    end
end



end