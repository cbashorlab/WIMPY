function out = chophat(reads_correct, positions, truncationcriteria, retain)

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

