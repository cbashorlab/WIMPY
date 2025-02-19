function [nums, locs] = FASTar(pregions, ref, step, bw)

nums = zeros(size(pregions));
locs = cell(size(pregions));

for j = 1:length(pregions)
    x = cell2mat(pregions(j));
    a1 = [];
    if ~contains(x, 'X')
    for i = 1:length(ref) - step
        a = ref(i:i+step);
        y = strfind(x, a);
        if ~isempty(y)
            a1(end+1:end+length(y)) = y;
        end
    end
    if ~isempty(a1)
        [f, x] = ksdensity(unique(a1), 'bandwidth', bw);
        a1 = x(islocalmax(f));
        locs(j) = {a1};
        nums(j) = length(a1);
    end
    end
end
