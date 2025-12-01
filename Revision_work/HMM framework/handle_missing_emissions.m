function [embd, hard_lab, tme_assoc, deltat_new, prop_lab] = handle_missing_emissions(embd, hard_lab, tme_assoc, deltat, prop_lab)
startidx = min(find(~isnan(hard_lab)));
embd = embd(prop_lab ~= 1, :);
hard_lab = hard_lab(prop_lab ~= 1, :);
tme_assoc = tme_assoc(prop_lab ~= 1, :);

keep = ~prop_lab;
keep = keep(startidx:end);
dt_new = [];

for i = 1:length(deltat)
    if keep(i)
        last_valid = i;
        dt_new(end+1) = deltat(i);
    else
        dt_new(last_valid) = dt_new(last_valid) + deltat(i);
        dt_new(i) = 0;
    end
    i
end

deltat_new = dt_new(~dt_new==0)';

end