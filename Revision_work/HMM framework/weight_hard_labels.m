function hard_labs = weight_hard_labels(props, celline)
% Load class weights
t = readtable('../brainslice_manuscript_repo/Revision_work/class_proportions.csv');
idx = contains(t.Var1, celline);

weights = table2array(t(idx,2:end));

% Take soft labels, weigh by class prevalence

[frames, cells] = size(props{1});

hard_labs = nan(size(props{1}));

for i=1:frames
    for j=1:cells
        soft_labs = cell2mat(cellfun(@(x) x(i,j), props(1:6), 'UniformOutput',false));
        soft_labs_weighted = soft_labs .* weights;
        [~,loc] = max(soft_labs_weighted);
        if(all(~isnan(soft_labs_weighted)))
            hard_labs(i,j) = loc;
        end
    end
end
end