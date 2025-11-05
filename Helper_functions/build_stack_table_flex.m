function t = build_stack_table_flex(stackfile,metadata)

% If stackfile is .txt file
try
    t=readtable(stackfile,'ReadVariableNames',false,'delimiter','\t');
    t=t.Var1;
catch
    % If stackfile is cell array
    t=table(stackfile);
    t.Properties.VariableNames = {'file_path'};
end

experiment_id = cellfun(@(x) extract_id(x), t.file_path, 'UniformOutput',false);
t.experiment_id = experiment_id;

if(nargin > 1)
    lookup = zeros(size(t.experiment_id));  % preallocate
    
    for i = 1:numel(t.experiment_id)
        idx = find(contains(metadata.filename_experiment_id, t.experiment_id{i}), 1);
        if ~isempty(idx)
            lookup(i) = idx;
        else
            lookup(i) = 0;  % no match
        end
    end

    t = [t metadata(lookup,2:end)];
end


end

function id = extract_id(str)
match = regexp(str, '/([^/]+)_roi', 'tokens');
id = match{1}{1};
end

