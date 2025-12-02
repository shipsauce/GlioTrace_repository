function [sequences, hard_labels, propagated_labels, tme_features, delta_ts, comb_features] = convert_data_for_hmm(tbl)

cellines = unique(tbl.HGCC);
% Loop through the cellines
for i=1:length(cellines)
    hgcc = cellines{i};
    subtable = tbl(tbl.HGCC == string(hgcc),:);

    perturbations = unique(subtable.perturbation);

    % Loop through the perturbations
    for j=1:length(perturbations)
        pert = perturbations{j};
        subtable_2 = subtable(subtable.perturbation == string(pert),:);
        dosez = unique(subtable_2.dose);

        % Loop through the doses
        for k=1:length(dosez)
            dose_curr = dosez(k);
            subtable_3 = subtable_2(subtable_2.dose == dose_curr,:);
            sequences = {};
            hard_labels = {};
            features = {};
            idx = logical((tbl.HGCC == string(hgcc)) .* (tbl.perturbation == string(pert)) .* (tbl.dose == dose_curr));
            rowcount = 1;

            % Aggregate observations across ROIs
            for l=1:height(subtable_3)
                info = subtable_3.props{l};
                embds = mapCoords_embed(subtable_3.traxs{l}, subtable_3.trays{l}, ...
                    subtable_3.x_coords{l}, subtable_3.y_coords{l}, subtable_3.embeddings_long{l}, subtable_3.startidx(l));
                tme_info = info{end};
                hard_labs = weight_hard_labels(info, subtable_3.HGCC{1});
                
                labs = subtable_3.propagated_labels{l};
                nan_idx = sum(isnan(labs),1);
                labs = labs(:,nan_idx == 0);

                [embds, hard_labs, x_coords, y_coords, tme_info, deltat] = handle_missing_emissions_v2(embds, hard_labs, tme_info, subtable_3.traxs{l}, subtable_3.trays{l}, subtable_3.delta_t(l), labs);

                radius = 50;
                [localDensity, microgliaFraction, polarization, vesselCloseness] = construct_features(x_coords, y_coords, tme_info, subtable_3.vesselMasks{l}, deltat, radius);

                % Loop through observation sequences (cells)
                for n=1:width(tme_info)
                    embd = embds(:,n);
                    hard_lab = hard_labs(:,n);
                    f1 = localDensity(:,n);
                    f2 = microgliaFraction(:,n);
                    f3 = polarization(:,n);
                    f4 = vesselCloseness(:,n);
                    f5 = deltat(:,n);

                    deltat = repmat(subtable_3.delta_t(l), [1 sum(~isnan(hard_lab))-1])';
                    prop_lab = labs(~isnan(labs(:,n)),n);

                    start = min(find(~isnan(tme_assoc)));
                    stop = max(find(~isnan(tme_assoc)));
    
                    sequences{rowcount} = cell2mat(embd(start:stop));
                    hard_labels{rowcount} = hard_lab(start:stop);
                    
                    tme_features{rowcount} = tme_assoc(start:stop);
                    delta_ts{rowcount} = deltat(1:end-1);
                    % try
                    % comb_features{rowcount} = [[deltat(1:end-1); nan(1,1)] tme_assoc(start:stop)];
                    % catch
                    %     1
                    % end
                    rowcount = rowcount+1;

                end
            end
        end
    end
end

sequences = sequences';
hard_labels = hard_labels';

tme_features = tme_features';
delta_ts = delta_ts';
% comb_features = comb_features';
end