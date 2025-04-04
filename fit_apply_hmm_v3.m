function tbl = fit_apply_hmm_v3(tbl)
% ------------------------- Fit HMM -----------------------------------
% This function takes the morphological label sequences in the properties structure
% 'props' of each row in a stacktable and uses them to fit the parameters 
% of an HMM model using the Baum-Welch algorithm. The fitted state
% transition matrix A_est and initial state distribution pi_est is then
% used in the Viterbi algorithm to find the most likely sequence of hidden 
% states given an input sequence of observations (morphological labels).
% For each combination of celline, perturbation and dose, all sequences are
% aggregated and used to fit a model, which is then used to correct the
% same set of sequences.
%
% Input parameters:
% tbl - stacktable of ROIs
%
% Output parameters:
% tbl - stacktable with added variables;
%            A_est - fitted state transition matrix 
%            pi_est - fitted initial state distribution
%            Viterbi_paths - corrected observation sequence (cell array)
%            Viterbi_paths2 - corrected observation sequence (matrix;
%                       row=frame, col=cell, entries=morph. label)
%
% @authors: Madeleine Skeppås, Sven Nelander
% @date: 15012025

% Add empty columns to table
tbl.A_est = cell(height(tbl),1);
tbl.pi_est = cell(height(tbl),1);

N = 6; % Number of states
M = N; % Number of observation symbols 
p = 100; % Number of sequences
AI_quality = 0.75; % Fraction correct classes by the observation AI

% Parameters - model fitting 
minSeqLen = 5;
maxSeqLen = 20;
maxIter = 100;
tol = 1e-4;
mu = 2; % Strength of the diagonal prior
lambda = 2; % L1-regularization parameter for sparsity

B= [[28 4 4 7 0 0] % 6 classes
    [1 127 11 26 2 0]
    [1 0 60 0 1 0]
    [2 5 13 116 0 1]
    [0 7 4 10 152 2]
    [0 0 1 0 0 196]];

B = B ./ sum(B, 2); % Normalize

pi_init = ones(1, N) / N; % Uniform initial state distribution

cellines = unique(tbl.HGCC);

% Fit A_est and pi_est for every comb of celline and perturbation

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
            idx = logical((tbl.HGCC == string(hgcc)) .* (tbl.perturbation == string(pert)) .* (tbl.dose == dose_curr));

            for l=1:height(subtable_3)
                seqs = subtable_3.props{l};
                seqs = seqs{end-1};

                for n=1:width(seqs)
                    seq = seqs(:,n);
                    start = min(find(~isnan(seq)));
                    stop = max(find(~isnan(seq)));

                    sequences = [sequences; seq(start:stop)'];
                end
            end

            A_est_avg = zeros(6,6,1000);
            pi_est_avg = zeros(1,6,1000);

            % Run model fitting with 1000 random versions of A_init
            parfor m=1:1000
                % Set seed and generator
                rng(m,'twister')

                % Create random initial transition matrix
                A_init = rand(N)+5*eye(N);
                A_init = A_init ./ sum(A_init, 2); % Normalize

                % Run model fitting
                [A_est, pi_est] = fit_transition_model(sequences, A_init, B, pi_init, maxIter, tol, mu, lambda);
                
                A_est_avg(:,:,m) = A_est;
                pi_est_avg(:,:,m) = pi_est;
            end

            tbl.A_est(idx) = {mean(A_est_avg,3)};
            tbl.pi_est(idx) = {mean(pi_est_avg,3)};
        end
    end
    fprintf(['Fit HMM parameters for celline: ' hgcc '...\n'])
end

% ------------------------------ Clean tracks ----------------------------

% Add empty column to table
tbl.Viterbi_paths = cell(height(tbl),1);
tbl.Viterbi_paths2 = cell(height(tbl),1);

% Apply transition model to clean sequences
for i=1:height(tbl)
    sequences = {};
    
    stack = tbl(i,:);
    seqs = stack.props{:};
    seqs = seqs{end-1};

    for l=1:width(seqs)
        seq = seqs(:,l);
        start = min(find(~isnan(seq)));
        stop = max(find(~isnan(seq)));
    
        sequences = [sequences; seq(start:stop)'];
    end

    A = stack.A_est{:};
    pi = stack.pi_est{:};

    viterbi_paths = apply_transition_model(sequences, A, B, pi);
    
    viterbi_paths2 = apply_transition_model2(seqs, A, B, pi);

    tbl.Viterbi_paths(i) = {viterbi_paths};

    tbl.Viterbi_paths2(i) = {viterbi_paths2};
end

end