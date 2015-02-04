function [ r, p ] = main_correlation( dataset, difficulty, diff_rate, measure, parameter, permute )
%The all-in function to compute the correlation of a specific measure
%MEASURE with learning rates on the first day for all subjects (good
%subjects if the DATASET is the 3day dataset), and of a specific DIFFICULTY
%if the DATASET is the 6week dataset. The potential PARAMETERs besides
%DIFFUSION_RATE are defined in structural format in PARAMETER.
%
%  INPUT: dataset    - '3day' (1) or '6week' (2)
%         difficulty - 1 for ext or 2 for mod or 3 for min for '6week'
%                      dataset, or vector with values 1 - 3
%         diff_rate  - diffusion rate, 0 if not diffused
%         measure    - a scalar integer 1 - 7 corresponding to
%                         1. number of switches
%                         2. number of switches greater than a threshold
%                         3. number of brain signal differences > threshold
%                         4. weighted total switches after clustering
%                         5. weighted total brain signal differences
%                         6. approximate entropy of the clustering results
%                         7. sample entropy of the clustering results
%         parameter  - potential parameters
%                         .temporal - temporal normalizing factor,
%                                     used in measures 1,2,4,6,7
%                         .numClust - number of clusters,
%                                     used in measures 1,2,4,6,7
%                         .alpha    - threshold in measures 2,3
%                                     weight in measures 4,5
%                         .r        - tolerance in entropy measures 6,7
%                                     (typically 0.2 * std)
%         permute    - flag for permutation test (either subject are 
%                      permutated or for each subject, her activity 
%                      profiles are permutated)
%                         0 or undefined if not permuted
%                         .learningOrder  - 1 if subjects' learning rates
%                                           permuted 
%                         .timeOrder      - rng seed for time order
% OUTPUT: r          - correlation between measures and learning 
%         p          - p-values between measures and learning
%
%---------------------------------------------------------------------
% coded by Weiyu Huang,  whuang@seas.upenn.edu
% Ver 1 : Jan 14, 2014
%---------------------------------------------------------------------


%--------------------------------------------------------------------%
%------------------ 1. check necessary parameters -------------------%
%--------------------------------------------------------------------%

% Note: 1. if unnecessary parameters exist, set them to -1
%       2. temporal 0 meaning default
%       3. if r is not defined, 0.2 * std is the most commonly used one
%
if ~exist('permute', 'var')
    permute.learningOrder = 0;
end

switch measure
    case 1
        parameter.alpha = -1; parameter.r = -1;
        if(~isfield(parameter, 'temporal')); parameter.temporal = 0; end;
        if(~isfield(parameter, 'numClust')); parameter.numClust = 10; end;
    case {2, 4}
        parameter.r = -1;
        if(~isfield(parameter, 'temporal')); parameter.temporal = 0; end;
        if(~isfield(parameter, 'numClust')); parameter.numClust = 10; end;
        if(~isfield(parameter, 'alpha')); parameter.alpha = 1; end
    case {3, 5}
        parameter.numClust = -1; parameter.temporal = -1; parameter.r = -1;
        if(~isfield(parameter, 'alpha')); parameter.alpha = 1; end;
    case {6, 7}
        parameter.alpha = -1; 
        if(~isfield(parameter, 'temporal')); parameter.temporal = 0; end;
        if(~isfield(parameter, 'numClust')); parameter.numClust = 10; end;
        if(~isfield(parameter, 'r')); parameter.r = 0.2; end
end

%--------------------------------------------------------------------%
%--------------------- 2. compute correlations ----------------------%
%--------------------------------------------------------------------%

% Compute the number of parameter combinations we want to have (number of
% rows P in the result table) and number of columns D
%
N_temporal = numel(parameter.temporal);
N_numClust = numel(parameter.numClust);
N_alpha = numel(parameter.alpha);
N_r = numel(parameter.r);
P = N_temporal * N_numClust * N_alpha * N_r;
D = numel(difficulty);
r = zeros(P, D); p = zeros(P, D);

% Loop over the table and fill in one entry at each time by calling the
% subroutine
%
for P_idx = 1:P
    for D_idx = 1:D
        % P_idx are ordered first by r, then by alpha, then numClust,
        % finally temporal
        %
        crt.r = parameter.r( mod(P_idx - 1, N_r) + 1 );
        crt.alpha = parameter.alpha( floor( mod(P_idx - 1, N_r * N_alpha) / N_r ) + 1 );
        crt.numClust = parameter.numClust( floor( mod(P_idx - 1, ...
            N_r * N_alpha * N_numClust) / (N_r * N_alpha) ) + 1 );
        crt.temporal = parameter.temporal( floor( (P_idx - 1) / ...
            (N_r * N_alpha * N_numClust) ) + 1 );
        switch dataset
            case {'3day', 1}
                [r(P_idx, D_idx), p(P_idx, D_idx)] = ...
                    compute_correlation_3day(diff_rate, measure, crt, permute.learningOrder);
            case {'6week', 2}
                [r(P_idx, D_idx), p(P_idx, D_idx)] = ...
                    compute_correlation_6week(difficulty(D_idx), diff_rate, measure, crt, permute.learningOrder);
        end
    end
end

end

