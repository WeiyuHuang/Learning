function [ r, p ] = main_correlation( dataset, diff_rate, measure, parameter, permute )
%The all-in function to compute the correlation of a specific measure
%MEASURE with learning rates on the first day for all subjects (good
%subjects if the DATASET is the 3day dataset), and of a specific DIFFICULTY
%if the DATASET is the 6week dataset. The potential PARAMETERs besides
%DIFFUSION_RATE are defined in structural format in PARAMETER.
%
%  INPUT: dataset    - '3day' (1) or '6week' (2)
%         diff_rate  - diffusion rate, 0 if not diffused, may be a vector
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
% Ver 1  : Jan 14, 2014
% Ver 1.1: Feb 4, 2014
%---------------------------------------------------------------------


%--------------------------------------------------------------------%
%------------------ 1. check necessary parameters -------------------%
%--------------------------------------------------------------------%

% Note: 1. if unnecessary parameters exist, set them to -1
%       2. temporal 0 meaning default
%       3. if r is not defined, 0.2 * std is the most commonly used one
%
if ~exist('permute', 'var')
    permute = 0;
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

% Call subroutines.
N_diffrate = numel(diff_rate);
N_temporal = numel(parameter.temporal);
N_numClust = numel(parameter.numClust);
N_alpha = numel(parameter.alpha);
N_r = numel(parameter.r);
P = N_diffrate * N_temporal * N_numClust * N_alpha * N_r;
switch dataset
    case {'3day', 1}
        r = zeros(P, 1);
        p = zeros(P, 1);
    case {'6week', 2}
        r = zeros(P, 3);
        p = zeros(P, 3);
end

% Loop over the table and fill in one entry at each time by calling the
% subroutine
%
for P_idx = 1:P
    % P_idx are ordered first by r, then by alpha, then numClust,
    % then temporal, finally diffrate
    %
    crt.r = parameter.r( mod(P_idx - 1, N_r) + 1 );
    crt.alpha = parameter.alpha( floor( mod(P_idx - 1, N_r * N_alpha) / N_r ) + 1 );
    crt.numClust = parameter.numClust( floor( mod(P_idx - 1, ...
        N_r * N_alpha * N_numClust) / (N_r * N_alpha) ) + 1 );
    crt.temporal = parameter.temporal( floor( mod(P_idx - 1, ...
        N_r * N_alpha * N_numClust * N_temporal) / (N_r * N_alpha * N_numClust) ) + 1 );
    diff_rate_crt = diff_rate( floor( (P_idx - 1) / (N_r * N_alpha * N_numClust * N_temporal) ) + 1 );
    switch dataset
        case {'3day', 1}
            [r(P_idx, :), p(P_idx, :)] = compute_correlation_3day(diff_rate_crt, measure, crt, permute);
        case {'6week', 2}
            [r(P_idx, :), p(P_idx, :)] = compute_correlation_6week(diff_rate_crt, measure, crt, permute);
    end
end

end