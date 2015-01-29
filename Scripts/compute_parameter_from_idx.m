function [ crt ] = compute_parameter_from_idx( measure, parameter, idx )
%The function to compute specific parameters such as temporal, numClust,
%alpha, and r given a range of parameter want to try and the idx yielding
%the best correlation. Useful to find best parameter combinations.
%
%  INPUT: measure    - a scalar integer 1 - 7 corresponding to
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
%         idx        - idx in the all range ordered from temporal,
%                      numClust, alpha, and r
% OUTPUT: r          - correlation between measures and learning 
%         p          - p-values between measures and learning
%
%---------------------------------------------------------------------
% coded by Weiyu Huang,  whuang@seas.upenn.edu
% Ver 1 : Jan 15, 2014
%---------------------------------------------------------------------

%--------------------------------------------------------------------%
%------------------ 1. check necessary parameters -------------------%
%--------------------------------------------------------------------%

% Note: 1. if unnecessary parameters exist, set them to -1
%       2. temporal 0 meaning default
%       3. if r is not defined, 0.2 * std is the most commonly used one
%
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
%--------------- 2. compute best parameters given idx ---------------%
%--------------------------------------------------------------------%

N_temporal = numel(parameter.temporal);
N_numClust = numel(parameter.numClust);
N_alpha = numel(parameter.alpha);
N_r = numel(parameter.r);
crt.r = parameter.r( mod(idx - 1, N_r) + 1 );
crt.alpha = parameter.alpha( floor( mod(idx - 1, N_r * N_alpha) / N_r ) + 1 );
crt.numClust = parameter.numClust( floor( mod(idx - 1, ...
    N_r * N_alpha * N_numClust) / (N_r * N_alpha) ) + 1 );
crt.temporal = parameter.temporal( floor( (idx - 1) / ...
    (N_r * N_alpha * N_numClust) ) + 1 );

%--------------------------------------------------------------------%
%----------------- 3. remove unnecessary parameters -----------------%
%--------------------------------------------------------------------%

switch measure
    case 1
        crt = rmfield(crt, 'alpha'); 
        crt = rmfield(crt, 'r');
    case {2, 4}
        crt = rmfield(crt, 'r');
    case {3, 5}
        crt = rmfield(crt, 'temporal');
        crt = rmfield(crt, 'numClust');
        crt = rmfield(crt, 'r');
    case {6, 7}
        crt = rmfield(crt, 'alpha');
end

end

