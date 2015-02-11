function [ r, p ] = compute_correlation_3day(diff_rate, measure, crt_unglobal, permute)
%Compute the correlaiton for 3day learning set for the measure MEASURE with
%the specific parameter defined in CRT and brain signal diffused with
%diffusion rate DIFF_RATE. Here we do a permutation test such that the
%learning rate for subjects are permuted wrt random seed PERMUTE.
%
%  INPUT: diff_rate  - diffusion rate, 0 if not diffused
%         measure    - a scalar integer 1 - 7 corresponding to
%                         1. number of switches
%                         2. number of switches greater than a threshold
%                         3. number of brain signal differences > threshold
%                         4. weighted total switches after clustering
%                         5. weighted total brain signal differences
%                         6. approximate entropy of the clustering results
%                         7. sample entropy of the clustering results
%         crt        - potential parameters
%                         .temporal - temporal normalizing factor,
%                                     used in measures 1,2,4,6,7
%                         .numClust - number of clusters,
%                                     used in measures 1,2,4,6,7
%                         .alpha    - threshold in measures 2,3
%                                     weight in measures 4,5
%                         .r        - tolerance in entropy measures 6,7
%                                     (typically 0.2 * std)
%         permute    - 'crt' if subjects' learnign rates permuted with 
%                      current random generator OR an integer to denote the
%                      specific permutation using the integer as the random
%                      seed
% OUTPUT: r          - correlation between measures and learning 
%         p          - p-values between measures and learning
%
%---------------------------------------------------------------------
% coded by Weiyu Huang,  whuang@seas.upenn.edu
% Ver 1  : Jan 14, 2014
% Ver 1.1: Feb 4, 2014     % remove one for loop for subjects
%---------------------------------------------------------------------


%--------------------------------------------------------------------%
%---------------------- 1. load learning data -----------------------%
%--------------------------------------------------------------------%

temp = load('../../learning_v1_Data/3_day_learning/behavioral_data.mat');
learning = -temp.mvt_3sess(:, 1); clear mvt_3sess;
global goodones; goodones = temp.goodones;
global diff_rate_global; diff_rate_global = diff_rate;
global crt; crt = crt_unglobal;

% Initializing storage space for each subject.
%
measure_results = zeros(length(goodones), 1);

%--------------------------------------------------------------------%
%-------------------- 2. load necessary results ---------------------%
%--------------------------------------------------------------------%
%
% Storage examples:
% pairwise distance      - diffRate0.1.mat
%                          dist_diff
% clustering results     - diffRate0.1_temporal6.mat
%                          Z, dist_diff
% cross-cluster distance - diffRate0.1_temporal6_numClust10.mat
%                          CC, T
%

switch measure
    % For measures 3 and 5, we need pairwise distance between pairs of
    % time points, first check if pairwise distance for this diff_rate
    % has been computed and load it if so, otherwise compute it.
    % Actually, we only need diffused signal, but using pairwise
    % distance can share computation results as other measures.
    %
    case {3, 5}
        if exist(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat'], 'file')
            load(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat']);
        else
            data = load('../Data_matlab/3_day_data.mat');
            dist_diff = compute_pairwise_distance(data);
            save(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat'], 'dist_diff');
        end % if checking for pairwise distance
        
        % For measures 1, 6, and 7, we need not only pairwise distance
        % between pairs of time points, we also need the clustering results
        % after temporal smoothing. We proceed this by checking whether
        % clustering results exist, otherwise check if dist_diff results
        % exist, finally compute what remains to be computed.
        %
    case {1, 6, 7}
        if exist(['../Intermediate_Results/3day/diffRate', ...
                num2str(diff_rate_global), '_temporal', num2str(crt.temporal), '.mat'], 'file')
            load(['../Intermediate_Results/3day/diffRate', ...
                num2str(diff_rate_global), '_temporal', num2str(crt.temporal), '.mat'], 'Z');
        else
            % clustering results do not exist, but pairwise distance
            % matrix may exist, which is helpful. Same as measures 3
            % and 5.
            %
            if exist(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat'], 'file')
                load(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat']);
            else
                data = load('../Data_matlab/3_day_data.mat');
                dist_diff = compute_pairwise_distance(data);
                save(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat'], 'dist_diff');
            end % if for pairwise distance loading
            
            [Z, dist_diff_temporal] = compute_temporal_clustering(dist_diff);
            save(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '_temporal', num2str(crt.temporal), ...
                '.mat'], 'Z', 'dist_diff_temporal');
        end % if for clustering loading
        
        % For measures 2 and 4, we need not only pairwise distance as well
        % as clustering results, we also need pairwise average distance
        % between pairs of clusters. We proceed this by checking whether
        % pairwise distance between pairs of clusters exist, otherwise
        % check if clustering results exist, otherwise check if dist_diff
        % results exist, finally compute what remains to be computed.
        %
    case {2, 4}
        if exist(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '_temporal', num2str(crt.temporal),...
                '_numClust', num2str(crt.numClust), '.mat'], 'file')
            load(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '_temporal', num2str(crt.temporal),...
                '_numClust', num2str(crt.numClust), '.mat'], 'CC', 'T');
        else
            % Cross-cluster distance do not exist, but clustering
            % results may exist.
            %
            if exist(['../Intermediate_Results/3day/diffRate', ...
                    num2str(diff_rate_global), '_temporal', num2str(crt.temporal), '.mat'], 'file')
                load(['../Intermediate_Results/3day/diffRate', ...
                    num2str(diff_rate_global), '_temporal', num2str(crt.temporal), '.mat'], 'Z', 'dist_diff_temporal');
            else
                % Clustering results do not exist, but pairwise
                % distance between pairs of signals may exist.
                %
                if exist(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat'], 'file')
                    load(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat']);
                else
                    data = load('../Data_matlab/3_day_data.mat');
                    dist_diff = compute_pairwise_distance(data);
                    save(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '.mat'], 'dist_diff');
                end % if for pairwise distance loading
                
                [Z, dist_diff_temporal] = compute_temporal_clustering(dist_diff);
            end % if for clustering loading
            
            [T, CC] = compute_cluster_and_clusterwise_distance(Z, dist_diff_temporal);
            save(['../Intermediate_Results/3day/diffRate', num2str(diff_rate_global), '_temporal', num2str(crt.temporal),...
                '_numClust', num2str(crt.numClust), '.mat'], 'CC', 'T');
        end % if for cross-cluster distance loading
end % switch
    
%--------------------------------------------------------------------%
%----------------------- 3. compute measures ------------------------%
%--------------------------------------------------------------------%

% We use subject_id_idx to denote the position of subject_id in
% goodones, i.e. max(subject_id_idx) = 18 while max(subject_id) = 20.
%
for subject_id = goodones
    subject_id_idx = find(goodones == subject_id);
    switch measure
        case 1
            % Cluster with the numClust parameter.
            %
            T = cluster(Z{subject_id}, 'maxclust', crt.numClust);
            N = length(T);
            T_dist = diff(T);
            measure_results(subject_id_idx) = sum(~~T_dist) / N;
        case 2
            T_crt = T{subject_id};
            CC_crt = CC{subject_id};
            N = length(T_crt);
            temp = zeros(1, N - 1);
            for i = 1:N-1
                temp(i) = (CC_crt(T_crt(i), T_crt(i + 1)));
            end
            measure_results(subject_id_idx) = sum(temp > crt.alpha) / N;
        case 3
            N = length(dist_diff{subject_id});
            temp = diag(dist_diff{subject_id}, 1);
            measure_results(subject_id_idx) = sum(temp > crt.alpha) / N;
        case 4
            T_crt = T{subject_id};
            CC_crt = CC{subject_id};
            N = length(T_crt);
            temp = zeros(1, N - 1);
            for i = 1:N-1
                temp(i) = (CC_crt(T_crt(i), T_crt(i + 1)));
            end
            measure_results(subject_id_idx) = sum(temp .^ crt.alpha) / N;
        case 5
            N = length(dist_diff{subject_id});
            temp = diag(dist_diff{subject_id}, 1);
            measure_results(subject_id_idx) = sum(temp .^ crt.alpha) / N;
        case 6
            T = cluster(Z{subject_id}, 'maxclust', crt.numClust);
            sd = std(T);
            measure_results(subject_id_idx) = ...
                entropyMeasure(2, crt.r * sd, T, 1, 'ApproximateEntropy');
        case 7
            T = cluster(Z{subject_id}, 'maxclust', crt.numClust);
            sd = std(T);
            measure_results(subject_id_idx) = entropyMeasure(2, crt.r * sd, T, 1, 'SampleEntropy');
    end % switch
end % subject_id

%--------------------------------------------------------------------%
%---------------------- 3. compute correlation ----------------------%
%--------------------------------------------------------------------%

% Want to do this wrt to the desired subject permutation.
%
if numel(permute) == length(learning);
    order = permute;
elseif ~permute
    order = 1:length(learning);
elseif strcmp(permute, 'crt') 
    order = randperm(length(learning));
else
    rng(permute);
    order = randperm(length(learning));
end 

% Compute and output correlation.
%
[R, P] = corrcoef(measure_results, learning(order));
r = R(1, 2); p = P(1, 2);

end


function dist_diff = compute_pairwise_distance(data)
global goodones;
global diff_rate_global;
dist_diff = cell(20, 1);
for subject_id = goodones
    aggregate_network = data.network{subject_id};
    aggregate_signal = data.signal{subject_id};
    L = diag(sum(aggregate_network)) - aggregate_network;
    diffusion_factor = inv(eye(size(L)) + diff_rate_global * L);
    compute_diffusion_distance = @(sig) norm(sig * diffusion_factor, 2);
    temp = zeros(2000);
    for i = 1:2000
        for j = i + 1:2000
            temp(i, j) = compute_diffusion_distance(aggregate_signal(i, :)...
                - aggregate_signal(j, :));
        end
    end
    dist_diff{subject_id} = temp + temp';
end

end 


function [Z, dist_diff_temporal] = compute_temporal_clustering(dist_diff)

global goodones;
global crt;
% Conduct clustering after temporal smoothing.
%
Z = cell(20, 1);
dist_diff_temporal = cell(20, 1);
for subject_id = goodones
    dist_diff_crt = dist_diff{subject_id};
    N = length(dist_diff_crt);
    temp = 1:N;
    temp_mean = mean(mean(dist_diff_crt));
    temporal_proximity = abs(repmat(temp, N, 1) - repmat(temp', 1, N));
    temporal_proximity = temporal_proximity .* temp_mean ./ (crt.temporal); %%%
    dist_diff_temporal{subject_id} = dist_diff_crt + temporal_proximity;
    
    % Ward clustering upon the diffused and temporal smoothed
    % signals.
    %
    dist_diff_temporal_squareform = squareform(dist_diff_temporal{subject_id});
    Z{subject_id} = linkage(dist_diff_temporal_squareform, 'ward');
end

end


function [T, CC] = compute_cluster_and_clusterwise_distance(Z, dist_diff_temporal)

global goodones;
global crt;
% Compute cross-cluster distances.
%
T = cell(20, 1);
CC = cell(20, 1);
for subject_id = goodones
    T_crt = cluster(Z{subject_id}, 'maxclust', crt.numClust);
    dist_diff_temporal_crt = dist_diff_temporal{subject_id};
    CC_crt = zeros(crt.numClust);
    for i = 1:crt.numClust
        for j = i + 1:crt.numClust
            temp = dist_diff_temporal_crt(T_crt == i, T_crt == j);
            CC_crt(i, j) = mean(mean(temp));
        end
    end
    CC{subject_id} = CC_crt + CC_crt';
    T{subject_id} = T_crt;
end
end