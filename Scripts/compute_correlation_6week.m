function [ r, p ] = compute_correlation_6week(diff_rate, measure, crt_unglobal, permute)
%Compute the correlaiton for 6week learning set for the measure MEASURE
%with the specific parameter defined in CRT and brain signal diffused with
%diffusion rate DIFF_RATE. Here we do a permutation test such that the
%learning rate for subjects are permuted wrt random seed PERMUTE.
%
%  INPUT: diff_rate  - diffusion rate, 0 if not diffused, may be a vector
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
% Ver 1.1: Feb 4, 2014
%---------------------------------------------------------------------


%--------------------------------------------------------------------%
%---------------------- 1. load learning data -----------------------%
%--------------------------------------------------------------------%

if exist('../Intermediate_Results/6_week_learning_data.mat', 'file')
    load('../Intermediate_Results/6_week_learning_data.mat');
else
    load('../../learning_v1_Data/6_week_learning/all_home_scanner_weighted.mat');
    scan_id = 1;
    learning_rate = zeros(20,6);
    for subject_id = 1:20
        for module_id = 1:6
            f = scanTrain_expFit( all_home_scanner_weighted, subject_id, module_id, scan_id, 0);
            learning_rate(subject_id, module_id) = -f.b;
        end
    end
    % We have 6 modules in 3 difficulties. However, the brain
    % signals are only categorized into 3 difficulties. So we take
    % the mean learning rates for each of the two modules
    % corresponding to each of the 3 difficulties as the learning
    % rate for such difficulty.
    %
    ext_learning = mean(learning_rate(:,1:2), 2);
    mod_learning = mean(learning_rate(:,3:4), 2);
    min_learning = mean(learning_rate(:,5:6), 2);
    learning = [ext_learning, mod_learning, min_learning];
    save('../Intermediate_Results/6_week_learning_data.mat', 'learning');
end % end for checking whether learning filed exist
global goodones; goodones = 1:20;
global diff_rate_global; diff_rate_global = diff_rate;
global crt; crt = crt_unglobal;

% Initializing storage space for each subject.
%
measure_results = zeros(length(goodones), 3);
    
%--------------------------------------------------------------------%
%-------------------- 2. load necessary results ---------------------%
%--------------------------------------------------------------------%
%
% Storage examples:
% pairwise distance      - diffRate0.1.mat
%                          dist_diff
% clustering results     - diffRate0.1_temporal6.mat
%                          Z, dist_diff
% cross-cluster distance - diffRate0.1_temporal6_numClust10.
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
        if exist(['../Intermediate_Results/6week/diffRate', num2str(diff_rate), '.mat'], 'file')
            load(['../Intermediate_Results/6week/diffRate', num2str(diff_rate), '.mat']);
        else
            data = load('../Data_matlab/6_week_data.mat');
            dist_diff = compute_pairwise_distance(data);
            save(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '.mat'], 'dist_diff');
        end % if checking for pairwise distance
        
        % For measures 1, 6, and 7, we need not only pairwise distance
        % between pairs of time points, we also need the clustering results
        % after temporal smoothing. We proceed this by checking whether
        % clustering results exist, otherwise check if dist_diff results
        % exist, finally compute what remains to be computed.
        %
    case {1, 6, 7}
        if exist(['../Intermediate_Results/6week/diffRate', ...
                num2str(diff_rate_global), '_temporal', num2str(crt.temporal), '.mat'], 'file')
            load(['../Intermediate_Results/6week/diffRate', ...
                num2str(diff_rate_global), '_temporal', num2str(crt.temporal), '.mat'], 'Z');
        else
            % clustering results do not exist, but pairwise distance
            % matrix may exist, which is helpful. Same as measures 3
            % and 5.
            %
            if exist(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '.mat'], 'file')
                load(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '.mat']);
            else
                data = load('../Data_matlab/6_week_data.mat');
                dist_diff = compute_pairwise_distance(data);
                save(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '.mat'], 'dist_diff');
            end % if for pairwise distance loading
            
            [Z, dist_diff_temporal] = compute_temporal_clustering(dist_diff);
            save(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '_temporal',...
                num2str(crt.temporal), '.mat'], 'Z', 'dist_diff_temporal');
        end % if for clustering loading
        
        % For measures 2 and 4, we need not only pairwise distance as well
        % as clustering results, we also need pairwise average distance
        % between pairs of clusters. We proceed this by checking whether
        % pairwise distance between pairs of clusters exist, otherwise
        % check if clustering results exist, otherwise check if dist_diff
        % results exist, finally compute what remains to be computed.
        %
    case {2, 4}
        if exist(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '_temporal', num2str(crt.temporal),...
                '_numClust', num2str(crt.numClust), '.mat'], 'file')
            load(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '_temporal', num2str(crt.temporal),...
                '_numClust', num2str(crt.numClust), '.mat'], 'CC', 'T');
        else
            % Cross-cluster distance do not exist, but clustering
            % results may exist.
            %
            if exist(['../Intermediate_Results/6week/diffRate', ...
                    num2str(diff_rate_global), '_temporal', num2str(crt.temporal), '.mat'], 'file')
                load(['../Intermediate_Results/6week/diffRate', ...
                    num2str(diff_rate_global), '_temporal', num2str(crt.temporal), '.mat'], 'Z', 'dist_diff_temporal');
            else
                % Clustering results do not exist, but pairwise
                % distance between pairs of signals may exist.
                %
                if exist(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '.mat'], 'file')
                    load(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '.mat']);
                else
                    data = load('../Data_matlab/6_week_data.mat');
                    dist_diff = compute_pairwise_distance(data);
                    save(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '.mat'], 'dist_diff');
                end % if for pairwise distance loading
                
                [Z, dist_diff_temporal] = compute_temporal_clustering(dist_diff);
            end % if for clustering loading
            
            [T, CC] = compute_cluster_and_clusterwise_distance(Z, dist_diff_temporal);
            save(['../Intermediate_Results/6week/diffRate', num2str(diff_rate_global), '_temporal', num2str(crt.temporal),...
                '_numClust', num2str(crt.numClust), '.mat'], 'CC', 'T');
        end % if for cross-cluster distance loading
end % switch

%--------------------------------------------------------------------%
%----------------------- 3. compute measures ------------------------%
%--------------------------------------------------------------------%
r = zeros(1, 3);
p = zeros(1, 3);
for subject_id = goodones
    for type_id = 1:3
        switch measure
            case 1
                % Cluster with the numClust parameter.
                %
                T = cluster(Z{subject_id, type_id}, 'maxclust', crt.numClust);
                N = length(T);
                T_dist = diff(T);
                measure_results(subject_id, type_id) = sum(~~T_dist) / N;
            case 2
                T_crt = T{subject_id, type_id};
                CC_crt = CC{subject_id, type_id};
                N = length(T_crt);
                temp = zeros(1, N - 1);
                for i = 1:N-1
                    temp(i) = (CC_crt(T_crt(i), T_crt(i + 1)));
                end
                measure_results(subject_id, type_id) = sum(temp > crt.alpha) / N;
            case 3
                N = length(dist_diff{subject_id, type_id});
                temp = diag(dist_diff{subject_id, type_id}, 1);
                measure_results(subject_id, type_id) = sum(temp > crt.alpha) / N;
            case 4
                T_crt = T{subject_id, type_id};
                CC_crt = CC{subject_id, type_id};
                N = length(T_crt);
                temp = zeros(1, N - 1);
                for i = 1:N-1
                    temp(i) = (CC_crt(T_crt(i), T_crt(i + 1)));
                end
                measure_results(subject_id, type_id) = sum(temp .^ crt.alpha) / N;
            case 5
                N = length(dist_diff{subject_id, type_id});
                temp = diag(dist_diff{subject_id, type_id}, 1);
                measure_results(subject_id, type_id) = sum(temp .^ crt.alpha) / N;
            case 6
                T = cluster(Z{subject_id, type_id}, 'maxclust', crt.numClust);
                sd = std(T);
                measure_results(subject_id, type_id) = ...
                    entropyMeasure(2, crt.r * sd, T, 1, 'ApproximateEntropy');
            case 7
                T = cluster(Z{subject_id, type_id}, 'maxclust', crt.numClust);
                sd = std(T);
                measure_results(subject_id, type_id) = entropyMeasure(2, crt.r * sd, T, 1, 'SampleEntropy');
        end % switch
    end % type_id
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
for type_id = 1:3
    [R, P] = corrcoef(measure_results(:, type_id), learning(order, type_id));
    r(type_id) = R(1, 2); p(type_id) = P(1, 2);
end

end

%--------------------------------------------------------------------%
%-------------------------- A. subroutines --------------------------%
%--------------------------------------------------------------------%


function dist_diff = compute_pairwise_distance(data)
global goodones;
global diff_rate_global;
dist_diff = cell(20, 3);
for subject_id = goodones
    for type_id = 1:3
        aggregate_network = data.network{subject_id, type_id};
        aggregate_signal = data.signal{subject_id, type_id};
        L = diag(sum(aggregate_network)) - aggregate_network;
        diffusion_factor = inv(eye(size(L)) + diff_rate_global * L);
        compute_diffusion_distance = @(sig) norm(sig * diffusion_factor, 2);
        totalTime = size(data.signal{subject_id, type_id}, 1);
        temp = zeros(totalTime);
        for i = 1:totalTime
            for j = i + 1:totalTime
                temp(i, j) = compute_diffusion_distance(aggregate_signal(i, :) - aggregate_signal(j, :));
            end
        end
        dist_diff{subject_id, type_id} = temp + temp';
    end
end

end 


function [Z, dist_diff_temporal] = compute_temporal_clustering(dist_diff)
global goodones;
global crt;
% Conduct clustering after temporal smoothing.
%
Z = cell(20, 3);
dist_diff_temporal = cell(20, 3);
for subject_id = goodones
    for type_id = 1:3
        dist_diff_crt = dist_diff{subject_id, type_id};
        N = length(dist_diff_crt);
        temp = 1:N;
        temp_mean = mean(mean(dist_diff_crt));
        temporal_proximity = abs(repmat(temp, N, 1) - repmat(temp', 1, N));
        temporal_proximity = temporal_proximity .* temp_mean ./ (crt.temporal); %%%
        dist_diff_temporal{subject_id, type_id} = dist_diff_crt + temporal_proximity;
        
        % Ward clustering upon the diffused and temporal smoothed
        % signals.
        %
        dist_diff_temporal_squareform = squareform(dist_diff_temporal{subject_id, type_id});
        Z{subject_id, type_id} = linkage(dist_diff_temporal_squareform, 'ward');
    end
end

end


function [T, CC] = compute_cluster_and_clusterwise_distance(Z, dist_diff_temporal)
global goodones;
global crt;
% Compute cross-cluster distances.
%
T = cell(20, 3);
CC = cell(20, 3);
for subject_id = goodones
    for type_id = 1:3
        T_crt = cluster(Z{subject_id, type_id}, 'maxclust', crt.numClust);
        dist_diff_temporal_crt = dist_diff_temporal{subject_id, type_id};
        CC_crt = zeros(crt.numClust);
        for i = 1:crt.numClust
            for j = i + 1:crt.numClust
                temp = dist_diff_temporal_crt(T_crt == i, T_crt == j);
                CC_crt(i, j) = mean(mean(temp));
            end
        end
        CC{subject_id, type_id} = CC_crt + CC_crt';
        T{subject_id, type_id} = T_crt;
    end
end

end

