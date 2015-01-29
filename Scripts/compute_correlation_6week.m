function [ r, p ] = compute_correlation_6week(type_id, diff_rate, measure, crt)
%Compute the correlaiton for 6week learning set for the measure MEASURE
%with the specific parameter defined in CRT and brain signal diffused with
%diffusion rate DIFF_RATE.
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
% OUTPUT: r          - correlation between measures and learning 
%         p          - p-values between measures and learning
%
%---------------------------------------------------------------------
% coded by Weiyu Huang,  whuang@seas.upenn.edu
% Ver 1 : Jan 14, 2014
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

% Initializing storage space for each subject.
%
measure_results = zeros(20, 1);
for subject_id = 1:20
    
%--------------------------------------------------------------------%
%-------------------- 2. load necessary results ---------------------%
%--------------------------------------------------------------------%
%
% Storage examples:
% pairwise distance      - subject1_type2_diffRate0.1.mat
%                          dist_diff
% clustering results     - subject1_type2_diffRate0.1_temporal6.mat
%                          Z, dist_diff
% cross-cluster distance - subject1_type2_diffRate0.1_temporal6_numClust10.
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
            if exist(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                    '_diffRate', num2str(diff_rate), '.mat'], 'file')
                load(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                    '_diffRate', num2str(diff_rate), '.mat']);
            else 
                dist_diff = compute_distDiff_loading(subject_id, type_id, diff_rate);
                
                % Save results.
                %
                save(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                    '_diffRate', num2str(diff_rate), '.mat'], 'dist_diff');
            end % if checking for pairwise distance
    
        % For measures 1, 6, and 7, we need not only pairwise distance
        % between pairs of time points, we also need the clustering results
        % after temporal smoothing. We proceed this by checking whether
        % clustering results exist, otherwise check if dist_diff results
        % exist, finally compute what remains to be computed.
        %
        case {1, 6, 7}
            if exist(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'], 'file')
                load(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'], 'Z');
            else
                % clustering results do not exist, but pairwise distance
                % matrix may exist, which is helpful. Same as measures 3
                % and 5.
                %
                if exist(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                        '_diffRate', num2str(diff_rate), '.mat'], 'file')
                    load(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                        '_diffRate', num2str(diff_rate), '.mat']);
                else
                    dist_diff = compute_distDiff_loading(subject_id, type_id, diff_rate);
                    
                    save(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                        '_diffRate', num2str(diff_rate), '.mat'], 'dist_diff');
                end % if for pairwise distance loading
                
                % Conduct clustering after temporal smoothing.
                %
                N = length(dist_diff);
                temp = 1:N;
                temp_mean = mean(mean(dist_diff));
                temporal_proximity = abs(repmat(temp, N, 1) - repmat(temp', 1, N));
                temporal_proximity = temporal_proximity .* temp_mean ./ (crt.temporal); %%%
                dist_diff_temporal = dist_diff + temporal_proximity;
                
                % Ward clustering upon the diffused and temporal smoothed
                % signals.
                %
                dist_diff_temporal_squareform = squareform(dist_diff_temporal);
                Z = linkage(dist_diff_temporal_squareform, 'ward');
                save(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), ...
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
            if exist(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal),...
                    '_numClust', num2str(crt.numClust), '.mat'], 'file')
                load(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal),...
                    '_numClust', num2str(crt.numClust), '.mat'], 'CC', 'T');
            else
                % Cross-cluster distance do not exist, but clustering
                % results may exist.
                %
                if exist(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'], 'file')
                    load(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'],...
                    'Z', 'dist_diff_temporal');
                else
                    % Clustering results do not exist, but pairwise
                    % distance between pairs of signals may exist.
                    %
                    if exist(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                            '_diffRate', num2str(diff_rate), '.mat'], 'file')
                        load(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                            '_diffRate', num2str(diff_rate), '.mat']);
                    else
                        dist_diff = compute_distDiff_loading(subject_id, type_id, diff_rate);
                        
                        save(['../Intermediate_Results/6week/subject', num2str(subject_id), ...
                            '_diffRate', num2str(diff_rate), '.mat'], 'dist_diff');
                    end % if for pairwise distance loading
                    
                    N = length(dist_diff);
                    temp = 1:N;
                    temp_mean = mean(mean(dist_diff));
                    temporal_proximity = abs(repmat(temp, N, 1) - repmat(temp', 1, N));
                    temporal_proximity = temporal_proximity .* temp_mean ./ (crt.temporal); %%%
                    dist_diff_temporal = dist_diff + temporal_proximity;
                    
                    dist_diff_temporal_squareform = squareform(dist_diff_temporal);
                    Z = linkage(dist_diff_temporal_squareform, 'ward');
                    save(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                        num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'], ...
                        'Z', 'dist_diff_temporal');
                end % if for clustering loading
                
                % Compute cross-cluster distances.
                %
                T = cluster(Z, 'maxclust', crt.numClust);
                N = length(T);
                CC = zeros(crt.numClust);
                for i = 1:crt.numClust
                    for j = i + 1:crt.numClust
                        temp = dist_diff_temporal(T == i, T == j);
                        CC(i, j) = mean(mean(temp));
                    end
                end
                CC = CC + CC';
                save(['../Intermediate_Results/6week/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal),...
                    '_numClust', num2str(crt.numClust), '.mat'], 'CC', 'T');
            end % if for cross-cluster distance loading
    end % switch
    
%--------------------------------------------------------------------%
%----------------------- 3. compute measures ------------------------%
%--------------------------------------------------------------------%
    switch measure
        case 1
            % Cluster with the numClust parameter.
            %
            T = cluster(Z, 'maxclust', crt.numClust);
            N = length(T);
            T_dist = diff(T);
            measure_results(subject_id) = sum(~~T_dist) / N;
        case 2
            N = length(T);
            temp = zeros(1, N - 1);
            for i = 1:N-1
                temp(i) = (CC(T(i), T(i + 1)));
            end
            measure_results(subject_id) = sum(temp > crt.alpha) / N;
        case 3
            N = length(dist_diff);
            temp = diag(dist_diff, 1);
            measure_results(subject_id) = sum(temp > crt.alpha) / N;
        case 4
            N = length(T);
            temp = zeros(1, N - 1);
            for i = 1:N-1
                temp(i) = (CC(T(i), T(i + 1)));
            end
            measure_results(subject_id) = sum(temp .^ crt.alpha) / N;
        case 5
            N = length(dist_diff);
            temp = diag(dist_diff, 1);
            measure_results(subject_id) = sum(temp .^ crt.alpha) / N;
        case 6
            T = cluster(Z, 'maxclust', crt.numClust);
            sd = std(T);
            measure_results(subject_id) = ...
                entropyMeasure(2, crt.r * sd, T, 1, 'ApproximateEntropy');
        case 7
            T = cluster(Z, 'maxclust', crt.numClust);
            sd = std(T);
            measure_results(subject_id) = entropyMeasure(2, crt.r * sd, T, 1, 'SampleEntropy');
    end % switch
end % for subject_id

%--------------------------------------------------------------------%
%---------------------- 3. compute correlation ----------------------%
%--------------------------------------------------------------------%

[R, P] = corrcoef(measure_results, learning(:, type_id));
r = R(1, 2); p = P(1, 2);

end

%--------------------------------------------------------------------%
%----------------- A. subroutine compute dist_diff ------------------%
%--------------------------------------------------------------------%

function dist_diff = compute_distDiff_loading(subject_id, type_id, diff_rate)
    
% Take the signal as the signals for the scanning period of the first
% session and the network by averaing all recorded coherence information
% over that time.
%
act = load(['../../learning_v1_Data/6_week_learning/subject',...
    num2str(subject_id), 'wav_sess1_bytype.mat']);
network = load(['../../learning_v1_Data/6_week_learning/subject',...
    num2str(subject_id),'coh_sess1_bytype.mat']);

% Select the appropriate learning intensity.
%
switch type_id % (1 - EXT, 2 - MOD, 3 - MIN)
    case 1
        crt_act_blocks = act.Af;
        crt_network_blocks = network.Af;
    case 2
        crt_act_blocks = act.Ao;
        crt_network_blocks = network.Ao;
    case 3
        crt_act_blocks = act.Ar;
        crt_network_blocks = network.Ar;
end

% Aggregate brain activity profiles
%
aggregate_signal = [];
aggregate_length = [];
for block = 1:numel(crt_act_blocks)
    crt_act = crt_act_blocks{block};
    aggregate_signal = [aggregate_signal; crt_act];
    aggregate_length = [aggregate_length; size(crt_act, 1)];
end % block

% Normalize the data such that the norm of brain signal over all regions
% for any time is 1.
%
for t1 = 1:sum(aggregate_length)
    aggregate_signal(t1,:) = aggregate_signal(t1,:) ./ norm(aggregate_signal(t1,:));
end

% Aggregate brain networks over the time for the session.
%
aggregate_network = zeros(112);
for block = 1:numel(crt_network_blocks)
    crt_network = crt_network_blocks{block};
    aggregate_network = aggregate_network + crt_network;
end
aggregate_network = aggregate_network ./ numel(crt_network_blocks);

% Compute diffusion distances between pairs of time indexes of brain
% activities (smoothing over network).
%
L = diag(sum(aggregate_network)) - aggregate_network;
diffusion_factor = inv(eye(size(L)) + diff_rate * L);
compute_diffusion_distance = @(sig) norm(sig * diffusion_factor, 2);
dist_diff = zeros(sum(aggregate_length));
for i = 1:sum(aggregate_length)
    for j = i + 1:sum(aggregate_length)
        dist_diff(i, j) = compute_diffusion_distance(aggregate_signal(i, :)...
            - aggregate_signal(j, :));
    end
end
dist_diff = dist_diff + dist_diff';

end % end of subroutine A.

