function [ r, p ] = compute_correlation_3day(diff_rate, measure, crt)
%Compute the correlaiton for 3day learning set for the measure MEASURE with
%the specific parameter defined in CRT and brain signal diffused with
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

load('../../learning_v1_Data/3_day_learning/behavioral_data.mat');
        learning = -mvt_3sess(:, 1); clear mvt_3sess;

% Initializing storage space for each subject.
%
measure_results = zeros(length(goodones), 1);
for subject_id = goodones
    
%--------------------------------------------------------------------%
%-------------------- 2. load necessary results ---------------------%
%--------------------------------------------------------------------%
%
% Storage examples:
% pairwise distance      - subject1_diffRate0.1.mat
%                          dist_diff
% clustering results     - subject1_diffRate0.1_temporal6.mat
%                          Z, dist_diff
% cross-cluster distance - subject1_diffRate0.1_temporal6_numClust10.mat
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
            if exist(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
                    '_diffRate', num2str(diff_rate), '.mat'], 'file')
                load(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
                    '_diffRate', num2str(diff_rate), '.mat']);
            else 
                data = load(['../../learning_v1_Data/3_day_learning/subject', ...
                    num2str(subject_id),'_activity_coherence_3day.mat']);
                
                % Take the signal as the  first 2000 time points
                % corresponding to the first day and take the average
                % network over the three day (caveat).
                %
                aggregate_signal = data.ACT(:,1:2000)';
                aggregate_network = data.COH_AVG;
                
                % Normalize the data such that the norm of brain signal
                % over all regions for any time is 1.
                %
                for t = 1:size(aggregate_signal, 1)
                    aggregate_signal(t, :) = aggregate_signal(t, :) ./ norm( aggregate_signal(t, :) );
                end
                
                % compute diffusion distances between pairs of time indexes of brain
                % activities (smoothing over network).
                %
                L = diag(sum(aggregate_network)) - aggregate_network;
                diffusion_factor = inv(eye(size(L)) + diff_rate * L);
                compute_diffusion_distance = @(sig) norm(sig * diffusion_factor, 2);
                dist_diff = zeros(2000);
                for i = 1:2000
                    for j = i + 1:2000
                        dist_diff(i, j) = compute_diffusion_distance(aggregate_signal(i, :)...
                            - aggregate_signal(j, :));
                    end
                end
                dist_diff = dist_diff + dist_diff';
                
                % Save results.
                %
                save(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
                    '_diffRate', num2str(diff_rate), '.mat'], 'dist_diff');
            end % if checking for pairwise distance
    
        % For measures 1, 6, and 7, we need not only pairwise distance
        % between pairs of time points, we also need the clustering results
        % after temporal smoothing. We proceed this by checking whether
        % clustering results exist, otherwise check if dist_diff results
        % exist, finally compute what remains to be computed.
        %
        case {1, 6, 7}
            if exist(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'], 'file')
                load(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'], 'Z');
            else
                % clustering results do not exist, but pairwise distance
                % matrix may exist, which is helpful. Same as measures 3
                % and 5.
                %
                if exist(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
                        '_diffRate', num2str(diff_rate), '.mat'], 'file')
                    load(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
                        '_diffRate', num2str(diff_rate), '.mat']);
                else
                    data = load(['../../learning_v1_Data/3_day_learning/subject', ...
                        num2str(subject_id),'_activity_coherence_3day.mat']);
                    aggregate_signal = data.ACT(:,1:2000)';
                    aggregate_network = data.COH_AVG;
                    
                    for t = 1:size(aggregate_signal, 1)
                        aggregate_signal(t, :) = aggregate_signal(t, :) ./ norm( aggregate_signal(t, :) );
                    end
                    
                    L = diag(sum(aggregate_network)) - aggregate_network;
                    diffusion_factor = inv(eye(size(L)) + diff_rate * L);
                    compute_diffusion_distance = @(sig) norm(sig * diffusion_factor, 2);
                    dist_diff = zeros(2000);
                    for i = 1:2000
                        for j = i + 1:2000
                            dist_diff(i, j) = compute_diffusion_distance(aggregate_signal(i, :)...
                                - aggregate_signal(j, :));
                            dist_diff(j, i) = dist_diff(i, j);
                        end
                    end
                 
                    save(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
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
                save(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
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
            if exist(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal),...
                    '_numClust', num2str(crt.numClust), '.mat'], 'file')
                load(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal),...
                    '_numClust', num2str(crt.numClust), '.mat'], 'CC', 'T');
            else
                % Cross-cluster distance do not exist, but clustering
                % results may exist.
                %
                if exist(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'], 'file')
                    load(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal), '.mat'],...
                    'Z', 'dist_diff_temporal');
                else
                    % Clustering results do not exist, but pairwise
                    % distance between pairs of signals may exist.
                    %
                    if exist(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
                            '_diffRate', num2str(diff_rate), '.mat'], 'file')
                        load(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
                            '_diffRate', num2str(diff_rate), '.mat']);
                    else
                        data = load(['../../learning_v1_Data/3_day_learning/subject', ...
                            num2str(subject_id),'_activity_coherence_3day.mat']);
                        aggregate_signal = data.ACT(:,1:2000)';
                        aggregate_network = data.COH_AVG;
                        
                        for t = 1:size(aggregate_signal, 1)
                            aggregate_signal(t, :) = aggregate_signal(t, :) ./...
                                norm( aggregate_signal(t, :) );
                        end
                        
                        L = diag(sum(aggregate_network)) - aggregate_network;
                        diffusion_factor = inv(eye(size(L)) + diff_rate * L);
                        compute_diffusion_distance = @(sig) norm(sig * diffusion_factor, 2);
                        dist_diff = zeros(2000);
                        for i = 1:2000
                            for j = i + 1:2000
                                dist_diff(i, j) = compute_diffusion_distance(aggregate_signal(i, :)...
                                    - aggregate_signal(j, :));
                                dist_diff(j, i) = dist_diff(i, j);
                            end
                        end
                 
                        save(['../Intermediate_Results/3day/subject', num2str(subject_id), ...
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
                    save(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
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
                save(['../Intermediate_Results/3day/subject', num2str(subject_id), '_diffRate', ...
                    num2str(diff_rate), '_temporal', num2str(crt.temporal),...
                    '_numClust', num2str(crt.numClust), '.mat'], 'CC', 'T');
            end % if for cross-cluster distance loading
    end % switch
    
%--------------------------------------------------------------------%
%----------------------- 3. compute measures ------------------------%
%--------------------------------------------------------------------%

    % We use subject_id_idx to denote the position of subject_id in
    % goodones, i.e. max(subject_id+idx) = 18 while max(subject_id) = 20.
    %
    subject_id_idx = find(goodones == subject_id);
    switch measure;
        case 1
            % Cluster with the numClust parameter.
            %
            T = cluster(Z, 'maxclust', crt.numClust);
            N = length(T);
            T_dist = diff(T);
            measure_results(subject_id_idx) = sum(~~T_dist) / N;
        case 2
            N = length(T);
            temp = zeros(1, N - 1);
            for i = 1:N-1
                temp(i) = (CC(T(i), T(i + 1)));
            end
            measure_results(subject_id_idx) = sum(temp > crt.alpha) / N;
        case 3
            N = length(dist_diff);
            temp = diag(dist_diff, 1);
            measure_results(subject_id_idx) = sum(temp > crt.alpha) / N;
        case 4
            N = length(T);
            temp = zeros(1, N - 1);
            for i = 1:N-1
                temp(i) = (CC(T(i), T(i + 1)));
            end
            measure_results(subject_id_idx) = sum(temp .^ crt.alpha) / N;
        case 5
            N = length(dist_diff);
            temp = diag(dist_diff, 1);
            measure_results(subject_id_idx) = sum(temp .^ crt.alpha) / N;
        case 6
            T = cluster(Z, 'maxclust', crt.numClust);
            sd = std(T);
            measure_results(subject_id_idx) = ...
                entropyMeasure(2, crt.r * sd, T, 1, 'ApproximateEntropy');
        case 7
            T = cluster(Z, 'maxclust', crt.numClust);
            sd = std(T);
            measure_results(subject_id_idx) = entropyMeasure(2, crt.r * sd, T, 1, 'SampleEntropy');
    end % switch
end % for subject_id

%--------------------------------------------------------------------%
%---------------------- 3. compute correlation ----------------------%
%--------------------------------------------------------------------%

[R, P] = corrcoef(measure_results, learning);
r = R(1, 2); p = P(1, 2);

end