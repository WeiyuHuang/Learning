function [ numSwitch ] = averageNetwork_numSwitch( ...
    subject_id, diffusion_rate, parameter, clusterMethod, temporalFactor )
% compute number of switches at first scan training for a given subject of
% particular type of training module (1 - EXT, 2 - MOD, 3 - MIN), use
% diffusion distances by averaging all brain networks in the first scan
% training with diffusion_rate

if exist(['../pairwise_timeindex_difference/subject',num2str(subject_id),...
        '_firstTrial_tempralFactor', num2str(temporalFactor), ...
        '_alpha', num2str(diffusion_rate), '.mat'], 'file')
    
    % already computed.
    %
    load(['../pairwise_timeindex_difference/subject',num2str(subject_id),...
        '_firstTrial_tempralFactor', num2str(temporalFactor), ...
        '_alpha', num2str(diffusion_rate), '.mat'], 'Z');
    N = length(Z);
    
elseif exist(['../pairwise_timeindex_difference/subject',num2str(subject_id),...
        '_firstTrial_alpha', num2str(diffusion_rate), '.mat'], 'file')
    
    % partially computed
    %
    load(['../pairwise_timeindex_difference/subject',num2str(subject_id),...
        '_firstTrial_alpha', num2str(diffusion_rate), '.mat'], 'dist_diff');
    N = length(dist_diff);
    temp = 1:N;
    temp_mean = mean(mean(dist_diff));
    temporal_proximity = abs(repmat(temp, N, 1) - repmat(temp', 1, N));
    temporal_proximity = temporal_proximity .* temp_mean ./ (N / temporalFactor);
    dist_diff_temporal = dist_diff + temporal_proximity;
    
    % clustering
    %
    dist_diff_temporal_squareform = squareform(dist_diff_temporal);
    Z = linkage(dist_diff_temporal_squareform, 'ward');
    
    save(['../pairwise_timeindex_difference/subject',num2str(subject_id),...
        '_firstTrial_tempralFactor', num2str(temporalFactor), ...
        '_alpha', num2str(diffusion_rate), '.mat'], 'Z', 'dist_diff_temporal');
    
else
    
    % load data
    %
    data = load(['../Activity_Connectivity/subject',num2str(subject_id),'_activity_coherence_3day.mat']);
    
    aggregate_signal = data.ACT(:,1:2000)';
    aggregate_network = data.COH_AVG;
%     
    % time zscore
    %
    aggregate_signal = aggregate_signal - repmat(mean(aggregate_signal), 2000, 1);
    % region zscore
%     %
%     for r = 1:size(aggregate_signal, 2)
%         aggregate_signal(:, r) = zscore(aggregate_signal(:, r));
%     end
    % time normalizing
    %
    for r = 1:size(aggregate_signal, 1)
        aggregate_signal(r, :) = aggregate_signal(r, :) ./ norm( aggregate_signal(r, :) );
    end
    
    % compute diffusion distances between pairs of time indexes of brain
    % activities (smoothing over network)
    %
    L = diag(sum(aggregate_network)) - aggregate_network;
    diffusion_factor = inv(eye(size(L)) + diffusion_rate * L);
    compute_diffusion_distance = @(sig) norm(sig * diffusion_factor, 2);
    dist_diff = zeros(length(aggregate_signal));
    for i = 1:length(aggregate_signal)
        for j = i+1:length(aggregate_signal)
            dist_diff(i,j) = compute_diffusion_distance(aggregate_signal(i,:)...
                - aggregate_signal(j,:));
            dist_diff(j,i) = dist_diff(i,j);
        end
    end
    
    % temporal smoothing
    %
    N = length(dist_diff);
    temp = 1:N;
    temp_mean = mean(mean(dist_diff));
    temporal_proximity = abs(repmat(temp, N, 1) - repmat(temp', 1, N));
    temporal_proximity = temporal_proximity .* temp_mean ./ (N / temporalFactor);
    dist_diff_temporal = dist_diff + temporal_proximity;
    
    % clustering
    %
    dist_diff_temporal_squareform = squareform(dist_diff_temporal);
    Z = linkage(dist_diff_temporal_squareform, 'ward');
    save(['../pairwise_timeindex_difference/subject',num2str(subject_id),...
        '_firstTrial_tempralFactor', num2str(temporalFactor), ...
        '_alpha', num2str(diffusion_rate), '.mat'], 'Z', 'dist_diff_temporal', 'dist_diff');
end

switch lower(clusterMethod) % see cluster
    case 'maxclust'
        T = cluster(Z, 'Maxclust', parameter(1));
    case {'cutoff', 'cutoff-default'}
        T = cluster(Z, 'Cutoff', parameter(1));
    case 'cutoff-depth'
        T = cluster(Z, 'Cutoff', parameter(1), 'Depth', parameter(2));
    case 'cutoff-distance'
        T = cluster(Z, 'Cutoff', parameter(1), 'Criterion','distance');
    case {'cutoff-depth-distance', 'cutoff-distance-depth'}
        T = cluster(Z, 'Cutoff', parameter(1), 'Depth', parameter(2), ...
            'Criterion','distance');
end

% average number of changes
% 
T_dist = diff(T);
numSwitch = (N - sum(~T_dist)) ./ N;
end

