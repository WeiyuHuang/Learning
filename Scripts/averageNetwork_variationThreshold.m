function [ shifted_variation_distance, weighted_shifted_variation_distance ] =...
    averageNetwork_variationThreshold(...
    subject_id, diffusion_rate, threshold )
% compute number of shifted variation less than threshold at first scan
% training for a given subject of particular type of training module (1 -
% EXT, 2 - MOD, 3 - MIN), use diffusion distances by averaging all brain
% networks in the first scan training with diffusion_rate



if exist(['../quantifiedVariation/subject',num2str(subject_id),'_firstTrial_alpha',...
        num2str(diffusion_rate), '.mat'], 'file')
    % already computed.
    %
    load(['../quantifiedVariation/subject',num2str(subject_id),'_firstTrial_alpha',...
        num2str(diffusion_rate), '.mat'], 'signal_diff');
    
else
    % load data
    %
    data = load(['../Activity_Connectivity/subject',num2str(subject_id),'_activity_coherence_3day.mat']);
    
    aggregate_signal = data.ACT(:,1:2000)';
    aggregate_network = data.COH_AVG;
    
    % time zscore
    %
    aggregate_signal = aggregate_signal - repmat(mean(aggregate_signal), 2000, 1);
    
    % time normalizing
    %
    for r = 1:size(aggregate_signal, 1)
        aggregate_signal(r, :) = aggregate_signal(r, :) ./ norm( aggregate_signal(r, :) );
    end
    
    % compute diffused brain signal using the underlying brain-brain
    % network (smoothing over network)
    %
    L = diag(sum(aggregate_network)) - aggregate_network;
    diffusion_factor = inv(eye(size(L)) + diffusion_rate * L);
    compute_diffusion_signal = @(sig) sig * diffusion_factor;
    signal_diff = zeros(size(aggregate_signal));
    for i = 1:size(signal_diff, 1)
        signal_diff(i,:) = compute_diffusion_signal(aggregate_signal(i,:));
    end
    
    % save results
    %
    save(['../quantifiedVariation/subject',num2str(subject_id),'_firstTrial_alpha',...
        num2str(diffusion_rate), '.mat'], 'signal_diff');
end

% quantify variation using shifted variation distance
% 
temp = zeros(size(signal_diff, 1)-1, 1);
for i = 2:size(signal_diff, 1)
    temp(i-1) = norm(signal_diff(i) - signal_diff(i-1), 2);
end
shifted_variation_distance = sum(temp > threshold) ./ length(temp);

weighted_shifted_variation_distance = sum(temp(temp > threshold)) ./ length(temp);

end

