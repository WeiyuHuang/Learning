
function [ shifted_variation_distance ] = averageNetwork_windowVariation(...
    subject_id, diffusion_rate, GAMMA, WINDOW )
% compute number of switches at first scan training for a given subject of
% particular type of training module (1 - EXT, 2 - MOD, 3 - MIN), use
% diffusion distances by averaging all brain networks in the first scan
% training with diffusion_rate


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

T = size(signal_diff, 1);
shifted_variation_distance = 0;
for w = 1:WINDOW
    temp = zeros(T - w, 1);
    for t_idx = 1:(T - w)
        t = t_idx + w;
        temp(t_idx) = norm(signal_diff(t) - signal_diff(t-w), 2) .^ 2;
    end
    shifted_variation_distance = shifted_variation_distance + sum(temp) .* (GAMMA .^ (w-1));
end
shifted_variation_distance = shifted_variation_distance ./ (T-1);

end
