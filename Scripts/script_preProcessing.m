% pre-processing scripts

%% pre processing for 3 day dataset

load ../../learning_v1_Data/3_day_learning/behavioral_data.mat
signal = cell(length(goodones), 1);
network = cell(length(goodones), 1);

for subject_id = goodones
    data = load(['../../learning_v1_Data/3_day_learning/subject', num2str(subject_id), '_activity_coherence_3day.mat']);
    temp = data.ACT(:,1:2000)';
    for t = 1:2000
        temp(t, :) = temp(t, :) ./ norm( temp(t, :) );
    end
    signal{subject_id} = temp;
    network{subject_id} = data.COH_AVG;
end

save('../Data_matlab/3_day_data', 'signal', 'network');

%% pre processing for 3 day dataset

goodones = 1:20;
signal = cell(length(goodones), 3);
network = cell(length(goodones), 3);

for subject_id = goodones
    for type_id = 1:3
        
        act_crt = load(['../../learning_v1_Data/6_week_learning/subject', num2str(subject_id), 'wav_sess1_bytype.mat']);
        network_crt = load(['../../learning_v1_Data/6_week_learning/subject', num2str(subject_id), 'coh_sess1_bytype.mat']);
        
        % Select the appropriate learning intensity.
        %
        switch type_id % (1 - EXT, 2 - MOD, 3 - MIN)
            case 1
                crt_act_blocks = act_crt.Af;
                crt_network_blocks = network_crt.Af;
            case 2
                crt_act_blocks = act_crt.Ao;
                crt_network_blocks = network_crt.Ao;
            case 3
                crt_act_blocks = act_crt.Ar;
                crt_network_blocks = network_crt.Ar;
        end
        
        % Aggregate brain activity profiles
        %
        aggregate_signal = [];
        aggregate_length = [];
        for block = 1:numel(crt_act_blocks)
            temp = crt_act_blocks{block};
            aggregate_signal = [aggregate_signal; temp];
            aggregate_length = [aggregate_length; size(temp, 1)];
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
        
        % Save results.
        signal{subject_id, type_id} = aggregate_signal;
        network{subject_id, type_id} = aggregate_network;
        
    end % type_id
end % subject_id

save('../Data_matlab/6_week_data', 'signal', 'network');
