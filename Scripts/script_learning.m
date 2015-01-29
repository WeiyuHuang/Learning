% script to have similar results as we have in 6-week learning
%

%% compute switches using sihifted variation distance but for a window

diffusion_rate = 10;
GAMMA_vector = [0.05 0.1 0.15];
WINDOW_vector = [1:10];
num_trial = length(GAMMA_vector) * length(WINDOW_vector);
corr_switchLearning_results = zeros(num_trial, 1);
switch_results = zeros(num_trial, 1);
goodones = [1 2 3 4 5 7 8 9 10 11 12 14 15 16 17 18 19 20];

for trial_idx = 1:num_trial
    GAMMA_idx = floor( (trial_idx - 1) / length(WINDOW_vector)) + 1;
    WINDOW_idx = mod(trial_idx - 1, length(WINDOW_vector)) + 1;
    GAMMA = GAMMA_vector(GAMMA_idx);
    WINDOW = WINDOW_vector(WINDOW_idx);
    switches = zeros(length(goodones), 1);
    
    for subject_id_idx = 1:length(goodones)
        subject_id = goodones(subject_id_idx);
        temp =  averageNetwork_windowVariation(...
            subject_id, diffusion_rate, GAMMA, WINDOW );
        switches(subject_id_idx) = temp;
    end
    switch_results(trial_idx) = mean(switches(:));
    
    r = corrcoef(switches(:,1), firstDay_learning);
    corr_switchLearning_results(trial_idx,1) = r(1,2);

end

%% compute switches by counting number of adjacent time difference higher than a threshold

diffusion_rate = 10;
threshold_range = logspace(-4, 0, 100);
corr_switchLearning_results = zeros(length(threshold_range), 1);
switch_results = zeros(length(threshold_range), 1);
goodones = [1 2 3 4 5 7 8 9 10 11 12 14 15 16 17 18 19 20];

for threshold_idx = 1:length(threshold_range)
    switches = zeros(length(goodones), 1);
    
    for subject_id_idx = 1:length(goodones)
        subject_id = goodones(subject_id_idx);
        [ temp, ~ ] = averageNetwork_variationThreshold(...
            subject_id, diffusion_rate, threshold_range(threshold_idx) );
        switches(subject_id_idx) = temp;
    end
    
    switch_results(threshold_idx) = mean(switches(:));
    
    r = corrcoef(switches, firstDay_learning);
    corr_switchLearning_results(threshold_idx) = r(1,2);
end



%% compute switches using numClust and counting switches

diffusion_rate = 10;
numClust_range = 2:20;
temporalFactor_range = 1:8;
corr_switchLearning_results = zeros(length(numClust_range) * length(temporalFactor_range), 1);
switch_results = zeros(length(numClust_range) * length(temporalFactor_range), 1);
goodones = [1 2 3 4 5 7 8 9 10 11 12 14 15 16 17 18 19 20];

for trial_idx = 1:length(numClust_range) * length(temporalFactor_range)
    temporalFactor = temporalFactor_range(mod((trial_idx-1), length(temporalFactor_range)) + 1);
    numClust = numClust_range(ceil(trial_idx / length(temporalFactor_range)));
    switches = zeros(length(goodones), 1);
    
    for subject_id_idx = 1:length(goodones)
        subject_id = goodones(subject_id_idx);
        temp =  averageNetwork_numSwitch( ...
            subject_id, diffusion_rate, numClust, 'maxclust', temporalFactor );
        switches(subject_id_idx) = temp;
    end
    switch_results(trial_idx) = mean(switches(:));
    
    r = corrcoef(switches(:,1), firstDay_learning);
    corr_switchLearning_results(trial_idx,1) = r(1,2);
end

%% compute switches using numClust with entropy (approximate or sample) measures

diffusion_rate = 10;
numClust_range = 2:10;
temporalFactor_range = 3:8;
entropyMethod = 'appEntropy';
corr_switchLearning_results = zeros(length(numClust_range) * length(temporalFactor_range), 1);
switch_results = zeros(length(numClust_range) * length(temporalFactor_range), 1);
goodones = [1 2 3 4 5 7 8 9 10 11 12 14 15 16 17 18 19 20];

for trial_idx = 1:length(numClust_range) * length(temporalFactor_range)
    temporalFactor = temporalFactor_range(mod((trial_idx-1), length(temporalFactor_range)) + 1);
    numClust = numClust_range(ceil(trial_idx / length(temporalFactor_range)));
    switches = zeros(length(goodones), 1);
    
    for subject_id_idx = 1:length(goodones)
        subject_id = goodones(subject_id_idx);
        temp =  averageNetwork_entropy( ...
            subject_id, diffusion_rate, numClust, 0.2, entropyMethod, temporalFactor );
        switches(subject_id_idx) = temp;
    end
    switch_results(trial_idx) = mean(switches(:));
    
    r = corrcoef(switches(:,1), firstDay_learning);
    corr_switchLearning_results(trial_idx,1) = r(1,2);
end

%% learning data

goodones = [1 2 3 4 5 7 8 9 10 11 12 14 15 16 17 18 19 20];

behaviorData = load(['../Behavior/behavioral_data.mat'], 'mvt_3sess');

firstDay_learning = - behaviorData.mvt_3sess(:,1);

%% define a function to find best numClust and TemporalFactor

trial_idx = find(corr_switchLearning_results > 0.8 * max(corr_switchLearning_results));

[ optimal_numClust, optimal_temporalFactor ] = ...
    optimal_numClust_temporalFactor( trial_idx, temporalFactor_range, numClust_range )

