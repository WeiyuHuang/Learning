function [ optimal_numClust, optimal_temporalFactor ] = ...
    optimal_numClust_temporalFactor( trial_idx, temporalFactor_range, numClust_range )

optimal_temporalFactor = temporalFactor_range(mod((trial_idx-1), length(temporalFactor_range)) + 1);
optimal_numClust = numClust_range(ceil(trial_idx / length(temporalFactor_range)));

end

