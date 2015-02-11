%% trial with new functions defined
%
% first verify: compute_correlation_3day

diff_rate = 1;
measure = 1;
crt.numClust = 10;
crt.alpha = -1;
crt.r = -1;
crt.temporal = 600;
[ r, p ] = compute_correlation_3day(diff_rate, measure, crt)

%% find best parameter combination

idx = find(r > 0.9 * max(r));
crtParameter = compute_parameter_from_idx(measure, parameter, idx)

%% 3day with measure 1

clear all; close all; clc;
diff_rate = 10;
measure = 1;
permute = 0;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 3);
[r, p] = main_correlation( '3day', diff_rate, measure, parameter, 0 );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '_permute', num2str(permute), '.mat'], 'r', 'p');

%% 3day with measure 2

clear all; close all; clc;
diff_rate = 10;
measure = 2;
permute = 0;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = logspace(-4, 0, 20);
[r, p] = main_correlation( '3day', diff_rate, measure, parameter, 0 );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '_permute', num2str(permute), '.mat'], 'r', 'p');

%% 3day with measure 3

clear all; close all; clc;
diff_rate = 0;
measure = 3;
permute = 0;
parameter.alpha = logspace(-4, 0, 100);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '_permute', num2str(permute), '.mat'], 'r', 'p');

%% 3day with measure 4

clear all; close all; clc;
diff_rate = 0;
measure = 4;
permute = 0;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = linspace(0.1, 4, 20);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '_permute', num2str(permute), '.mat'], 'r', 'p');

%% 3day with measure 5

clear all; close all; clc;
diff_rate = 0;
measure = 5;
permute = 0;
parameter.alpha = linspace(0.1, 4, 100);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '_permute', num2str(permute), '.mat'], 'r', 'p');

%% 3day with measure 6

clear all; close all; clc;
diff_rate = 0;
measure = 6;
permute = 0;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '_permute', num2str(permute), '.mat'], 'r', 'p');

%% 3day with measure 7

clear all; close all; clc;
diff_rate = 0;
measure = 7;
permute = 0;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '_permute', num2str(permute), '.mat'], 'r', 'p');


%% 3day with measure 1 permutation

clear all; close all; clc;
diff_rate = 10;
measure = 1;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 1:10;
parameter.temporal = linspace(100, 1000, 10);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter, permute );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 2 permutation

clear all; close all; clc;
diff_rate = 10;
measure = 2;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = logspace(-4, 0, 20);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter, permute );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 3 permutation

clear all; close all; clc;
diff_rate = 10;
measure = 3;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.alpha = logspace(-4, 0, 100);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter, permute );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 4 permutation

clear all; close all; clc;
diff_rate = 10;
measure = 4;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = linspace(0.1, 4, 20);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter, permute );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 5 permutation

clear all; close all; clc;
diff_rate = 10;
measure = 5;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.alpha = linspace(0.1, 4, 100);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter, permute );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 6 permutation

clear all; close all; clc;
diff_rate = 10;
measure = 6;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter, permute );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 7 permutation

clear all; close all; clc;
diff_rate = 10;
measure = 7;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter, permute );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');


%% compute all measures given a permutation

clear all; close all;
for permute_seed = 100
    if permute_seed == 100
        permute = 0;
    else
        rng(permute_seed); permute = randperm(18);
    end
    r_all = cell(7, 1);
    p_all = cell(7, 1);
    diff_rate = 0;
    
    for measure = 1:7
        switch measure
            case 1
                parameter.numClust = 2:11;
                parameter.temporal = linspace(100, 1000, 10);
            case 2
                parameter.numClust = 2:11;
                parameter.temporal = linspace(100, 1000, 10);
                parameter.alpha = logspace(-4, 0, 10);
            case 3
                parameter.alpha = logspace(-4, 0, 100);
            case 4
                parameter.numClust = 2:11;
                parameter.temporal = linspace(100, 1000, 10);
                parameter.alpha = linspace(0.1, 4, 10);
            case 5
                parameter.alpha = linspace(0.1, 4, 100);
            case {6, 7}
                parameter.numClust = 2:11;
                parameter.temporal = linspace(100, 1000, 10);
        end
        [r_all{measure}, p_all{measure}] = main_correlation( '3day', diff_rate, measure, parameter, permute );
    end
    
    if permute_seed == 100
        save(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat'], 'r_all', 'p_all');
    else
        save(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permute_seed), '.mat'],...
            'r_all', 'p_all');
    end
end

%% compute median

r_median = [];
r_max = [];
for measure = 1:7
    temp = r_all{measure};
    r_median = [r_median, median(temp(~isnan(temp)))];
    r_max = [r_max, max(temp(~isnan(temp)))];
end
disp([r_median', r_max'])

%% boxplot everything

permuteSeed_all = [1 2 3 4 5 13 33 233 330 333];
everything = zeros(10, 7);
for permuteSeed = permuteSeed_all
    idx = find(permuteSeed_all == permuteSeed);
    load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
    for measure = 1:7
        temp = r_all{measure};
        everything(idx, measure) = median(temp(~isnan(temp)));
    end
end

load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
noPermute = zeros(7, 1);
for measure = 1:7
    temp = r_all{measure};
    noPermute(measure) = median(temp(~isnan(temp)));
end

boxplot(everything); grid on; hold on;
plot(1:7, noPermute);
