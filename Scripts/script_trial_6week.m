
%% median correlation 

median_results = zeros(1, 3);
max_results = zeros(1, 3);
for i = 1:3
    r_temp = r(:, i);
    median_results(i) = median(r_temp(~isnan(r_temp)));
    max_results(i) = max(r_temp(~isnan(r_temp)));
end
disp([max_results; median_results]);

%% find best parameter combination
crtParameter = cell(1, 3);
for i = 1:3
    r_temp = r(:, i);
    idx = find(r_temp >= 0.7 * max(r_temp));
    crtParameter{i} = compute_parameter_from_idx(measure, parameter, idx);
    crtParameter{i}.difficulty = i;
    disp(crtParameter{i});
end

%% 6week with measure 1

clear all; close all; clc;
diff_rate = 10;
measure = 1;
permute = 0;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
[r, p] = main_correlation( '6week', diff_rate, measure, parameter, 0 );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 2

clear all; close all; clc;
diff_rate = 0;
measure = 2;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = logspace(-4, 0, 20);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 3

clear all; close all; clc;
diff_rate = 0;
measure = 3;
parameter.alpha = logspace(-4, 0, 100);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 4

clear all; close all; clc;
diff_rate = 0;
measure = 4;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = linspace(0.1, 4, 20);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 5

clear all; close all; clc;
diff_rate = 0;
measure = 5;
parameter.alpha = linspace(0.1, 100, 200);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 6

clear all; close all; clc;
diff_rate = 0;
measure = 6;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 7

clear all; close all; clc;
diff_rate = 0;
measure = 7;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');


%% 6week with measure 1 subject permutation

clear all; close all; clc;
diff_rate = 10;
measure = 1;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 1:10;
parameter.temporal = linspace(100, 1000, 10);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter, permute );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 2 subject permutation

clear all; close all; clc;
diff_rate = 10;
measure = 2;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = logspace(-4, 0, 20);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter, permute );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 3 subject permutation

clear all; close all; clc;
diff_rate = 10;
measure = 3;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.alpha = logspace(-4, 0, 100);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter, permute );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 4 subject permutation

clear all; close all; clc;
diff_rate = 10;
measure = 4;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = linspace(0.1, 4, 20);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter, permute );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 5 subject permutation

clear all; close all; clc;
diff_rate = 10;
measure = 5;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.alpha = linspace(0.1, 100, 200);
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter, permute );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 6 subject permutation

clear all; close all; clc;
diff_rate = 10;
measure = 6;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter, permute );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% 6week with measure 7 subject permutation

clear all; close all; clc;
diff_rate = 10;
measure = 7;
rngSeed = 13; rng(rngSeed); permute.learningOrder = 1;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '6week', 1:3, diff_rate, measure, parameter, permute );

%%
save(['../Results/6week_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_rngSeed', num2str(rngSeed), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% compute all measures given a permutation

% clear all; close all;
for permute_seed = 100
    if permute_seed == 100
        permute = 0;
    else
        rng(permute_seed); permute = randperm(20);
    end
    r_all = cell(7, 3);
    p_all = cell(7, 3);
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
                parameter.alpha = logspace(-3, 2, 100);
            case {6, 7}
                parameter.numClust = 2:11;
                parameter.temporal = linspace(100, 1000, 10);
        end
        [r, p] = main_correlation( '6week', diff_rate, measure, parameter, permute );
        for type_id = 1:3
            r_all{measure, type_id} = r(:, type_id);
            p_all{measure, type_id} = p(:, type_id);
        end
    end
    
    if permute_seed == 100
        save(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat'], 'r_all', 'p_all');
    else
        save(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permute_seed), '.mat'],...
            'r_all', 'p_all');
    end
end

%% compute median

r_median = zeros(7, 3);
r_max = zeros(7, 3);
for measure = 1:7
    for type_id = 1:3
        temp = r_all{measure, type_id};
        r_median(measure, type_id) = median(temp(~isnan(temp)));
        r_max(measure, type_id) = max(temp(~isnan(temp)));
    end
end
disp([r_median])

%% boxplot everything

diff_rate = 10;
permuteSeed_all = [1 2 3 4 5 6 7 8 9 10 13];
everything = zeros(10, 21);
for permuteSeed = permuteSeed_all
    for type_id = 1:3
        idx = find(permuteSeed_all == permuteSeed);
        load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
        for measure = 1:7
            temp = r_all{measure, type_id};
            everything(idx, measure + (type_id - 1) * 7) = median(temp(~isnan(temp)));
        end
    end
end

load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
noPermute = zeros(7, 3);
for measure = 1:7
    for type_id = 1:3
        temp = r_all{measure, type_id};
        noPermute(measure, type_id) = median(temp(~isnan(temp)));
    end
end

diff_rate = 0;
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
noDiffusion = zeros(7, 3);
for measure = 1:7
    for type_id = 1:3
        temp = r_all{measure, type_id};
        noDiffusion(measure, type_id) = median(temp(~isnan(temp)));
    end
end

boxplot(everything); grid on; hold on;
for type_id = 1:3
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noPermute(:, type_id), 'b');
end
for type_id = 1:3
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noDiffusion(:, type_id), 'r');
end

%% load results and set proper names

load('3day_allMeasures_diffRate10.mat');
r_3day = r_all; p_3day = p_all;

load('6week_allMeasures_diffRate10');
r_6w = r_all; p_62 = p_all;
clear r_all p_all;

%% range for measure 1, 2, 4

range{1} = [11:12, 21:22, 31:32]; % temporal 200-400, numClust 2-4, for measure 167, PASSED

range{2} = [101:104, 111:114, 201:204, 211:214, 301:304, 311:314]; % PASSED

range{3} = [40:60]; % PASSED

% range{4} = [101:120, 201:220, 301:320];
range{4} = [101:102, 111:112, 201:202, 211:212, 301:302, 311:312]; % PASSED

%% test these ranges on measure 1,6,7

% figure (1)
% r_6w{7,3}(range1)


subplot(341), boxplot(r_3day{1,1}(range1)); grid on;
subplot(342), boxplot(r_6w{1,1}(range1)); grid on;
subplot(343), boxplot(r_6w{1,2}(range1)); grid on;
subplot(344), boxplot(r_6w{1,3}(range1)); grid on;

subplot(345), boxplot(r_3day{1,1}(range1)); grid on;
subplot(346), boxplot(r_6w{1,1}(range1)); grid on;
subplot(347), boxplot(r_6w{1,2}(range1)); grid on;
subplot(348), boxplot(r_6w{1,3}(range1)); grid on;

subplot(349), boxplot(r_3day{1,1}(range1)); grid on;
subplot(3,4,10), boxplot(r_6w{1,1}(range1)); grid on;
subplot(3,4,11), boxplot(r_6w{1,2}(range1)); grid on;
subplot(3,4,12), boxplot(r_6w{1,3}(range1)); grid on;

%% imagesc on specific measures to tune range

measure = 4;

imagesc([r_3day{measure}(range{measure}), r_6w{measure,1}(range{measure}),...
    r_6w{measure,2}(range{measure}), r_6w{measure,3}(range{measure})]);

%% boxplot to make sure working

measure = 4;

boxplot([r_3day{measure}(range{measure}), r_6w{measure,1}(range{measure}),...
    r_6w{measure,2}(range{measure}), r_6w{measure,3}(range{measure})]);




