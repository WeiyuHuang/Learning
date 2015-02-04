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
parameter.numClust = 1:10;
parameter.temporal = linspace(100, 1000, 10);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 2

clear all; close all; clc;
diff_rate = 0;
measure = 2;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.alpha = logspace(-4, 0, 20);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)),...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 3

clear all; close all; clc;
diff_rate = 0;
measure = 3;
parameter.alpha = logspace(-4, 0, 100);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 4

clear all; close all; clc;
diff_rate = 0;
measure = 4;
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
    '.mat'], 'r', 'p');

%% 3day with measure 5

clear all; close all; clc;
diff_rate = 0;
measure = 5;
parameter.alpha = linspace(0.1, 4, 100);
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_alpha', num2str(parameter.alpha(1)), 'to', num2str(parameter.alpha(end)),...
    'n', num2str(length(parameter.alpha)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 6

clear all; close all; clc;
diff_rate = 0;
measure = 6;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');

%% 3day with measure 7

clear all; close all; clc;
diff_rate = 0;
measure = 7;
parameter.numClust = 2:11;
parameter.temporal = linspace(100, 1000, 10);
parameter.r = 0.2;
[r, p] = main_correlation( '3day', -1, diff_rate, measure, parameter );

%%
save(['../Results/3day_diff_rate', num2str(diff_rate), '_measure', num2str(measure), ...
    '_temporal', num2str(parameter.temporal(1)), 'to', num2str(parameter.temporal(end)), ...
    '_numClust', num2str(parameter.numClust(1)), 'to', num2str(parameter.numClust(end)), ...
    '.mat'], 'r', 'p');


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

