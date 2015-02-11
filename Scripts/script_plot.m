% script for plotting results

%% 1 - 28 boxplot of median ranging over all parameters.
% boxes are learning permutations with diffusion, red is no diffusion
% without permutation, blue is diffused without permutation 

% First, boxplot for permutations with diffusion
%
diff_rate = 10;
permuteSeed_3day = [1 2 3 4 5 6 13 33 233 330 333];
permuteSeed_6w = [1 2 3 4 5 6 7 8 9 10 13];
everything = zeros(11, 28);
for permuteSeed = permuteSeed_3day
    idx = find(permuteSeed_3day == permuteSeed);
    load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
    for measure = 1:7
        temp = r_all{measure, 1};
        everything(idx, measure) = median(temp(~isnan(temp)));
    end
end
for permuteSeed = permuteSeed_6w
    for type_id = 1:3
        idx = find(permuteSeed_6w == permuteSeed);
        load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
        for measure = 1:7
            temp = r_all{measure, type_id};
            everything(idx, 7 + measure + (type_id - 1) * 7) = median(temp(~isnan(temp)));
        end
    end
end

% Second, blue line for diffusion without permutation
%
noPermute = zeros(7, 4);
load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = 1:7
    temp = r_all{measure, 1};
    noPermute(measure, 1) = median(temp(~isnan(temp)));
end
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = 1:7
    for type_id = 1:3
        temp = r_all{measure, type_id};
        noPermute(measure, type_id+1) = median(temp(~isnan(temp)));
    end
end

% Third, red line for no permutation, no diffusion
%
diff_rate = 0;
noDiffusion = zeros(7, 4);
load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = 1:7
    temp = r_all{measure, 1};
    noDiffusion(measure, 1) = median(temp(~isnan(temp)));
end
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = 1:7
    for type_id = 1:3
        temp = r_all{measure, type_id};
        noDiffusion(measure, type_id+1) = median(temp(~isnan(temp)));
    end
end

% Finally, plot
boxplot(everything); grid on; hold on;
for type_id = 1:4
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noPermute(:, type_id), 'b');
end
for type_id = 1:4
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noDiffusion(:, type_id), 'r');
end

%% 2 - 28 boxplot of max ranging over all parameters.
% boxes are learning permutations with diffusion, red is no diffusion
% without permutation, blue is diffused without permutation 

% First, boxplot for permutations with diffusion
%
diff_rate = 10;
permuteSeed_3day = [1 2 3 4 5 6 13 33 233 330 333];
permuteSeed_6w = [1 2 3 4 5 6 7 8 9 10 13];
everything = zeros(11, 28);
for permuteSeed = permuteSeed_3day
    idx = find(permuteSeed_3day == permuteSeed);
    load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
    for measure = 1:7
        temp = r_all{measure, 1};
        everything(idx, measure) = max(temp(~isnan(temp)));
    end
end
for permuteSeed = permuteSeed_6w
    for type_id = 1:3
        idx = find(permuteSeed_6w == permuteSeed);
        load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
        for measure = 1:7
            temp = r_all{measure, type_id};
            everything(idx, 7 + measure + (type_id - 1) * 7) = max(temp(~isnan(temp)));
        end
    end
end

% Second, blue line for diffusion without permutation
%
noPermute = zeros(7, 4);
load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = 1:7
    temp = r_all{measure, 1};
    noPermute(measure, 1) = max(temp(~isnan(temp)));
end
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = 1:7
    for type_id = 1:3
        temp = r_all{measure, type_id};
        noPermute(measure, type_id+1) = max(temp(~isnan(temp)));
    end
end

% Third, red line for no permutation, no diffusion
%
diff_rate = 0;
noDiffusion = zeros(7, 4);
load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = 1:7
    temp = r_all{measure, 1};
    noDiffusion(measure, 1) = max(temp(~isnan(temp)));
end
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = 1:7
    for type_id = 1:3
        temp = r_all{measure, type_id};
        noDiffusion(measure, type_id+1) = max(temp(~isnan(temp)));
    end
end

% Finally, plot
boxplot(everything); grid on; hold on;
for type_id = 1:4
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noPermute(:, type_id), 'b');
end
for type_id = 1:4
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noDiffusion(:, type_id), 'r');
end


%% 3 - 24 boxplot of median ranging over some parameters
% boxes are learning permutations with diffusion, red is no diffusion
% without permutation, blue is diffused without permutation 

% Zero, define parameter range
range{1} = [11:12, 21:22, 31:32]; 
range{6} = [11:12, 21:22, 31:32]; 
range{7} = [11:12, 21:22, 31:32]; 
range{2} = [101:104, 111:114, 201:204, 211:214, 301:304, 311:314];
range{3} = [40:60];
range{4} = [101:102, 111:112, 201:202, 211:212, 301:302, 311:312];

% First, boxplot for permutations with diffusion
%
diff_rate = 10;
permuteSeed_3day = [1 2 3 4 5 6 13 33 233 330 333];
permuteSeed_6w = [1 2 3 4 5 6 7 8 9 10 13];
everything = zeros(11, 28);
for permuteSeed = permuteSeed_3day
    idx = find(permuteSeed_3day == permuteSeed);
    load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
    for measure = [1:4 6 7]
        temp = r_all{measure, 1}(range{measure});
        everything(idx, measure) = median(temp(~isnan(temp)));
    end
end
for permuteSeed = permuteSeed_6w
    for type_id = 1:3
        idx = find(permuteSeed_6w == permuteSeed);
        load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
        for measure = [1:4 6 7]
            temp = r_all{measure, type_id}(range{measure});
            everything(idx, 7 + measure + (type_id - 1) * 7) = median(temp(~isnan(temp)));
        end
    end
end

% Second, blue line for diffusion without permutation
%
noPermute = zeros(7, 4);
load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = [1:4 6 7]
    temp = r_all{measure, 1}(range{measure});
    noPermute(measure, 1) = median(temp(~isnan(temp)));
end
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = [1:4 6 7]
    for type_id = 1:3
        temp = r_all{measure, type_id}(range{measure});
        noPermute(measure, type_id+1) = median(temp(~isnan(temp)));
    end
end

% Third, red line for no permutation, no diffusion
%
diff_rate = 0;
noDiffusion = zeros(7, 4);
load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = [1:4 6 7]
    temp = r_all{measure, 1}(range{measure});
    noDiffusion(measure, 1) = median(temp(~isnan(temp)));
end
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = [1:4 6 7]
    for type_id = 1:3
        temp = r_all{measure, type_id}(range{measure});
        noDiffusion(measure, type_id+1) = median(temp(~isnan(temp)));
    end
end

% Finally, plot
boxplot(everything); grid on; hold on;
for type_id = 1:4
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noPermute(:, type_id), 'b');
end
for type_id = 1:4
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noDiffusion(:, type_id), 'r');
end



%% 4 - 24 boxplot of max ranging over some parameters
% boxes are learning permutations with diffusion, red is no diffusion
% without permutation, blue is diffused without permutation 

% Zero, define parameter range
range{1} = [11:12, 21:22, 31:32]; 
range{6} = [11:12, 21:22, 31:32]; 
range{7} = [11:12, 21:22, 31:32]; 
range{2} = [101:104, 111:114, 201:204, 211:214, 301:304, 311:314];
range{3} = [40:60];
range{4} = [101:102, 111:112, 201:202, 211:212, 301:302, 311:312];

% First, boxplot for permutations with diffusion
%
diff_rate = 10;
permuteSeed_3day = [1 2 3 4 5 6 13 33 233 330 333];
permuteSeed_6w = [1 2 3 4 5 6 7 8 9 10 13];
everything = zeros(11, 28);
for permuteSeed = permuteSeed_3day
    idx = find(permuteSeed_3day == permuteSeed);
    load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
    for measure = [1:4 6 7]
        temp = r_all{measure, 1}(range{measure});
        everything(idx, measure) = max(temp(~isnan(temp)));
    end
end
for permuteSeed = permuteSeed_6w
    for type_id = 1:3
        idx = find(permuteSeed_6w == permuteSeed);
        load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '_permuteSeed', num2str(permuteSeed), '.mat'])
        for measure = [1:4 6 7]
            temp = r_all{measure, type_id}(range{measure});
            everything(idx, 7 + measure + (type_id - 1) * 7) = max(temp(~isnan(temp)));
        end
    end
end

% Second, blue line for diffusion without permutation
%
noPermute = zeros(7, 4);
load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = [1:4 6 7]
    temp = r_all{measure, 1}(range{measure});
    noPermute(measure, 1) = max(temp(~isnan(temp)));
end
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = [1:4 6 7]
    for type_id = 1:3
        temp = r_all{measure, type_id}(range{measure});
        noPermute(measure, type_id+1) = max(temp(~isnan(temp)));
    end
end

% Third, red line for no permutation, no diffusion
%
diff_rate = 0;
noDiffusion = zeros(7, 4);
load(['../Results/3day_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = [1:4 6 7]
    temp = r_all{measure, 1}(range{measure});
    noDiffusion(measure, 1) = max(temp(~isnan(temp)));
end
load(['../Results/6week_allMeasures_diffRate', num2str(diff_rate), '.mat']);
for measure = [1:4 6 7]
    for type_id = 1:3
        temp = r_all{measure, type_id}(range{measure});
        noDiffusion(measure, type_id+1) = max(temp(~isnan(temp)));
    end
end

% Finally, plot
boxplot(everything); grid on; hold on;
for type_id = 1:4
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noPermute(:, type_id), 'b');
end
for type_id = 1:4
    plot((type_id - 1)*7+1:(type_id - 1)*7+7, noDiffusion(:, type_id), 'r');
end

