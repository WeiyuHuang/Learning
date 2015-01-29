function [ f ] = scanTrain_expFit( MT_data, subject_id, module_id, scan_id, plot_figure)
%Exponential fit for subject_id (1~20), module_id (1~6), on scan_id(1~4)
%using only scan training data, plot exp fit result, if plot_figure == 1

% only for subject and module of interest
%
temp = MT_data{subject_id, module_id};

% only consider scan tranining data for that day
%
idx = (temp(:,2) == 99) & (temp(:,3) == scan_id);
temp = temp(idx, :);

% define x
%
x = temp(:,1);
x = x - min(x) - 1; % set start as 0.

% define y
%
y = temp(:, 4);

% exponential fit
%
f = fit(x,y,'exp1');
if (plot_figure)
    figure; plot(f,x,y)
end

end

