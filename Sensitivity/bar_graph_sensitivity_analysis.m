% Clear workspace and command window
clear;
clc;

% Load sensitivity analysis variables from a MAT file
load('Sobol_sensitivity.mat')

% Close all existing figures
close all;

% Define labels for parameters
TestL = {'k_{FL1}', 'k_{LF1}', 'k_{LS1}', '\tau_{RBC}', 'k_{FL2}', 'k_{LF2}', 'k_{LS2}'};

% Process sensitivity ratios for the first set of variables
ww = 1;
QQ = srd(:, 2:end) ./ vY(2:end);
QQ = abs(QQ);

% Compute mean sensitivity ratios for every 10 data points
for nd = 1:10:size(QQ, 2)
    QQm(:, ww) = mean(QQ(:, nd:nd+9), 2);
    ww = ww + 1;
end

% Plot sensitivity ratios for the first set of variables
figure();
bar(QQm', 'stacked');

xlabel('Time (h)');
ylabel('First-order indices');
ylim([0 1]);
set(gca, 'FontSize', 24);

% Process sensitivity ratios for the second set of variables
ww = 1;
GG = sqd(:, 2:end) ./ vYY(2:end);
GG = abs(GG);

% Compute mean sensitivity ratios for every 10 data points
for nd = 1:10:size(GG, 2)
    GGm(:, ww) = mean(GG(:, nd:nd+9), 2);
    ww = ww + 1;
end

% Plot sensitivity ratios for the second set of variables
figure();
bar(GGm', 'stacked');
% Plot details
xlabel('Time (h)');
ylabel('First-order indices');
ylim([0 1]);
set(gca, 'FontSize', 24);

% Process sensitivity ratios for the third set of variables
ww = 1;
XX = slld(:, 2:end) ./ vYYY(2:end);
XX = abs(XX);

% Compute mean sensitivity ratios for every 10 data points
for nd = 1:10:size(XX, 2)
    XXm(:, ww) = mean(XX(:, nd:nd+9), 2);
    ww = ww + 1;
end

% Plot sensitivity ratios for the third set of variables
figure();
bar(XXm', 'stacked');
legend(TestL, 'Location', 'best');
xlabel('Time (h)');
ylabel('First-order indices');
ylim([0 1]);
set(gca, 'FontSize', 24);