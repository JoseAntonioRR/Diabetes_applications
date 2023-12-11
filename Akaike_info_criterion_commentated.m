% Clear workspace and close all figures
clear;
close all;

% Load relevant data files
load("three_var_sqsg.mat");
load("four_var_sqsg.mat");
load("six_var_sqsg.mat");
load("ten_var_sqsg.mat");
load("all_common.mat");
load("two_var_sqsg.mat");

%% Calculate Akaike Information Criterion (AIC)

% Initialize vectors to store AIC values
v = zeros(1, 3);
w = zeros(1, 3);

% Calculate AIC values for different variable sets
v(1) = 60 * 10 * log(three_variable_sqsigma) + 2 * (3 * 60 + 7);
v(2) = 60 * 10 * log(four_variable_sqsigma) + 2 * (4 * 60 + 6);
v(3) = 60 * 10 * log(six_variable_sqsigma) + 2 * (6 * 60 + 4);
v(4) = 60 * 10 * log(ten_variable_sqsigma) + 2 * 10 * 60;
v(5) = 60 * 10 * log(fval) + 2 * 10;
v(6) = 60 * 10 * log(dab) + 2 * (2 * 60 + 8);

% Calculate AIC values with penalty term (P)
w = v;

% Define penalty term P for each variable set
P = zeros(1, 6);
P(1) = 2 * (3 * 60 + 7) * ((3 * 60 + 7) + 1) / (60 * 10 - (3 * 60 + 7) - 1);
P(2) = 2 * (4 * 60 + 6) * ((4 * 60 + 6) + 1) / (60 * 10 - (4 * 60 + 6) - 1);
P(3) = 2 * (6 * 60 + 4) * ((6 * 60 + 4) + 1) / (60 * 10 - (6 * 60 + 4) - 1);
P(4) = 2 * (10 * 60) * ((10 * 60) + 1) / (60 * 10 - (10 * 60) - 1);
P(5) = 2 * (10) * ((10) + 1) / (60 * 10 - 10 - 1);
P(6) = 2 * (2 * 60 + 8) * ((2 * 60 + 8) + 1) / (60 * 10 - (2 * 60 + 8) - 1);

% Add penalty term to AIC values
w = w + P;

%% Plot Results

% Plot AIC values with and without penalty
figure();
box on;
hold on;
w(4) = abs(w(4)); % Ensure non-negative value for better visualization
h = [5, 6, 1:4];
hm = [5, 6, 1:3];
plot(1:6, v(h), '-*', 1:5, w(hm), '-d', 'LineWidth', 1.8);

% Add labels, legend, and save the figure
legend('Standard', 'Low sample', 'Location', 'best');
ylabel('Akaike Information Criterion');
xticks([1 2 3 4 5 6]);
xticklabels({'Variables = 1', 'Variables = 2', 'Variables = 3', 'Variables = 4', 'Variables = 6', 'Variables = 10'});
hold off;
saveas(gcf, 'akaike_graph.pdf');
h = [5, 1:4];
saveas(gcf, 'akaike_graph_complete.pdf');

% Plot AIC values with penalty only
figure();
box on;
hold on;
plot(1:5, v(h), '-*', 'LineWidth', 1.8);

% Add labels, legend, and save the figure
legend('Standard', 'Low sample', 'Location', 'best');
ylabel('Akaike Information Criterion');
xticks([1 2 3 4 5]);
xticklabels({'Variables = 1', 'Variables = 3', 'Variables = 4', 'Variables = 6', 'Variables = 10'});
hold off;
saveas(gcf, 'akaike_graph.pdf');