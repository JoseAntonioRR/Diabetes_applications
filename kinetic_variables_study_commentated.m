% Set the number of iterations
n = 200;

% Load common variables
load('common_var.mat');

% Define global variable SET_U
global SET_U;
SET_U = U;

% Initialize variable for storing relevant indices
jj = zeros(1, 4);

% Define global variable for storing relevant variables
global variab_rev;
variab_rev = zeros(size(X, 1), 5);

% Read relevant patient data from Excel file
T = readtable('relevan_patients.xlsx');
Mt = T{1:60, 15:37};
X = Mt;

%% Optimization Setup

% Set options for optimization
options = optimset('Display', 'off', 'Tolfun', 1e-12, 'Algorithm', 'active-set');

% Initialize matrices for optimization constraints
A = [];
B = [];
Aeq = [];
Beq = [];
beq = [];

% Define lower and upper bounds for optimization variables
lb = [5 * 10^(-6), 5 * 10^(-5), 0.1, 0.03133, 0.1 * 0.0025, 0.1 * 0.0025, 1 / (90 * 24), 0, -0.2, -0.2];
ub = [10^(-4), 10^(-3), 0.5, 3, 10 * 0.0025, 10 * 0.0025, 1 / (20 * 24), 0.95, 0.2, 0.2];

% Initialize iteration counter
ii = 1;

% Set up variables for optimization
vec_d = [10, 8, 7, 5, 4, 3];
jj(:, 1) = SET_U(vec_d);
jj(:, 2) = vec_d;

% Remove variables from lb, ub, and x0 that correspond to jj
for count_it = length(vec_d):-1:1
    lb(vec_d(count_it)) = [];
    ub(vec_d(count_it)) = [];
    x0(vec_d(count_it)) = [];
end

% Perform optimization for each row in X
while ii <= size(X, 1)
    % Perform optimization using fmincon
    [U, fval] = fmincon(@(x) myfun(x, n, X, ii, jj), x0, A, B, Aeq, beq, lb, ub, [], options);
    rqw(ii, :) = [U, fval];
    ii = ii + 1;
end

% Save the optimized variables and indices
save('save_variables.mat', "rqw", "jj");