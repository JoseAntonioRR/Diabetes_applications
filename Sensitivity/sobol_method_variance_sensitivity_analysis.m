% Clear workspace
close all;
clear;
clc;

% Set random number generator
rng 4614;

% Number of repetitions
Rep = 10000;

% Generate Sobol sequences
p1 = sobolset(7, 'Skip', 2e4);
R1 = net(p1, Rep);
p2 = sobolset(7, 'Skip', 50);
R2 = net(p2, Rep);

% Initialize parameter vector
CCx = zeros(Rep, 7);

% Set fixed parameter values
fixed_params = [8.08764074487497e-06; 0.179281376301474; 0.0117079616864615; 427.024266568773; 0.000163848513250089; 1.60469308927761; 0.00677482687633872];

% Populate parameter matrix with fixed values
for gg = 1:Rep
    CCx(gg, :) = fixed_params;
end

% Perturb the parameters using Sobol sequences
V1 = R1 * 3 .* CCx / 2 + CCx / 2;
V2 = R2 * 3 .* CCx / 2 + CCx / 2;

% Generate normal distributions for parameter values
pd = makedist('Normal', 'mu', 52, 'sigma', 25);
t = truncate(pd, 1, 100);

% Sample from the distribution for x1 and x2
x1 = icdf(t, R1(:, 4));
x2 = icdf(t, R2(:, 4));

% Construct matrices A and B
A = [V1(:, 1:3), x1, V1(:, 5:7)];
B = [V2(:, 1:3), x2, V2(:, 5:7)];

% Initialize 3D matrix Ab for interpolation
Ab = zeros(size(A, 1), size(A, 2), size(A, 2));

% Fill Ab matrix
Ab(:, 2:7, 1) = A(:, 2:7);
Ab(:, 1, 1) = B(:, 1);

for ii = 2:1:size(Ab, 3)
    Ab(:, 1:ii-1, ii) = A(:, 1:ii-1);
    Ab(:, ii+1:end, ii) = A(:, ii+1:end);
    Ab(:, ii, ii) = B(:, ii);
end

% Set initial conditions for ODE solver
x = [8.08764074487497e-06; 0.179281376301474; 0.0117079616864615; 427.024266568773; 0.000163848513250089; 1.60469308927761; 0.00677482687633872];
Gm = 180;

D = [-x(1)*Gm-(x(2)+x(3)+1/x(4)) -x(1)*Gm -x(1)*Gm;...
    -x(5)*Gm -x(5)*Gm-(x(6)+x(7)+1/x(4)) -x(5)*Gm;...
    x(3) x(7) -1/x(4)];

M = [-x(1)*Gm*100; -x(5)*Gm*100; 0];
y0(:, 1) = linsolve(D, M);

% Define time span for ODE solver
tspan = [-1000 24];
xq = 0:0.1:tspan(2);

% Initialize matrices for interpolation results
fa_u = zeros(Rep, length(xq));
fb_u = zeros(Rep, length(xq));
fa_v = zeros(Rep, length(xq));
fb_v = zeros(Rep, length(xq));
fa_w = zeros(Rep, length(xq));
fb_w = zeros(Rep, length(xq));
fab_u = zeros(Rep, length(xq), size(Ab, 3));
fab_v = zeros(Rep, length(xq), size(Ab, 3));
fab_w = zeros(Rep, length(xq), size(Ab, 3));

% Perform ODE integration and interpolation
for jj = 1:size(A, 1)
    k1 = A(jj, :);
    [t1, y1] = ode45(@(t, y) vdp1(t, y, k1), tspan, y0);
    k2 = B(jj, :);
    [t2, y2] = ode45(@(t, y) vdp1(t, y, k2), tspan, y0);
    
    fa_u(jj, :) = interp1(t1, y1(:, 1), xq);
    fb_u(jj, :) = interp1(t2, y2(:, 1), xq);
    fa_v(jj, :) = interp1(t1, y1(:, 2), xq);
    fb_v(jj, :) = interp1(t2, y2(:, 2), xq);
    fa_w(jj, :) = interp1(t1, y1(:, 3), xq);
    fb_w(jj, :) = interp1(t2, y2(:, 3), xq);
    
    for qq = 1:size(A, 2)
        k3 = Ab(jj, :, qq);
        [t3, y3] = ode45(@(t, y) vdp1(t, y, k3), tspan, y0);
        fab_u(jj, :, qq) = interp1(t3, y3(:, 1), xq);
        fab_v(jj, :, qq) = interp1(t3, y3(:, 2), xq);
        fab_w(jj, :, qq) = interp1(t3, y3(:, 3), xq);
    end
end

% Compute sensitivity ratios
srd = zeros(7, length(xq));
for i = 1:7
    srd(i, :) = mean(fb_u .* (fab_u(:,:,i) - fa_u));
end

sqd = zeros(7, length(xq));
for i = 1:7
    sqd(i, :) = mean(fb_v .* (fab_v(:,:,i) - fa_v));
end

slld = zeros(7, length(xq));
for i = 1:7
    slld(i, :) = mean(fb_w .* (fab_w(:,:,i) - fa_w));
end

% Calculate variance for each variable
vY = var(fa_u, 0, 1);
vYY = var(fa_v, 0, 1);
vYYY = var(fa_w, 0, 1);

% Save values in .mat

save('Sobol_sensitivity.mat')