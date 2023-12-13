%% FLOQUET STABILITY  CALCULATOR

% Clear workspace, command window, and close all figures
clear all
clc
close all

% Load data for kinetic variables
load('UdefporPen.mat')

% Rename variables for clarity
U = fU;
m = U(:,1);
i = 0;
h = zeros(3, 10201);
q = h;
e = h;
j = h;
r = 0;

% Loop to compute eigenvalues for different values of 'a' and 'b'
for a = 50:1:200
    b = 50;
    while (b + a) < 251
        i = i + 1;

        % Solve ODE for three different initial conditions
        [t1, y1] = ode45(@(x, y) labile_sim_homogeus(x, y, @(t) a - b * cos(0.25 * pi * t), m), [0 8], [1; 0; 0]);
        h(1, i) = y1(size(t1, 1), 1);
        h(2, i) = y1(size(t1, 1), 2);
        h(3, i) = y1(size(t1, 1), 3);

        [t2, y2] = ode45(@(x, y) labile_sim_homogeus(x, y, @(t) a - b * cos(0.25 * pi * t), m), [0 8], [0; 1; 0]);
        q(1, i) = y2(size(t2, 1), 1);
        q(2, i) = y2(size(t2, 1), 2);
        q(3, i) = y2(size(t2, 1), 3);

        [t3, y3] = ode45(@(x, y) labile_sim_homogeus(x, y, @(t) a - b * cos(0.25 * pi * t), m), [0 8], [0; 0; 1]);
        j(1, i) = y3(size(t3, 1), 1);
        j(2, i) = y3(size(t3, 1), 2);
        j(3, i) = y3(size(t3, 1), 3);

        % Increment 'b' for the next iteration
        b = b + 1;
        A = [h(:, i), q(:, i), j(:, i)];
        e(:, i) = eig(A);
    end
end

% Extract relevant data
fnr = i;
fa = abs(e(1, :));
fb = abs(e(2, :));
fc = abs(e(3, :));

% Plot results
figure(1)
plot(50:1:50 + i - 1, fa(1:fnr), '--')
xlabel('Constant "b"')
ylabel('Modulus of Eigenvalue 1 ')

figure(2)
plot(50:1:50 + i - 1, fb(1:fnr), '--')
xlabel('Constant "b"')
ylabel('Modulus of Eigenvalue 2 ')

figure(3)
plot(50:1:50 + i - 1, fc(1:fnr), '--')
xlabel('Constant "b"')
ylabel('Modulus of Eigenvalue 3 ')

% Display maximum and minimum eigenvalues
disp('Maximum Eigenvalue 1:')
disp(max(fa))
disp('Maximum Eigenvalue 2:')
disp(max(fb))
disp('Minimum Eigenvalue 1:')
disp(min(fa))
disp('Minimum Eigenvalue 2:')
disp(min(fb))

% Define the labile_sim function
function out = labile_sim_homogeus(t, h, G, m)
    out(1) = m(1) * G(t) * (-h(1) - h(2) - h(3)) - (m(2) + m(3) + 1 / m(4)) * h(1);
    out(2) = m(5) * G(t) * (-h(1) - h(2) - h(3)) - (m(6) + m(7) + 1 / m(4)) * h(2);
    out(3) = m(3) * h(1) + m(7) * h(2) - h(3) / m(4);
    out = [out(1); out(2); out(3)];
end





