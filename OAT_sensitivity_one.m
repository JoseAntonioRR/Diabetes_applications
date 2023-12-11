%% INITIALIZE CODE
clear
clc
close all

% Points in which the solution is calculated
n = 200;
%   Common variables used
load('common_var.mat')
global SET_U;
SET_U = U;

% Save labile hemoglobin results
global oat_values
T = readtable('relevan_patients.xlsx');
Mt = T{1:60, 15:37};
X = Mt;

oat_values = zeros(size(X, 1), 5);

%% Code for oat sensivility

% Mean as initial point

% Minimized with fmincon 

% Introduction of glucose profiles
options = optimset('Display', 'off', 'Tolfun', 1e-12, 'Algorithm', 'active-set');

kinetic_var = zeros(600, 11);

% Minimization options
A = [];
B = [];
Aeq = [];
Beq = [];
beq = [];
jj = 1;
count_ij = 1;

% Get minimized results for One At a Time sensitivity
% For each variable
while jj <= size(SET_U, 2)
    ii = 1;   
    % For each patient
    while ii <= size(X, 1)
        lb = [5 * 10^(-6), 5 * 10^(-5), 0.1, 0.03133, 0.1 * 0.0025, 0.1 * 0.0025, 1/(90 * 24), 0, -0.2, -0.2];
        ub = [10^(-4), 10^(-3), 0.5, 3, 10 * 0.0025, 10 * 0.0025, 1/(20 * 24), 0.95, 0.2, 0.2]; 

        [U, fval] = fmincon(@(x)myfun(x, n, X, ii, jj), x0, A, B, Aeq, beq, lb, ub, [], options);
        
        % Save each variable and error
        kinetic_var(count_ij, :) = [U, fval];
        kinetic_var(count_ij, jj) = SET_U(jj);
        ii = ii + 1;
        count_ij = count_ij + 1;
    end
    jj = jj + 1;
end

ii = 1;

% Plot obtained error
for mnf = 1:60:541
    hold on
    plot((mnf/60) + 1, mean(kinetic_var(mnf:mnf+59, 11)), '*')
    hold off
end

save('decision.mat')

%% Runge kutta solver that calcualtes the difference with empirical values
function fobm = myfun(x, n, X, ii, lm)
    global oat_values
    global SET_U;

    cc = spline(0:0.5:2, X(ii, 1:5:21));
    x(lm) = SET_U(lm);

    % Hemoglobin function value
    tsl = zeros(5, 1);
    tsp = tsl;

    l = 2 / n;

    % Both labile types
    hlabk = x(8) * X(ii, 3);
    hlabq = (1 - x(8)) * X(ii, 3) + x(9);

    % Glycated hemoglobin values

    hglic = X(ii, 2) + x(10);


    % Initial values
    tsl(1) = hlabk + hlabq - X(ii, 3);
    tsp(1) = hglic - X(ii, 2);
    oat_values(ii, 1) = hlabk + hlabq;

    % Runge-Kutta 4th Order method
    for im = 0:3
        for i = 2:n/4
    
    kl1=l*(x(1)*((ppval(cc,(i-1+n*im/4)*l))*(100-hlabq-hlabk-hglic))-(x(3)+x(5)+x(7))*hlabk);
    
    ql1=l*(x(2)*((ppval(cc,(i-1+n*im/4)*l))*(100-hlabq-hlabk-hglic))-(x(4)+x(6)+x(7))*hlabq);
    
    hgl1=x(5)*hlabk+x(6)*hlabq-x(7)*hglic;
%-----
    
    kl2=l*(x(1)*((ppval(cc,(i-1+n*im/4)*l+0.5))*(100-(hlabk+0.5*kl1)-(hlabq+0.5*ql1)-hglic-0.5*hgl1))-...
    (x(3)+x(5)+x(7))*(hlabk+0.5*kl1));

    ql2=l*(x(2)*((ppval(cc,(i-1+n*im/4)*l+0.5))*(100-(hlabk+0.5*kl1)-(hlabq+0.5*ql1)-hglic-0.5*hgl1))-...
    (x(4)+x(6)+x(7))*(hlabq+0.5*ql1));

    hgl2=x(5)*(hlabk+0.5*kl1)+x(6)*(hlabq+0.5*ql1)-x(7)*(hglic+0.5*hgl1);
%-----
    kl3=l*(x(1)*((ppval(cc,(i-1+n*im/4+0.5)*l))*(100-(hlabk+kl2*0.5)-(hlabq+ql2*0.5)-hglic-0.5*hgl2))-...
    (x(3)+x(5)+x(7))*(hlabk+kl2*0.5));

    ql3=l*(x(2)*((ppval(cc,(i-1+n*im/4+0.5)*l))*(100-(hlabk+kl2*0.5)-(hlabq+ql2*0.5)-hglic-0.5*hgl2))-...
    (x(4)+x(6)+x(7))*(hlabq+ql2*0.5));

    hgl3=x(5)*(hlabk+0.5*kl2)+x(6)*(hlabq+0.5*ql2)-x(7)*(hglic+0.5*hgl2);
%-----
    kl4=l*(x(1)*((ppval(cc,(i+n*im/4)*l))*(100-(hlabk+kl3)-(hlabq+ql3)-(hglic+hgl3)))-...
    (x(3)+x(5)+x(7))*(hlabk+kl3));

    ql4=l*(x(2)*((ppval(cc,(i+n*im/4)*l))*(100-(hlabk+kl3)-(hlabq+ql3)-(hglic+hgl3)))-...
    (x(4)+x(6)+x(7))*(hlabq+ql3));

    hgl4=x(5)*(hlabk+kl3)+x(6)*(hlabq+ql3)-x(7)*(hglic+hgl3);
%-----

    

    hlabk =hlabk+kl1/6+kl2/3+kl3/3+kl4/6;
    
    hlabq =hlabq+ql1/6+ql2/3+ql3/3+ql4/6;

    hglic =hglic+hgl1/6+hgl2/3+hgl3/3+hgl4/6;
    %refreshed valeus   

end
tsl(im+2)=X(ii,3+5*(1+im))-hlabk-hlabq;
tsp(im+2)=X(ii,2+5*(1+im))-hglic;
oat_values(ii,im+2)=hlabk+hlabq;
i=i+1;

    kl1=l*(x(1)*((ppval(cc,(i-1+n*im/4)*l))*(100-hlabq-hlabk-hglic))-(x(3)+x(5)+x(7))*hlabk);
    
    ql1=l*(x(2)*((ppval(cc,(i-1+n*im/4)*l))*(100-hlabq-hlabk-hglic))-(x(4)+x(6)+x(7))*hlabq);
    
    hgl1=x(5)*hlabk+x(6)*hlabq-x(7)*hglic;
%-----
    
    kl2=l*(x(1)*((ppval(cc,(i-1+n*im/4)*l+0.5))*(100-(hlabk+0.5*kl1)-(hlabq+0.5*ql1)-hglic-0.5*hgl1))-...
    (x(3)+x(5)+x(7))*(hlabk+0.5*kl1));

    ql2=l*(x(2)*((ppval(cc,(i-1+n*im/4)*l+0.5))*(100-(hlabk+0.5*kl1)-(hlabq+0.5*ql1)-hglic-0.5*hgl1))-...
    (x(4)+x(6)+x(7))*(hlabq+0.5*ql1));

    hgl2=x(5)*(hlabk+0.5*kl1)+x(6)*(hlabq+0.5*ql1)-x(7)*(hglic+0.5*hgl1);
%-----
    kl3=l*(x(1)*((ppval(cc,(i-1+n*im/4+0.5)*l))*(100-(hlabk+kl2*0.5)-(hlabq+ql2*0.5)-hglic-0.5*hgl2))-...
    (x(3)+x(5)+x(7))*(hlabk+kl2*0.5));

    ql3=l*(x(2)*((ppval(cc,(i-1+n*im/4+0.5)*l))*(100-(hlabk+kl2*0.5)-(hlabq+ql2*0.5)-hglic-0.5*hgl2))-...
    (x(4)+x(6)+x(7))*(hlabq+ql2*0.5));

    hgl3=x(5)*(hlabk+0.5*kl2)+x(6)*(hlabq+0.5*ql2)-x(7)*(hglic+0.5*hgl2);
%-----
    kl4=l*(x(1)*((ppval(cc,(i+n*im/4)*l))*(100-(hlabk+kl3)-(hlabq+ql3)-(hglic+hgl3)))-...
    (x(3)+x(5)+x(7))*(hlabk+kl3));

    ql4=l*(x(2)*((ppval(cc,(i+n*im/4)*l))*(100-(hlabk+kl3)-(hlabq+ql3)-(hglic+hgl3)))-...
    (x(4)+x(6)+x(7))*(hlabq+ql3));

    hgl4=x(5)*(hlabk+kl3)+x(6)*(hlabq+ql3)-x(7)*(hglic+hgl3);
%-----
    hlabk =hlabk+kl1/6+kl2/3+kl3/3+kl4/6;
    
    hlabq =hlabq+ql1/6+ql2/3+ql3/3+ql4/6;
    
    hglic =hglic+hgl1/6+hgl2/3+hgl3/3+hgl4/6;

end

% Finalized R-K method

        % Total error
        sr1 = sum(tsl.^2) + sum(tsp.^2);

        % Error of the approximation of our simulation
        fobm = sr1;
 end




