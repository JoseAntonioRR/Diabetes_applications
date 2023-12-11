%%INITIAL COMMON KINETIC VARIABLES FOR ALL PATIENTS

clear all
clc
close all
% Points in which the solution is calculated
n=200; 

global common_par_LAbile_hem
global chelab
global che_g
T = readtable('relevan_patients.xlsx');
Mt=T{1:60,15:37};
X = Mt;

% Variable that saves labile hemoglobin results

common_par_LAbile_hem=zeros(size(X,1),5);

%%    

% Mean as initial point

% Minimized with fmincon 

% Introduction of glucose profiles

options = optimset('Display','off','Tolfun',1e-12,'Algorithm','active-set');

ii=1;
A = []; B = []; Aeq = []; Beq = []; 
beq = [];  
 
     lb = [5*10^(-6) 5*10^(-5) 0.1 0.03133 0.1*0.0025 0.1*0.0025 1/(90*24) 0 -1 -0.2];%/50;
 ub = [10^(-4) 10^(-3) 0.5 3 10*0.0025 10*0.0025 1/(20*24) 0.95 1 0.2];%*50;  
 x0 = [5.2*10^(-5) 5.2*10^(-4) 0.3133 0.3133 0.0025 0.0025 1/(53*24) 0.7 0 0];
[U,fval] = fmincon(@(x)myfun(x,n,X),x0,A,B,Aeq,beq,lb,ub,[],options);
%Percentage errors of labile glycated and total hemoglobin
rel_err_perc=100*fval/sum(X(:,2:5:23)+X(:,3:5:23),'all');
rel_err_labi=100*chelab/sum(X(:,3:5:23),'all');
rel_err_glic=100*che_g/sum(X(:,2:5:22),'all');


%save variables that will be used elsewhere
save('common_var.mat')
%% Runge kutta solver that calcualtes the difference with empirical values

function fobm = myfun(x,n,X)
global chelab
global che_g

ii=1;
fobm=0;
chelab=0;
che_g=0;
while ii<=size(X,1)
global common_par_LAbile_hem
cc = spline(0:0.5:2, X(ii, 1:5:21));
    x(lm) = SET_U(lm);

    % Hemoglobin function value
    tsl = zeros(5, 1);
    tsp = tsl;

    l = 2 / n;

    % Both labile types
    hlabk = x(8) * X(ii, 3);
    hlabq = (1 - x(8)) * X(ii, 3) + x(9);

    hglic = X(ii, 2) + x(10);

    % Glycated hemoglobin values

    % Initiated values
    tsl(1) = hlabk + hlabq - X(ii, 3);
    tsp(1) = hglic - X(ii, 2);
    common_par_LAbile_hem(ii,1)=hlabk+hlabq;
    % Runge-Kutta 4th Order method
for im=0:3
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
    %actualizamos el valor de las integrales   

end
tsl(im+2)=X(ii,3+5*(1+im))-hlabk-hlabq;
tsp(im+2)=X(ii,2+5*(1+im))-hglic;
common_par_LAbile_hem(ii,im+2)=hlabk+hlabq;
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

sr1=sum(tsl.^2)+sum(tsp.^2);
chelab=chelab+sum(tsl.^2);
che_g=che_g+sum(tsp.^2);

% Error
fobm=fobm+sr1;
ii=ii+1;
end

end

