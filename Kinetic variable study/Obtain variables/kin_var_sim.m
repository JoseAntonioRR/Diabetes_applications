function fobm = kin_var_sim(x,n,X,ii,jj)
global variab_rev
global SET_U;
y=zeros(1,10);
for num_it=1:length(jj)
    y(jj(num_it,2))=jj(num_it,1);
end
for num_it=1:length(y)
    count_x=1;
    if y(num_it)==0
            y(num_it)=x(count_x); 
            count_x=count_x+1;
    end
end
cc = spline(0:0.5:2,X(ii,1:5:21));
%valor adimensional de la función hemoglobina
tsl=[0;0;0;0;0];
tsp=tsl;
l=2/n;
%valore hemoglobina labil
hlabk=jj(3)*X(ii,3);
hlabq=(1-jj(3))*X(ii,3)+jj(2);

hglic=X(ii,2)+jj(1);
%valor hemoglobina glicada

%incializamos el valore de las integrales
tsl(1)=hlabk+hlabq-X(ii,3);
tsp(1)=hglic-X(ii,2);
variab_rev(ii,1)=hlabk+hlabq;
%método Runge-Kutta Orden 4
for im=0:3
for i = 2:n/4
    
    kl1=l*(x(1)*((ppval(cc,(i-1+n*im/4)*l))*(100-hlabq-hlabk-hglic))-(jj(5)+x(3)+x(4))*hlabk);
    
    ql1=l*(x(2)*((ppval(cc,(i-1+n*im/4)*l))*(100-hlabq-hlabk-hglic))-(jj(4)+jj(6)+x(4))*hlabq);
    
    hgl1=x(3)*hlabk+jj(6)*hlabq-x(4)*hglic;
%-----
    
    kl2=l*(x(1)*((ppval(cc,(i-1+n*im/4)*l+0.5))*(100-(hlabk+0.5*kl1)-(hlabq+0.5*ql1)-hglic-0.5*hgl1))-...
    (jj(5)+x(3)+x(4))*(hlabk+0.5*kl1));

    ql2=l*(x(2)*((ppval(cc,(i-1+n*im/4)*l+0.5))*(100-(hlabk+0.5*kl1)-(hlabq+0.5*ql1)-hglic-0.5*hgl1))-...
    (jj(4)+jj(6)+x(4))*(hlabq+0.5*ql1));

    hgl2=x(3)*(hlabk+0.5*kl1)+jj(6)*(hlabq+0.5*ql1)-x(4)*(hglic+0.5*hgl1);
%-----
    kl3=l*(x(1)*((ppval(cc,(i-1+n*im/4+0.5)*l))*(100-(hlabk+kl2*0.5)-(hlabq+ql2*0.5)-hglic-0.5*hgl2))-...
    (jj(5)+x(3)+x(4))*(hlabk+kl2*0.5));

    ql3=l*(x(2)*((ppval(cc,(i-1+n*im/4+0.5)*l))*(100-(hlabk+kl2*0.5)-(hlabq+ql2*0.5)-hglic-0.5*hgl2))-...
    (jj(4)+jj(6)+x(4))*(hlabq+ql2*0.5));

    hgl3=x(3)*(hlabk+0.5*kl2)+jj(6)*(hlabq+0.5*ql2)-x(4)*(hglic+0.5*hgl2);
%-----
    kl4=l*(x(1)*((ppval(cc,(i+n*im/4)*l))*(100-(hlabk+kl3)-(hlabq+ql3)-(hglic+hgl3)))-...
    (jj(5)+x(3)+x(4))*(hlabk+kl3));

    ql4=l*(x(2)*((ppval(cc,(i+n*im/4)*l))*(100-(hlabk+kl3)-(hlabq+ql3)-(hglic+hgl3)))-...
    (jj(4)+jj(6)+x(4))*(hlabq+ql3));

    hgl4=x(3)*(hlabk+kl3)+jj(6)*(hlabq+ql3)-x(4)*(hglic+hgl3);
%-----

    

    hlabk =hlabk+kl1/6+kl2/3+kl3/3+kl4/6;
    
    hlabq =hlabq+ql1/6+ql2/3+ql3/3+ql4/6;

    hglic =hglic+hgl1/6+hgl2/3+hgl3/3+hgl4/6;
    %actualizamos el valor de las integrales   

end
tsl(im+2)=X(ii,3+5*(1+im))-hlabk-hlabq;
tsp(im+2)=X(ii,2+5*(1+im))-hglic;
variab_rev(ii,im+2)=hlabk+hlabq;
i=i+1;

    kl1=l*(x(1)*((ppval(cc,(i-1+n*im/4)*l))*(100-hlabq-hlabk-hglic))-(jj(5)+x(3)+x(4))*hlabk);
    
    ql1=l*(x(2)*((ppval(cc,(i-1+n*im/4)*l))*(100-hlabq-hlabk-hglic))-(jj(4)+jj(6)+x(4))*hlabq);
    
    hgl1=x(3)*hlabk+jj(6)*hlabq-x(4)*hglic;
%-----
    
    kl2=l*(x(1)*((ppval(cc,(i-1+n*im/4)*l+0.5))*(100-(hlabk+0.5*kl1)-(hlabq+0.5*ql1)-hglic-0.5*hgl1))-...
    (jj(5)+x(3)+x(4))*(hlabk+0.5*kl1));

    ql2=l*(x(2)*((ppval(cc,(i-1+n*im/4)*l+0.5))*(100-(hlabk+0.5*kl1)-(hlabq+0.5*ql1)-hglic-0.5*hgl1))-...
    (jj(4)+jj(6)+x(4))*(hlabq+0.5*ql1));

    hgl2=x(3)*(hlabk+0.5*kl1)+jj(6)*(hlabq+0.5*ql1)-x(4)*(hglic+0.5*hgl1);
%-----
    kl3=l*(x(1)*((ppval(cc,(i-1+n*im/4+0.5)*l))*(100-(hlabk+kl2*0.5)-(hlabq+ql2*0.5)-hglic-0.5*hgl2))-...
    (jj(5)+x(3)+x(4))*(hlabk+kl2*0.5));

    ql3=l*(x(2)*((ppval(cc,(i-1+n*im/4+0.5)*l))*(100-(hlabk+kl2*0.5)-(hlabq+ql2*0.5)-hglic-0.5*hgl2))-...
    (jj(4)+jj(6)+x(4))*(hlabq+ql2*0.5));

    hgl3=x(3)*(hlabk+0.5*kl2)+jj(6)*(hlabq+0.5*ql2)-x(4)*(hglic+0.5*hgl2);
%-----
    kl4=l*(x(1)*((ppval(cc,(i+n*im/4)*l))*(100-(hlabk+kl3)-(hlabq+ql3)-(hglic+hgl3)))-...
    (jj(5)+x(3)+x(4))*(hlabk+kl3));

    ql4=l*(x(2)*((ppval(cc,(i+n*im/4)*l))*(100-(hlabk+kl3)-(hlabq+ql3)-(hglic+hgl3)))-...
    (jj(4)+jj(6)+x(4))*(hlabq+ql3));

    hgl4=x(3)*(hlabk+kl3)+jj(6)*(hlabq+ql3)-x(4)*(hglic+hgl3);
%-----
    hlabk =hlabk+kl1/6+kl2/3+kl3/3+kl4/6;
    
    hlabq =hlabq+ql1/6+ql2/3+ql3/3+ql4/6;
    
    hglic =hglic+hgl1/6+hgl2/3+hgl3/3+hgl4/6;

end

%finalizamos el método R-K

% y el método del trapecio
%disp(tsl)
%calculamos nuestros valores ideales
sr1=sum(tsl.^2)+sum(tsp.^2);
%disp('sr1')
%disp(sr1)
%disp('sr2')
%disp(sr2)
%fprintf('hlab es %1.7f mm\n',Q1)
%fprintf('hglic es %1.7f mm\n',Q2)

%obtenemos el error relativo de la aproximación de nuestra simulación
fobm=sr1;
end
