%% INITIALIZE CODE

clear
clc
close all
% Points in which the solution is calculated

load('num_variable.mat')
%   Variables used

%%
n=200;

%variables used 

tem_rqw=zeros(60,10);

%variables are selected as either fixed or subject dependent

%1) subject dependant

for count_p=1:size(rqw,2)
    tem_rqw(:,jj(count_p,2))=rqw(:,count_p);
end

%2) fixed

for count_p=1:size(tem_rqw,2)
    if max(tem_rqw(:,count_p))==0
        tem_rqw(:,count_p)=SET_U(count_p);
    end
end
rqw=tem_rqw;

%Patients that will be plotted
ploted_patients=[3;4;5;6;12;13;16;20;28;33;37;40;45;54;56;60;Inf];

overallM1=0;
overallM2=0;
overallI1=100;
overallI2=100;

cnof4s=1;
cnof1s=1;

%range of graphs
load('max_and_mins.mat')
%%

T = readtable('relevan_patients.xlsx');
Mt=T{1:60,15:37};
X = Mt;


for ii=1:size(X,1)
i=1;
cc = spline(0:0.5:2,X(ii,1:5:21));
%initial values
tsl=[0;0;0;0;0];
tsp=tsl;
l=2/n;
    % Hemoglobin function value

    % Both labile types
hlabk(ii,i)=rqw(ii,8)*X(ii,3);
hlabq(ii,i)=(1-rqw(ii,8))*X(ii,3)+rqw(ii,9);

% Glycated hemoglobin values

hglic(ii,i)=X(ii,2)+rqw(ii,10);


    % Initial values

important_Hlbcgxd(ii,1)=hlabk(ii,i)+hlabq(ii,i);
    % Runge-Kutta 4th Order method
for i = 2:n
    
    kl1=l*(rqw(ii,1)*((ppval(cc,(i-1)*l))*(100-hlabq(ii,i-1)-hlabk(ii,i-1)-hglic(ii,i-1)))-(rqw(ii,3)+rqw(ii,5)+rqw(ii,7))*hlabk(ii,i-1));
    
    ql1=l*(rqw(ii,2)*((ppval(cc,(i-1)*l))*(100-hlabq(ii,i-1)-hlabk(ii,i-1)-hglic(ii,i-1)))-(rqw(ii,4)+rqw(ii,6)+rqw(ii,7))*hlabq(ii,i-1));
    
    hgl1=rqw(ii,5)*hlabk(ii,i-1)+rqw(ii,6)*hlabq(ii,i-1)-rqw(ii,7)*hglic(ii,i-1);
%-----
    
    kl2=l*(rqw(ii,1)*((ppval(cc,(i-0.5)*l+0.5))*(100-(hlabk(ii,i-1)+0.5*kl1)-(hlabq(ii,i-1)+0.5*ql1)-hglic(ii,i-1)-0.5*hgl1))-...
    (rqw(ii,3)+rqw(ii,5)+rqw(ii,7))*(hlabk(ii,i-1)+0.5*kl1));

    ql2=l*(rqw(ii,2)*((ppval(cc,(i-0.5)*l+0.5))*(100-(hlabk(ii,i-1)+0.5*kl1)-(hlabq(ii,i-1)+0.5*ql1)-hglic(ii,i-1)-0.5*hgl1))-...
    (rqw(ii,4)+rqw(ii,6)+rqw(ii,7))*(hlabq(ii,i-1)+0.5*ql1));

    hgl2=rqw(ii,5)*(hlabk(ii,i-1)+0.5*kl1)+rqw(ii,6)*(hlabq(ii,i-1)+0.5*ql1)-rqw(ii,7)*(hglic(ii,i-1)+0.5*hgl1);
%-----
    kl3=l*(rqw(ii,1)*((ppval(cc,(i-0.5)*l))*(100-(hlabk(ii,i-1)+kl2*0.5)-(hlabq(ii,i-1)+ql2*0.5)-hglic(ii,i-1)-0.5*hgl2))-...
    (rqw(ii,3)+rqw(ii,5)+rqw(ii,7))*(hlabk(ii,i-1)+kl2*0.5));

    ql3=l*(rqw(ii,2)*((ppval(cc,(i-0.5)*l))*(100-(hlabk(ii,i-1)+kl2*0.5)-(hlabq(ii,i-1)+ql2*0.5)-hglic(ii,i-1)-0.5*hgl2))-...
    (rqw(ii,4)+rqw(ii,6)+rqw(ii,7))*(hlabq(ii,i-1)+ql2*0.5));

    hgl3=rqw(ii,5)*(hlabk(ii,i-1)+0.5*kl2)+rqw(ii,6)*(hlabq(ii,i-1)+0.5*ql2)-rqw(ii,7)*(hglic(ii,i-1)+0.5*hgl2);
%-----
    kl4=l*(rqw(ii,1)*((ppval(cc,(i)*l))*(100-(hlabk(ii,i-1)+kl3)-(hlabq(ii,i-1)+ql3)-(hglic(ii,i-1)+hgl3)))-...
    (rqw(ii,3)+rqw(ii,5)+rqw(ii,7))*(hlabk(ii,i-1)+kl3));

    ql4=l*(rqw(ii,2)*((ppval(cc,(i)*l))*(100-(hlabk(ii,i-1)+kl3)-(hlabq(ii,i-1)+ql3)-(hglic(ii,i-1)+hgl3)))-...
    (rqw(ii,4)+rqw(ii,6)+rqw(ii,7))*(hlabq(ii,i-1)+ql3));

    hgl4=rqw(ii,5)*(hlabk(ii,i-1)+kl3)+rqw(ii,6)*(hlabq(ii,i-1)+ql3)-rqw(ii,7)*(hglic(ii,i-1)+hgl3);
%-----

    

    hlabk(ii,i) =hlabk(ii,i-1)+kl1/6+kl2/3+kl3/3+kl4/6;
    
    hlabq(ii,i) =hlabq(ii,i-1)+ql1/6+ql2/3+ql3/3+ql4/6;

    hglic(ii,i) =hglic(ii,i-1)+hgl1/6+hgl2/3+hgl3/3+hgl4/6;



end
end


ii=1;
%Plotting section
while ii<=size(X,1)

fig = figure;
cc = spline(0:0.5:2,X(ii,1:5:21));

t1=l:l:2;
psf=hlabk(ii,:)+hlabq(ii,:);
xconf = [t1 t1(end:-1:1)] ;         
yconf = [psf+0.1 psf(end:-1:1)-0.1];



hold on
box on
  yyaxis left
  ylabel('Glucose (mg dl^{-1})','FontSize',16)
  ylim([60 140]);
  yticks([60:20:140])
 

h(3)=plot(t1,ppval(cc,t1),'b-.');
h(4)=plot(0:0.5:2,X(ii,1:5:21),'bo','MarkerFaceColor',[ 0.5843 0.8157 0.9882],'MarkerSize',12);
  ax = gca;
  ax.FontSize = 16;
M1=max(ppval(cc,t1));
I1=min(ppval(cc,t1));
yyaxis right
ylabel('HbA_{1d} (%)','FontSize',16)
set(gca,'FontSize',12)
 p = fill(xconf,yconf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';   
alpha(.5)
  h(1)=plot(t1,hlabk(ii,:)+hlabq(ii,:),'-','color',[0.8500, 0.3250, 0.0980]);
  h(2)=plot(0:0.5:2,X(ii,3:5:23),'o','color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor','#EDB120','MarkerSize',12);
  M2=max(hlabk(ii,:)+hlabq(ii,:));
  I2=min(hlabk(ii,:)+hlabq(ii,:));

ax = gca;

xlabel('Time (h)','FontSize',16)
xlim([-0.01 2.01])
xticks([0:0.5:2])

%find neccessary range of graphs
if M1>overallM1
    overallM1=M1;
end

if M2>overallM2
    overallM2=M2;
end

if I1<overallI1
    overallI1=I1;
end

if I2<overallI2
    overallI2=I2;
end

Annotation.LegendInformation.IconDisplayStyle = 'off';
hold off
drawnow
    frame = getframe(fig);
    im{ii} = frame2im(frame);
    ii=ii+1;
 
hold on
box on
ii=ii-1;
  yyaxis left
  ylabel('Glucose (mg dl^{-1})','FontSize',16)
  overallM1=260;
  ylim([overallI1 overallM1]);
  yticks([overallI1:30:overallM1])

h(3)=plot(t1,ppval(cc,t1),'b-.');
h(4)=plot(0:0.5:2,X(ii,1:5:21),'bo','MarkerFaceColor',[ 0.5843 0.8157 0.9882],'MarkerSize',12);
M1=max(ppval(cc,t1));
yyaxis right
ylabel('HbA_{1d} (%)')
ylim([0.499 3.001]);
 yticks([0.5:0.5:3])
 p = fill(xconf,yconf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';   
alpha(.5)
  h(1)=plot(t1,hlabk(ii,:)+hlabq(ii,:),'--','color',[0.8500, 0.3250, 0.0980]);
  h(2)=plot(0:0.5:2,X(ii,3:5:23),'o','color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor','#EDB120','MarkerSize',12);
  M2=max(hlabk(ii,:)+hlabq(ii,:));


 a = get(gca,'YTickLabel');
  set(gca,'YTickLabel',a,'FontSize',20)
xlabel('Time (h)')
xlim([-0.01 2.01])
ii=ii+1;
saveas(gcf,sprintf('FIG_true_max_%d.pdf',ii-1))

end
figure()
plot(t1,hlabk(round(ii/2),:),'-o',t1,hlabq(round(ii/2),:),'--')

legend('slow','fast')


filename = 'save_as.gif';

for idx = 1:size(X,1)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

%%
close all

cn_mg=1;
ii=1;
while ii<=size(X,1)
    
% Save in single figure
if ii==ploted_patients(cn_mg)
disp(ii)
cc = spline(0:0.5:2,X(ii,1:5:21));


bb=subplot(4,4,cn_mg);

cn_mg=cn_mg+1;


t1=l:l:2;
psf=hlabk(ii,:)+hlabq(ii,:);
xconf = [t1 t1(end:-1:1)] ;         
yconf = [psf+0.1 psf(end:-1:1)-0.1];

%figure


hold on
box on
  of1s=[1;5;9;13;Inf];

  yyaxis left
if cn_mg-1==of1s(cnof1s)
    cnof1s=cnof1s+1;

  ylabel('Glucose (mg/dl)','FontSize',7)
  
  
  
end
yticks([60:40:260]);
ylim([59.99 260.01]);

h(3)=plot(t1,ppval(cc,t1),'b-');
h(4)=plot(0:0.5:2,X(ii,1:5:21),'bo','MarkerFaceColor',[ 0.5843 0.8157 0.9882],'MarkerSize',4.5);
  ax = gca;
  ax.FontSize = 16;
M1=max(ppval(cc,t1));
I1=min(ppval(cc,t1));
of4s=[4;8;12;16;Inf];
yyaxis right
if cn_mg-1==of4s(cnof4s)
    cnof4s=cnof4s+1;
ylabel('HbA_{1d} (%)','FontSize',7)
set(gca,'FontSize',12)



end
 p = fill(xconf,yconf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';   
alpha(.5)
  h(1)=dashline(t1,hlabk(ii,:)+hlabq(ii,:),2,0.5,2,0.5,'color',[0.8500, 0.3250, 0.0980]);
  h(2)=plot(0:0.5:2,X(ii,3:5:23),'o','color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor','#EDB120','MarkerSize',4.5);
  M2=max(hlabk(ii,:)+hlabq(ii,:));
  I2=min(hlabk(ii,:)+hlabq(ii,:));

  ylim([0.9 3.001]);
  yticks(1:0.5:3)
ax = gca;
xlim([-0.01 2.01])
if cn_mg>13
xlabel('Time (h)','FontSize',7)
xlim([-0.01 2.01])

end
change_1=0;
change_2=0;
qmm={'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI'};
txt = qmm(cn_mg-1);
text(0.05,2.82,txt,'FontName','Helvetica')
meangly=mean(X(ii,2:5:22));
etc_gly=X(ii,2:5:22);
S_gly=std(X(ii,2:5:22));
if cn_mg==8
    change_1=0;
    change_2=-1;
end
text_math="HbA_{1c} \in [" + round(min(etc_gly),2)  + "," + round(max(etc_gly),2) +"]";
text(0.45+change_1,2.67+change_2,text_math,'FontName','Helvetica','FontSize',6,'Interpreter','tex')
change_1=0;
change_2=0;
xticks([0:0.5:2])
 a = get(gca,'YTickLabel');
  set(gca,'YTickLabel',a,'FontSize',7)
if M1>overallM1
    overallM1=M1;
end

if M2>overallM2
    overallM2=M2;
end

if I1<overallI1
    overallI1=I1;
end

if I2<overallI2
    overallI2=I2;
end

Annotation.LegendInformation.IconDisplayStyle = 'off';

hold off

    

hold on
box on
[mod(cn_mg-1,5)/4,floor((cn_mg-1)/5)/4];

end
ii=ii+1;


end


saveas(gcf,'single_graph_mosaic.pdf')
