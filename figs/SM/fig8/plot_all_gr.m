clc;clear;close all;
restoredefaultpath;
pos_x0=0.06;
pos_y0=0.12;
length_x=0.27;
length_y=0.78;
load("data\mymap.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley\paircorrel');
pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);
load("data\gr_36_24_waa0.6.mat");

Cb1=parula(250);
Cb1=Cb1(end:-1:1,:);

Cb2=mycmap;
Cb=[Cb1;Cb2];

value_range=[min(real(gr_values)),max(real(gr_values))];
for cq=1:length(gr_values)
    
    % cq
    % gr_values(cq)
    coord=gr_points(cq,:);
     
    Color=determine_color2(value_range,real(gr_values(cq)),Cb);
    plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',10);hold on;%colorbar;

end
axis([coord0(1),coord0(1)+Dx*norm(ax),coord0(2),coord0(2)+3*Dy*norm(ay)])

xlabel('$x/a_M$','interpreter','latex','fontsize',13)
ylabel('$y/a_M$','interpreter','latex','fontsize',13)


% yticks([coord0(2)+ norm(ay)*0.5  coord0(2)+ norm(ay)*1.5  coord0(2)+ norm(ay)*2.5  coord0(2)+ norm(ay)*3.5  ])
% yticklabels({'0','1','2','3'})
% xticks([coord0(1)  coord0(1)+ norm(ay)  coord0(1)+ norm(ay)*2  coord0(1)+ norm(ay)*3  ])
% xticklabels({'0','1','2','3'})


yticks([-2*norm(ay),  0,  2*norm(ay), 4*norm(ay), 6*norm(ay), 8*norm(ay), 10*norm(ay)])
yticklabels({'-2','0','2','4','6','8','10'})
xticks([0, 2*norm(ay)   norm(ay)*4     norm(ay)*6 ])
xticklabels({'0','2','4','6'})

title("$\nu_{\rm{total}}=3+2/3, N_s=36$",'interpreter','latex','FontSize',12)



qucolor=[0.4660 0.6740 0.1880];
qucolor=[0.9290 0.6940 0.1250];
qucolor=[0.2 0.6 0.2];
p1 = [4*norm(ax) 6*norm(ay)];                         % First Point
p2 = [6*norm(ax) 6*norm(ay)];                         % Second Point
dp = p2-p1;                         % Difference
qu1=quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'MaxHeadSize',2,'Color',qucolor);

p1 = [4*norm(ax) 6*norm(ay)];                         % First Point
p2 = [3*norm(ax) 7.5*norm(ay)];                         % Second Point
dp = p2-p1;                         % Difference
qu2=quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'MaxHeadSize',2,'Color',qucolor);


%real space cluster boundary
p1=[0,6*norm(ay)];
p2=[6*norm(ax),9*norm(ay)];
p3=[6*norm(ax),3*norm(ay)];
p4=[0,0*norm(ay)];

plot([p1(1),p2(1)],[p1(2),p2(2)],'k--','LineWidth',2);hold on;
plot([p2(1),p3(1)],[p2(2),p3(2)],'k--','LineWidth',2);hold on;
plot([p3(1),p4(1)],[p3(2),p4(2)],'k--','LineWidth',2);hold on;
plot([p4(1),p1(1)],[p4(2),p1(2)],'k--','LineWidth',2);hold on;




%set(gca,'TickDir','out');
% annotation('doublearrow', x, y)
set(gca,'TickLength',[0.03, 0.01])



caxis(value_range);
Ax=gca;
colormap(gca,Cb);
cb=colorbar;
%ylabel(cb,'Pair correlation $g(\rm{r}-\rm{r}^{\prime})$','Rotation',270,'interpreter','latex','fontsize',13)

Ax.FontSize=12;
% axis off

%text(coord0(1)-norm(ax)*1.2,coord0(2)+Dy*norm(ay)*1.1,'(d)','interpreter','latex','fontsize',14)
%title("Pair correlation  $g(\rm{r})$",'interpreter','latex','FontSize',12)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+1/3,pos_y0,length_x,length_y];
subplot('Position',pos);


load("data\gr_28_14_waa0.6.mat");

Cb1=parula(250);
Cb1=Cb1(end:-1:1,:);

Cb2=mycmap;
Cb=[Cb1;Cb2];

value_range=[min(real(gr_values)),max(real(gr_values))];
for cq=1:length(gr_values)
    
    % cq
    % gr_values(cq)
    coord=gr_points(cq,:);
     
    Color=determine_color2(value_range,real(gr_values(cq)),Cb);
    plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',10);hold on;%colorbar;

end
axis([coord0(1),coord0(1)+Dx*norm(ax),coord0(2),coord0(2)+3*Dy*norm(ay)])

xlabel('$x/a_M$','interpreter','latex','fontsize',13)
%ylabel('$y/a_M$','interpreter','latex','fontsize',13)


% yticks([coord0(2)+ norm(ay)*0.5  coord0(2)+ norm(ay)*1.5  coord0(2)+ norm(ay)*2.5  coord0(2)+ norm(ay)*3.5  ])
% yticklabels({'0','1','2','3'})
% xticks([coord0(1)  coord0(1)+ norm(ay)  coord0(1)+ norm(ay)*2  coord0(1)+ norm(ay)*3  ])
% xticklabels({'0','1','2','3'})

yticks([0,  norm(ay)*2, norm(ay)*4, norm(ay)*6, norm(ay)*8, norm(ay)*10, norm(ay)*12, norm(ay)*14, norm(ay)*16])
yticklabels({'0','2','4','6','8','10','12','14','16'})
xticks([-2*norm(ay) 0*norm(ay)   norm(ay)*2   norm(ay)*4  ])
xticklabels({'-2','0','2','4'})

title("$\nu_{\rm{total}}=3+1/2, N_s=28$",'interpreter','latex','FontSize',12)


%set(gca,'TickDir','out');
% annotation('doublearrow', x, y)
set(gca,'TickLength',[0.03, 0.01])




qucolor=[0.4660 0.6740 0.1880];
qucolor=[0.9290 0.6940 0.1250];
qucolor=[0.2 0.6 0.2];
p1 = [2*norm(ax) 1*norm(ay)];                         % First Point
p2 = [4*norm(ax) 0*norm(ay)];                         % Second Point
dp = p2-p1;                         % Difference
qu1=quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'MaxHeadSize',2,'Color',qucolor);

p1 = [2*norm(ax) 1*norm(ay)];                         % First Point
p2 = [2*norm(ax) 3*norm(ay)];                         % Second Point
dp = p2-p1;                         % Difference
qu2=quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'MaxHeadSize',2,'Color',qucolor);


%real space cluster boundary
p1=[0,0];
p2=[6*norm(ax),-1*norm(ay)];
p3=[4*norm(ax),4*norm(ay)];
p4=[-2*norm(ax),5*norm(ay)];

plot([p1(1),p2(1)],[p1(2),p2(2)],'k--','LineWidth',2);hold on;
plot([p2(1),p3(1)],[p2(2),p3(2)],'k--','LineWidth',2);hold on;
plot([p3(1),p4(1)],[p3(2),p4(2)],'k--','LineWidth',2);hold on;
plot([p4(1),p1(1)],[p4(2),p1(2)],'k--','LineWidth',2);hold on;





caxis(value_range);
Ax=gca;
colormap(gca,Cb);
cb=colorbar;
%ylabel(cb,'Pair correlation $g(\rm{r}-\rm{r}^{\prime})$','Rotation',270,'interpreter','latex','fontsize',13)

Ax.FontSize=12;
% axis off

%text(coord0(1)-norm(ax)*1.2,coord0(2)+Dy*norm(ay)*1.1,'(d)','interpreter','latex','fontsize',14)
%title("Pair correlation  $g(\rm{r})$",'interpreter','latex','FontSize',12)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+2/3,pos_y0,length_x,length_y];
subplot('Position',pos);

load("data\gr_36_12_waa0.6.mat");

Cb1=parula(250);
Cb1=Cb1(end:-1:1,:);

Cb2=mycmap;
Cb=[Cb1;Cb2];

value_range=[min(real(gr_values)),max(real(gr_values))];
for cq=1:length(gr_values)
    
    % cq
    % gr_values(cq)
    coord=gr_points(cq,:);
     
    Color=determine_color2(value_range,real(gr_values(cq)),Cb);
    plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',10);hold on;%colorbar;

end
axis([coord0(1),coord0(1)+Dx*norm(ax),coord0(2),coord0(2)+3*Dy*norm(ay)])

xlabel('$x/a_M$','interpreter','latex','fontsize',13)
%ylabel('$y/a_M$','interpreter','latex','fontsize',13)


% yticks([coord0(2)+ norm(ay)*0.5  coord0(2)+ norm(ay)*1.5  coord0(2)+ norm(ay)*2.5  coord0(2)+ norm(ay)*3.5  ])
% yticklabels({'0','1','2','3'})
% xticks([coord0(1)  coord0(1)+ norm(ay)  coord0(1)+ norm(ay)*2  coord0(1)+ norm(ay)*3  ])
% xticklabels({'0','1','2','3'})

% yticks([-2*norm(ay),  0,  norm(ay)*2])
% yticklabels({'-2','0','2'})
yticks([-2*norm(ay),  0,  2*norm(ay), 4*norm(ay), 6*norm(ay), 8*norm(ay), 10*norm(ay)])
yticklabels({'-2','0','2','4','6','8','10'})
xticks([norm(ay)*0   norm(ay)*2   norm(ay)*4  norm(ay)*6 ])
xticklabels({'0','2','4','6'})

title("$\nu_{\rm{total}}=3+1/3, N_s=36$",'interpreter','latex','FontSize',12)



%set(gca,'TickDir','out');
% annotation('doublearrow', x, y)
set(gca,'TickLength',[0.03, 0.01])






qucolor=[0.4660 0.6740 0.1880];
qucolor=[0.9290 0.6940 0.1250];
qucolor=[0.2 0.6 0.2];
p1 = [4*norm(ax) 6*norm(ay)];                         % First Point
p2 = [5*norm(ax) 5.5*norm(ay)];                         % Second Point
dp = p2-p1;                         % Difference
qu1=quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'MaxHeadSize',2,'Color',qucolor);

p1 = [4*norm(ax) 6*norm(ay)];                         % First Point
p2 = [4*norm(ax) 7*norm(ay)];                         % Second Point
dp = p2-p1;                         % Difference
qu2=quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'MaxHeadSize',2,'Color',qucolor);


%real space cluster boundary
p1=[0,6*norm(ay)];
p2=[6*norm(ax),9*norm(ay)];
p3=[6*norm(ax),3*norm(ay)];
p4=[0,0*norm(ay)];
plot([p1(1),p2(1)],[p1(2),p2(2)],'k--','LineWidth',2);hold on;
plot([p2(1),p3(1)],[p2(2),p3(2)],'k--','LineWidth',2);hold on;
plot([p3(1),p4(1)],[p3(2),p4(2)],'k--','LineWidth',2);hold on;
plot([p4(1),p1(1)],[p4(2),p1(2)],'k--','LineWidth',2);hold on;




caxis(value_range);
Ax=gca;
colormap(gca,Cb);
cb=colorbar;
%ylabel(cb,'Pair correlation $g(\rm{r}-\rm{r}^{\prime})$','Rotation',270,'interpreter','latex','fontsize',13)

Ax.FontSize=12;
% axis off

%text(coord0(1)-norm(ax)*1.2,coord0(2)+Dy*norm(ay)*1.1,'(d)','interpreter','latex','fontsize',14)
%title("Pair correlation  $g(\rm{r})$",'interpreter','latex','FontSize',12)




set(gcf,'Position',[100 100 900 600])



%exportgraphics(gcf,'all_gr.eps','ContentType','vector')
exportgraphics(gcf,'all_gr.jpg')



function Color=determine_color2(value_range,value,cb)
assert(value_range(1)<=value);
assert(value<=value_range(2));
Values=linspace(value_range(1),value_range(2),size(cb,1));

for cc=1:size(cb,1)-1
    if (value>Values(cc))&&(value<Values(cc+1))
        Color=(cb(cc,:)+cb(cc+1,:))/2;
        break;
    elseif (value==Values(1))|(value==Values(end))
        Color=cb(cc,:);

    end
end
    

end


function Color=determine_color(value_range,value)
assert(value_range(1)<=value);
assert(value<=value_range(2));
    color1=[0,0,1];
    color2=[1,1,1];
    color3=[1,0,0];
    
    value1=value_range(1);
    value2=mean(value_range);
    value3=value_range(2);
    if value<value2
        da=value-value1;
        db=value2-value;
        Color=(da*color2+db*color1)/(da+db);
    elseif value==value2
        Color=color2;
    elseif value>value2
        da=value-value2;
        db=value3-value;
        Color=(da*color3+db*color2)/(da+db);
    end

end