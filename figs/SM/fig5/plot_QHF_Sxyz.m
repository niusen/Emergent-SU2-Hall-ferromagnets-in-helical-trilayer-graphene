clc;clear;close all;
restoredefaultpath;
pos_x0=0.05;
pos_y0=0.14;
length_x=0.25;
length_y=0.75;
%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\order');

markerlength=10;
Markerwidth=0.02;
LW=3;

waa_set=[0,0.2,0.4,0.6,0.8,0.85,0.9,0.95];

option=2;
%%%%%%%%%%%%%%%
pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);

for ww=1:length(waa_set)
    filenm="Sxyz_waa"+num2str(waa_set(ww))+"_option"+num2str(option)+".mat";
    load("data/"+filenm);
    values=sort(real(eig(Opx)));
    N=waa_set(ww);

    for cc=1:length(values)
        f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),(values(cc))*ones(1,markerlength),'r-','linewidth',LW);hold on;
    end
end

axis([-0.05,1,-8,8]);
set(gca,'fontsize',12)
xticks([0,0.2,0.4,0.6,0.8,1])
yticks([-8,-6,-4,-2,0,2,4,6,8])
xlabel('$w_{AA}/w_{AB}$','interpreter','latex','FontSize',16)
title("$S^x$",'interpreter','latex','FontSize',14)
set(gca,'TickLength',[0.03, 0.01])

%%%%%%%%%%%%%%%
pos=[pos_x0+1/3,pos_y0,length_x,length_y];
subplot('Position',pos);

for ww=1:length(waa_set)
    filenm="Sxyz_waa"+num2str(waa_set(ww))+"_option"+num2str(option)+".mat";
    load("data/"+filenm);
    values=sort(real(eig(Opy)));
    N=waa_set(ww);

    for cc=1:length(values)
        f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),(values(cc))*ones(1,markerlength),'r-','linewidth',LW);hold on;
    end
end

axis([-0.05,1,-8,8]);
set(gca,'fontsize',12)
xticks([0,0.2,0.4,0.6,0.8,1])
yticks([-8,-6,-4,-2,0,2,4,6,8])
xlabel('$w_{AA}/w_{AB}$','interpreter','latex','FontSize',16)
title("$S^y$",'interpreter','latex','FontSize',14)
set(gca,'TickLength',[0.03, 0.01])



%%%%%%%%%%%%%%%
pos=[pos_x0+2/3,pos_y0,length_x,length_y];
subplot('Position',pos);

for ww=1:length(waa_set)
    filenm="Sxyz_waa"+num2str(waa_set(ww))+"_option"+num2str(option)+".mat";
    load("data/"+filenm);
    values=sort(real(eig(Opz)));
    N=waa_set(ww);

    for cc=1:length(values)
        f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),(values(cc))*ones(1,markerlength),'r-','linewidth',LW);hold on;
    end
end

axis([-0.05,1,-8,8]);
set(gca,'fontsize',12)
xticks([0,0.2,0.4,0.6,0.8,1])
yticks([-8,-6,-4,-2,0,2,4,6,8])
xlabel('$w_{AA}/w_{AB}$','interpreter','latex','FontSize',16)
title("$S^z$",'interpreter','latex','FontSize',14)
set(gca,'TickLength',[0.03, 0.01])




set(gcf,'Position',[100 100 800 420])











