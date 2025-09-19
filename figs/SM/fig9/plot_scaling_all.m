clc;clear;close all;
restoredefaultpath;
pos_x0=0.07;
pos_y0=0.12;
length_x=0.25;
length_y=0.32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+2/3,pos_y0+0.5,length_x,length_y];
subplot('Position',pos);

Ees=[];
Nset=[];


%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley')
load("data/Es_24A_8_waa0.6.mat");
Ns=24;
plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(2)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(3)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(4)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(4)-E_all(1)];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27_HF_2valley')
load("data/Es_27_9_waa0.6.mat");
Ns=27;
plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(2)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(3)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(4)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(4)-E_all(1)];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N30A_HF_2valley')
load("data/Es_30A_10_waa0.6.mat");
Ns=30;
plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(2)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(3)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(4)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(4)-E_all(1)];
Nset=[Nset,Ns];


%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley')
load("data/Es_36_12_waa0.6.mat");
Ns=36;
plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(2)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(3)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(4)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(4)-E_all(1)];
Nset=[Nset,Ns];

range=2:4;
xs=1./Nset(range);
ys=Ees(range);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'b--');hold on;



axis([0,1/20,0,2]);
set(gca,'fontsize',12)
% xticks([0 1/36 1/30 1/27  1/24 ])
% xticklabels({'0','1/36','1/30', '1/27','1/24'})
xticks([1/36 1/30  1/27 1/24 ])
xticklabels({36','30','27','24'})
ax = gca; 
ax.XAxis.FontSize=9;
% xticks([0 1/36  1/24 ])
% xticklabels({'0','1/36','1/24'})
set(gca,'TickLength',[0.03, 0.01])
text(1/55,1.7,'Lowest excitations',  'Interpreter', 'latex', 'fontsize', 12,'color','b');
xlabel('$1/N_s$','Interpreter','latex','FontSize',12);
ylabel('$E_n-E_0$ (meV)','interpreter','latex','fontsize',12)
title("$\nu_{\rm{total}}=3+1/3$, FCI",'interpreter','latex','fontsize',13,'Color','r');

text(-0.01,2.2,'(c)','interpreter','latex','fontsize',14)

% text(0.024,0.4,'$N_g=\,\,\,\,3 \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,3 \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,3$','interpreter','latex','fontsize',12,'color','r')
% text(0.0215,0.3,'$3$ ground states','interpreter','latex','fontsize',12,'color','r')
text(0.014,0.5,'$N_g=3$','interpreter','latex','fontsize',12,'color','r')

text(0.014,-0.23,"$N_s=$",'interpreter','latex','fontsize',12,'color','k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+1/3,pos_y0+0.5,length_x,length_y];
subplot('Position',pos);

Ees=[];
Nset=[];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N16_HF_2valley')
load("data/Es_16_8_waa0.6.mat");
Ns=16;
%3,2,2,2
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(9)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:9
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(10)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(10)-E_all(1)];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N20B_HF_2valley')
load("data/Es_20_10_waa0.6.mat");
Ns=20;
%3,2,2,2
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(9)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:11
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(12)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(12)-E_all(1)];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley')
load("data/Es_24A_12_waa0.6.mat");
Ns=24;
%3,3,3,4
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(13)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:13
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(14)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(14)-E_all(1)];
Nset=[Nset,Ns];


%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley')
load("data/Es_28_14_waa0.6.mat");
Ns=28;
%3,4,4,4
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(15)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:15
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(16)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(16)-E_all(1)];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N32_HF_2valley')
load("data/Es_32_16_waa0.6.mat");
Ns=32;
%5,4,4,4
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(17)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:17
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(18)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(18)-E_all(1)];
Nset=[Nset,Ns];


range=2:5;
xs=1./Nset(range);
ys=Ees(range);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'b--');hold on;

%axis([0.0,1/16,0,8]);
axis([0,1/15,0,8.2]);
set(gca,'fontsize',12)
% xticks([0 1/36 1/30 1/27  1/24 ])
% xticklabels({'0','1/36','1/30', '1/27','1/24'})
yticks([0 4 8])
xticks([1/32  1/28  1/24 1/20 1/16])
xticklabels({32','28','24','20','16'})
ax = gca; 
ax.XAxis.FontSize=9;
% xticks([0 1/32   1/16])
% xticklabels({'0','1/32','1/16'})
%set(gca,'TickLength',[0.03, 0.01])

xlabel('$1/N_s$','Interpreter','latex','FontSize',12);
ylabel('$E_n-E_0$ (meV)','interpreter','latex','fontsize',12)
%title("Spectral gap",'interpreter','latex','fontsize',12)

text(-0.015,2.2*4,'(b)','interpreter','latex','fontsize',14)

title("$\nu_{\rm{total}}=3+1/2$, QHF",'interpreter','latex','fontsize',13,'Color','r');

text(0.015,2.5,"$N_g=17\,15 \,13\,\,\,\,11 \,\,\,\,\,\,\,\,\,\,9$",'interpreter','latex','fontsize',12,'color','r')

text(0.013,-0.9,"$N_s=$",'interpreter','latex','fontsize',12,'color','k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0,pos_y0+0.5,length_x,length_y];
subplot('Position',pos);

Ees=[];
Nset=[];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24E_HF_2valley')
load("data/Es_24E_16_waa0.6.mat");
Ns=24;
%3,3,3,4
plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(2)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(3)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(4)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(4)-E_all(1)];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27_HF_2valley')
load("data/Es_27_18_waa0.6.mat");

Ns=27;
plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(2)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(3)-E_all(1),'ro','MarkerSize',7);hold on;
plot(1/Ns,E_all(4)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(4)-E_all(1)];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley')
load("data/Es_36_24_waa0.6.mat");

Ns=36;
f1=plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
f1=plot(1/Ns,E_all(2)-E_all(1),'ro','MarkerSize',7);hold on;
f1=plot(1/Ns,E_all(3)-E_all(1),'ro','MarkerSize',7);hold on;
f2=plot(1/Ns,E_all(4)-E_all(1),'bo','MarkerSize',7);hold on;

Ees=[Ees,E_all(4)-E_all(1)];
Nset=[Nset,Ns];


range=1:3;
xs=1./Nset(range);
ys=Ees(range);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'b--');hold on;


    % h=legend([f1;f2],'Ground state','Excitation');
    % set(h,'Interpreter','latex','Location','northwest')
    % set(h,'fontsize',13)
    % set(h,'units','normalized');
    % %set(h,'NumColumns',2);
    % h.ItemTokenSize = [15,10];
    % %h.NumColumns=2;

axis([0,1/20,0,3]);
set(gca,'fontsize',12)
xticks([1/36  1/27  1/24 ])
xticklabels({'36','27','24'})
ax = gca; 
ax.XAxis.FontSize=9;
set(gca,'TickLength',[0.03, 0.01])

xlabel('$1/N_s$','Interpreter','latex','FontSize',12);
ylabel('$E_n-E_0$ (meV)','interpreter','latex','fontsize',12)
title("$\nu_{\rm{total}}=3+2/3$, QHC",'interpreter','latex','fontsize',13,'Color','r');


text(1/45,0.5,'$N_g=3$',  'Interpreter', 'latex', 'fontsize', 12,'color','r');
text(-0.01,2.2*1.5,'(a)','interpreter','latex','fontsize',14)

text(0.011,-0.33,"$N_s=$",'interpreter','latex','fontsize',12,'color','k')

%text(0.023,-0.75,"$1/N_s=$",'interpreter','latex','fontsize',12,'color','k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+2/3,pos_y0,length_x,length_y*0.95];
subplot('Position',pos);
restoredefaultpath;

Sqmaxs=[];
Nset=[];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley\pair_correl')
load("data/Sq_24_8_waa0.6.mat");
Ns=24;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27_HF_2valley\pair_correl')
load("data/Sq_27_9_waa0.6.mat");
Ns=27;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N30A_HF_2valley\pair_correl')
load("data/Sq_30_10_waa0.6");
Ns=30;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley\paircorrel')
load("data/Sq_36_12_waa0.6");
Ns=36;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];


range=1:4;
xs=1./Nset(range);
ys=Sqmaxs(range);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'r--');hold on;

axis([0,1/20,0,3e-3]);
set(gca,'fontsize',12)
% xticks([0 1/36 1/30 1/27  1/24 ])
% xticklabels({'0','1/36','1/30', '1/27','1/24'})
xticks([1/36 1/30  1/27 1/24 ])
xticklabels({36','30','27','24'})
ax = gca; 
ax.XAxis.FontSize=9;
% xticks([0 1/36  1/24 ])
% xticklabels({'0','1/36','1/24'})
set(gca,'TickLength',[0.03, 0.01])

xlabel('$1/N_s$','Interpreter','latex','FontSize',12);
ylabel('$\max\{S(\mathbf{q})\}/N_s$','interpreter','latex','fontsize',12)
%title("Spectra",'interpreter','latex','fontsize',12)

text(-0.01,(3e-3)*1.1,'(f)','interpreter','latex','fontsize',14)

% text(0.024,0.4,'$N_g=\,\,\,\,3 \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,3 \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,3$','interpreter','latex','fontsize',12,'color','r')
% text(0.0215,0.3,'$3$ ground states','interpreter','latex','fontsize',12,'color','r')
text(0.02,0.5,'$N_g=3$','interpreter','latex','fontsize',12,'color','r')

text(0.021,-0.23,"$N_s=$",'interpreter','latex','fontsize',12,'color','k')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+1/3,pos_y0,length_x,length_y*0.95];
subplot('Position',pos);
restoredefaultpath;

Sqmaxs=[];
Nset=[];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N16_HF_2valley\pair_correl')
load("data/Sq_16_8_waa0.6.mat");
Ns=16;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N20B_HF_2valley\paircorrel')
load("data/Sq_20_10_waa0.6.mat");
Ns=20;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley\pair_correl')
load("data/Sq_24_12_waa0.6");
Ns=24;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\paircorrel')
load("data/Sq_28_14_waa0.6");
Ns=28;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N32_HF_2valley\pair_correl')
load("data/Sq_32_16_waa0.6.mat");
Ns=32;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];


range=2:5;
xs=1./Nset(range);
ys=Sqmaxs(range);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'r--');hold on;

axis([0,1/15,0,3e-3]);
set(gca,'fontsize',12)
% xticks([0 1/36 1/30 1/27  1/24 ])
% xticklabels({'0','1/36','1/30', '1/27','1/24'})
xticks([1/32 1/28 1/24  1/20 1/16 ])
xticklabels({32','28', '24','20','16'})
ax = gca; 
ax.XAxis.FontSize=9;
% xticks([0 1/36  1/24 ])
% xticklabels({'0','1/36','1/24'})
set(gca,'TickLength',[0.03, 0.01])

xlabel('$1/N_s$','Interpreter','latex','FontSize',12);
ylabel('$\max\{S(\mathbf{q})\}/N_s$','interpreter','latex','fontsize',12)
%title("Spectra",'interpreter','latex','fontsize',12)

text(-0.015,(3e-3)*1.1,'(e)','interpreter','latex','fontsize',14)

% text(0.024,0.4,'$N_g=\,\,\,\,3 \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,3 \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,3$','interpreter','latex','fontsize',12,'color','r')
% text(0.0215,0.3,'$3$ ground states','interpreter','latex','fontsize',12,'color','r')
text(0.02,0.5,'$N_g=3$','interpreter','latex','fontsize',12,'color','r')

text(0.021,-0.23,"$N_s=$",'interpreter','latex','fontsize',12,'color','k')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,pos_y0,length_x,length_y*0.95];
subplot('Position',pos);
restoredefaultpath;


Sqmaxs=[];
Nset=[];


%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24E_HF_2valley\pair_correl')
load("data/Sq_24_16_waa0.6");
Ns=24;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27_HF_2valley\pair_correl')
load("data/Sq_27_18_waa0.6");
Ns=27;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley\paircorrel')
load("data/Sq_36_24_waa0.6.mat");
Ns=36;
Sq_values=sort(Sq_values);
plot(1/Ns,Sq_values(end)/Ns,'ro','MarkerSize',7);hold on;

Sqmaxs=[Sqmaxs,Sq_values(end)/Ns];
Nset=[Nset,Ns];


range=1:3;
xs=1./Nset(range);
ys=Sqmaxs(range);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'r--');hold on;

axis([0,1/20,0,3e-3]);
set(gca,'fontsize',12)
% xticks([0 1/36 1/30 1/27  1/24 ])
% xticklabels({'0','1/36','1/30', '1/27','1/24'})
xticks([1/36 1/27  1/24 ])
xticklabels({'36','27','24'})
ax = gca; 
ax.XAxis.FontSize=9;
% xticks([0 1/36  1/24 ])
% xticklabels({'0','1/36','1/24'})
set(gca,'TickLength',[0.03, 0.01])

xlabel('$1/N_s$','Interpreter','latex','FontSize',12);
ylabel('$\max\{S(\mathbf{q})\}/N_s$','interpreter','latex','fontsize',12)
%title("Spectra",'interpreter','latex','fontsize',12)

text(-0.01,(3e-3)*1.1,'(d)','interpreter','latex','fontsize',14)

% text(0.024,0.4,'$N_g=\,\,\,\,3 \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,3 \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,3$','interpreter','latex','fontsize',12,'color','r')
% text(0.0215,0.3,'$3$ ground states','interpreter','latex','fontsize',12,'color','r')
text(0.02,0.5,'$N_g=3$','interpreter','latex','fontsize',12,'color','r')

text(0.021,-0.23,"$N_s=$",'interpreter','latex','fontsize',12,'color','k')


%%%%%%%%%%%%%%%%
set(gcf,'Position',[100 100 850 380])




exportgraphics(gcf,'scaling_all.eps','ContentType','vector')



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