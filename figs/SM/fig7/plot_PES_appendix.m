clc;clear;close all;
restoredefaultpath;
%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N30A_HF_2valley\PES');
pos_x0=0.1;
pos_y0=0.12;
length_x=0.35;
length_y=0.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);
filenm='PES_theta1.44_waa0.3_epsilon4_Np10_NA_3_sec';
N=30;
PES_set=cell(N,1);
for cc=1:N
    fn=[filenm,num2str(cc),'.mat'];
    load(fn);
    PES_set{cc}=PES;
end

mu=8.8;
N_total=0;
for cc=1:length(PES_set)

es=PES_set{cc};
es=es(end:-1:1);
es=-log(es);
plot(cc-1,es,'k.');hold on;

N_total=N_total+length(find(mu>es));

end


%plot(0:length(PES_set),mu*ones(1,1+length(PES_set)),'r');hold on;
%title(['N_A=4,  counting=', num2str(N_total)])
title(['N_A=3'])


es_total=[];
for cc=1:length(PES_set)

es=PES_set{cc};
es=es(end:-1:1);
es=-log(es);
es_total=[es_total;es];

end

count=3250;
es_total=sort(es_total);
mu=(es_total(count)+es_total(count+1))/2;
plot(0:length(PES_set),mu*ones(1,1+length(PES_set)),'r','LineWidth',2);hold on;

set(gca,'fontsize',12)
xlabel('$k$','Interpreter','latex','FontSize',12);
ylabel('PES','interpreter','latex','fontsize',12)
title("$N_A=3$",'interpreter','latex','fontsize',12)

%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+0.5,pos_y0,length_x,length_y];
subplot('Position',pos);

filenm='PES_theta1.44_waa0.3_epsilon4_Np10_NA_4_sec';
N=30;
PES_set=cell(N,1);
for cc=1:N
    fn=[filenm,num2str(cc),'.mat'];
    load(fn);
    PES_set{cc}=PES;
end

mu=8.8;
N_total=0;
for cc=1:length(PES_set)

es=PES_set{cc};
es=es(end:-1:1);
es=-log(es);
plot(cc-1,es,'k.');hold on;

N_total=N_total+length(find(mu>es));

end


%plot(0:length(PES_set),mu*ones(1,1+length(PES_set)),'r');hold on;
%title(['N_A=4,  counting=', num2str(N_total)])
title(['N_A=4'])


es_total=[];
for cc=1:length(PES_set)

es=PES_set{cc};
es=es(end:-1:1);
es=-log(es);
es_total=[es_total;es];

end

es_total=sort(es_total);
mu=(es_total(17250)+es_total(17250+1))/2;
plot(0:length(PES_set),mu*ones(1,1+length(PES_set)),'r','LineWidth',2);hold on;

mu2=(es_total(24840)+es_total(24840+1))/2;
plot(0:length(PES_set),mu2*ones(1,1+length(PES_set)),'r','LineWidth',2);hold on;

set(gca,'fontsize',12)
xlabel('$k$','Interpreter','latex','FontSize',12);
ylabel('PES','interpreter','latex','fontsize',12)
title("$N_A=4$",'interpreter','latex','fontsize',12)



set(gcf,'Position',[100 100 450 550])


exportgraphics(gcf,'PES_appendix.eps','ContentType','vector')
