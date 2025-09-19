clc;clear;close all;
restoredefaultpath
pos_x0=0.07;
pos_y0=0.115;
length_x=0.24;
length_y=0.18;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0,2/3+pos_y0,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27');
Kg_set=[6+1,18+1,21+1];
Ng_set=[1,1,1];

waa_set=0:0.05:0.8;

E_global=zeros(10,length(waa_set));

Eg_Ka=[];
Eg_Kb=[];
Eg_Kc=[];
Es_outside=[];

epsilon=4;
for cc=1:length(waa_set)
    filenm="e1theta1.44_epsilon"+num2str(epsilon)+"_waa"+num2str(waa_set(cc))+"_Np18_flx_0_0.mat";
    load("Ns27/"+filenm);
    Eg=min(sortcell(En_set));
    E_all=sortcell(En_set);
    E_global(:,cc)=E_all(1:size(E_global,1));


    Eg_Ka=[Eg_Ka,En_set{Kg_set(1)}(1:8)];
    Eg_Kb=[Eg_Kb,En_set{Kg_set(2)}(1:8)];
    Eg_Kc=[Eg_Kc,En_set{Kg_set(3)}(1:8)];

    E_outside=[];
    for dd=1:length(En_set)
        if length(find(Kg_set==dd))==0
            
            E_outside=[E_outside;En_set{dd}(1:8)];
        end
    end
    E_outside=sort(E_outside);
    Es_outside=[Es_outside,E_outside];

end


Eg=[Eg_Ka(1:1,:);Eg_Kb(1:1,:);Eg_Kc(1:1,:);];
Eg=min(Eg,[],1);

Eexcitation=min([Eg_Ka(2,:);Eg_Kb(2,:);Eg_Kc(2,:);Es_outside(1,:)],[],1);


patch_x1=0;
patch_x2=0.6;
patch_y1=-2;
patch_y2=10;
fn=patch([patch_x1 patch_x2,patch_x2 patch_x1], [patch_y1,patch_y1,patch_y2,patch_y2], [0.8 0.8 0.8],'LineStyle','none');hold on;


f1=plot(waa_set,Eg_Ka(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kb(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kc(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f9=plot(waa_set,Eexcitation-Eg,'k--','LineWidth',1.5);hold on;

    h=legend([f1;f9],'Ground states','Excitations');
    set(h,'Interpreter','latex','Location','northeast')
    set(h,'fontsize',12)
    set(h,'units','normalized');
    %set(h,'NumColumns',2);
    h.ItemTokenSize = [15,10];
    
set(gca,'fontsize',12);
ylabel('E (meV)','Interpreter','latex','FontSize',15);
xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);

title(['(a) $\nu=2/3$'],'interpreter','latex','fontsize',15)
axis([0,0.8,-0,3])
xticks([0 0.2,0.4,0.6,0.8])
text(0.2,1,'QHC',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0,1/3+0.02+pos_y0,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27_HF_2valley');
Kg_set=[6+1,18+1,21+1];
Ng_set=[1,1,1];

waa_set=0:0.05:0.8;

E_global=zeros(10,length(waa_set));

Eg_Ka=[];
Eg_Kb=[];
Eg_Kc=[];
Es_outside=[];

epsilon=4;
for cc=1:length(waa_set)
    filenm="e1theta1.44_epsilon"+num2str(epsilon)+"_waa"+num2str(waa_set(cc))+"_Np18_flx_0_0.mat";
    load("Ns27_HF/"+filenm);
    Eg=min(sortcell(En_set));
    E_all=sortcell(En_set);
    E_global(:,cc)=E_all(1:size(E_global,1));


    Eg_Ka=[Eg_Ka,En_set{Kg_set(1)}(1:8)];
    Eg_Kb=[Eg_Kb,En_set{Kg_set(2)}(1:8)];
    Eg_Kc=[Eg_Kc,En_set{Kg_set(3)}(1:8)];

    E_outside=[];
    for dd=1:length(En_set)
        if length(find(Kg_set==dd))==0
            
            E_outside=[E_outside;En_set{dd}(1:8)];
        end
    end
    E_outside=sort(E_outside);
    Es_outside=[Es_outside,E_outside];

end


Eg=[Eg_Ka(1:1,:);Eg_Kb(1:1,:);Eg_Kc(1:1,:);];
Eg=min(Eg,[],1);

Eexcitation=min([Eg_Ka(2,:);Eg_Kb(2,:);Eg_Kc(2,:);Es_outside(1,:)],[],1);

patch_x1=0;
patch_x2=0.68;
patch_y1=-2;
patch_y2=10;
fn=patch([patch_x1 patch_x2,patch_x2 patch_x1], [patch_y1,patch_y1,patch_y2,patch_y2], [0.8 0.8 0.8],'LineStyle','none');hold on;



f1=plot(waa_set,Eg_Ka(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kb(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kc(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f9=plot(waa_set,Eexcitation-Eg,'k--','LineWidth',1.5);hold on;


    
set(gca,'fontsize',12);
ylabel('E (meV)','Interpreter','latex','FontSize',15);
xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);

title(['(d) $\nu_{\rm{total}}=3+2/3$'],'interpreter','latex','fontsize',15)
axis([0,0.8,-0,3])
xticks([0 0.2,0.4,0.6,0.8])
text(0.2,0.7,'QHC',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+1/3,2/3+pos_y0,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28\spectrum_scan_waa')
Kg_set=[0+1,1+1,11+1,15+1];
Ng_set=[3,4,4,4];

waa_set=0:0.05:0.95;

E_global=zeros(10,length(waa_set));

Eg_Ka=[];
Eg_Kb=[];
Eg_Kc=[];
Eg_Kd=[];
Es_outside=[];

epsilon=4;
for cc=1:length(waa_set)
    if waa_set(cc)==0
        waa_str="0.0";
    else
        waa_str=num2str(waa_set(cc));
    end
    N=28;
    En_set=cell(N);
    filenm0="_e1theta1.44_epsilon"+num2str(epsilon)+"_waa"+waa_str+"_Np14_flx_0_0.mat";
    for dd=1:N
    filenm="Ns28/"+"Kind"+num2str(dd)+filenm0;
    load(filenm);
    En_set{dd}=eu;
    end


    Eg=min(sortcell(En_set));
    E_all=sortcell(En_set);
    E_global(:,cc)=E_all(1:size(E_global,1));


    Eg_Ka=[Eg_Ka,En_set{Kg_set(1)}(1:8)];
    Eg_Kb=[Eg_Kb,En_set{Kg_set(2)}(1:8)];
    Eg_Kc=[Eg_Kc,En_set{Kg_set(3)}(1:8)];
    Eg_Kd=[Eg_Kd,En_set{Kg_set(4)}(1:8)];

    E_outside=[];
    for dd=1:length(En_set)
        if length(find(Kg_set==dd))==0
            
            E_outside=[E_outside;En_set{dd}(1:8)];
        end
    end
    E_outside=sort(E_outside);
    Es_outside=[Es_outside,E_outside];

end


Eg=[Eg_Ka(1:3,:);Eg_Kb(1:4,:);Eg_Kc(1:4,:);Eg_Kd(1:4,:)];
Eg=min(Eg,[],1);

Eexcitation=min([Eg_Ka(4,:);Eg_Kb(5,:);Eg_Kc(5,:);Eg_Kd(5,:);Es_outside(1,:)],[],1);

patch_x1=0;
patch_x2=0.85;
patch_y1=-2;
patch_y2=10;
fn=patch([patch_x1 patch_x2,patch_x2 patch_x1], [patch_y1,patch_y1,patch_y2,patch_y2], [0.8 0.8 0.8],'LineStyle','none');hold on;



f1=plot(waa_set,Eg_Ka(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f1=plot(waa_set,Eg_Ka(2,:)-Eg,'ro-','MarkerSize',6);hold on;
f1=plot(waa_set,Eg_Ka(3,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kb(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kb(2,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kb(3,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kb(4,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kc(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kc(2,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kc(3,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kc(4,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kd(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kd(2,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kd(3,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kd(4,:)-Eg,'ro-','MarkerSize',6);hold on;

f9=plot(waa_set,Eexcitation-Eg,'k--','LineWidth',1.5);hold on;


    
set(gca,'fontsize',12);
axis([0,1,-0,8])
xticks([0 0.2,0.4,0.6,0.8])
yticks([0,4,8])
xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);

title(['(b) $\nu=1/2$'],'interpreter','latex','fontsize',15)
text(0.2,2,'QHF',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+1/3,1/3+0.02+pos_y0,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\spectrum_scan_waa')
Kg_set=[0+1,1+1,11+1,15+1];
Ng_set=[3,4,4,4];

waa_set=0:0.05:0.95;

E_global=zeros(10,length(waa_set));

Eg_Ka=[];
Eg_Kb=[];
Eg_Kc=[];
Eg_Kd=[];
Es_outside=[];


for cc=1:length(waa_set)
    if waa_set(cc)==0
        waa_str="0.0";
    else
        waa_str=num2str(waa_set(cc));
    end
    N=28;
    En_set=cell(N);
    filenm0="_e1theta1.44_epsilon"+num2str(epsilon)+"_waa"+waa_str+"_Np14_flx_0_0.mat";
    for dd=1:N
    filenm="Ns28_HF/"+"Kind"+num2str(dd)+filenm0;
    load(filenm);
    En_set{dd}=eu;
    end


    Eg=min(sortcell(En_set));
    E_all=sortcell(En_set);
    E_global(:,cc)=E_all(1:size(E_global,1));


    Eg_Ka=[Eg_Ka,En_set{Kg_set(1)}(1:8)];
    Eg_Kb=[Eg_Kb,En_set{Kg_set(2)}(1:8)];
    Eg_Kc=[Eg_Kc,En_set{Kg_set(3)}(1:8)];
    Eg_Kd=[Eg_Kd,En_set{Kg_set(4)}(1:8)];

    E_outside=[];
    for dd=1:length(En_set)
        if length(find(Kg_set==dd))==0
            
            E_outside=[E_outside;En_set{dd}(1:8)];
        end
    end
    E_outside=sort(E_outside);
    Es_outside=[Es_outside,E_outside];

end


Eg=[Eg_Ka(1:3,:);Eg_Kb(1:4,:);Eg_Kc(1:4,:);Eg_Kd(1:4,:)];
Eg=min(Eg,[],1);

Eexcitation=min([Eg_Ka(4,:);Eg_Kb(5,:);Eg_Kc(5,:);Eg_Kd(5,:);Es_outside(1,:)],[],1);

patch_x1=0;
patch_x2=0.89;
patch_y1=-2;
patch_y2=10;
fn=patch([patch_x1 patch_x2,patch_x2 patch_x1], [patch_y1,patch_y1,patch_y2,patch_y2], [0.8 0.8 0.8],'LineStyle','none');hold on;



f1=plot(waa_set,Eg_Ka(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f1=plot(waa_set,Eg_Ka(2,:)-Eg,'ro-','MarkerSize',6);hold on;
f1=plot(waa_set,Eg_Ka(3,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kb(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kb(2,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kb(3,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kb(4,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kc(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kc(2,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kc(3,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kc(4,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kd(1,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kd(2,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kd(3,:)-Eg,'ro-','MarkerSize',6);hold on;
f3=plot(waa_set,Eg_Kd(4,:)-Eg,'ro-','MarkerSize',6);hold on;

f9=plot(waa_set,Eexcitation-Eg,'k--','LineWidth',1.5);hold on;


    
set(gca,'fontsize',12);
axis([0,1,-0,8])
xticks([0 0.2,0.4,0.6,0.8])
yticks([0,4,8])
xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);

title(['(e) $\nu_{\rm{total}}=3+1/2$'],'interpreter','latex','fontsize',15)
text(0.2,2,'QHF',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+2/3,2/3+pos_y0,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A')
Kg_set=[6+1,13+1,20+1];
Ng_set=[1,1,1];

waa_set=0:0.05:0.8;

E_global=zeros(10,length(waa_set));

Eg_Ka=[];
Eg_Kb=[];
Eg_Kc=[];
Es_outside=[];

epsilon=4;

for cc=1:length(waa_set)
    filenm="e1theta1.44_epsilon"+num2str(epsilon)+"_waa"+num2str(waa_set(cc))+"_Np8_flx_0_0.mat";
    load("Ns24/"+filenm);
    Eg=min(sortcell(En_set));
    E_all=sortcell(En_set);
    E_global(:,cc)=E_all(1:size(E_global,1));


    Eg_Ka=[Eg_Ka,En_set{Kg_set(1)}(1:8)];
    Eg_Kb=[Eg_Kb,En_set{Kg_set(2)}(1:8)];
    Eg_Kc=[Eg_Kc,En_set{Kg_set(3)}(1:8)];

    E_outside=[];
    for dd=1:length(En_set)
        if length(find(Kg_set==dd))==0
            
            E_outside=[E_outside;En_set{dd}(1:8)];
        end
    end
    E_outside=sort(E_outside);
    Es_outside=[Es_outside,E_outside];

end


Eg=[Eg_Ka(1:1,:);Eg_Kb(1:1,:);Eg_Kc(1:1,:)];
Eg=min(Eg,[],1);

Eexcitation=min([Eg_Ka(2,:);Eg_Kb(2,:);Eg_Kc(2,:);Es_outside(1,:)],[],1);

patch_x1=0;
patch_x2=0.35;
patch_y1=-2;
patch_y2=10;
fn=patch([patch_x1 patch_x2,patch_x2 patch_x1], [patch_y1,patch_y1,patch_y2,patch_y2], [0.8 0.8 0.8],'LineStyle','none');hold on;



f1=plot(waa_set,Eg_Ka(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kb(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kc(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f9=plot(waa_set,Eexcitation-Eg,'k--','LineWidth',1.5);hold on;


    

set(gca,'fontsize',12);
xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);

title(['(c) $\nu=1/3$'],'interpreter','latex','fontsize',15)
axis([0,0.8,-0,3])
xticks([0 0.2,0.4,0.6,0.8])
text(0.05,0.35,'FCI',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+2/3,1/3+0.02+pos_y0,length_x,length_y];
subplot('Position',pos);


%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley')
Kg_set=[6+1,13+1,20+1];
Ng_set=[1,1,1];

waa_set=0:0.05:0.8;

E_global=zeros(10,length(waa_set));

Eg_Ka=[];
Eg_Kb=[];
Eg_Kc=[];
Es_outside=[];



for cc=1:length(waa_set)
    filenm="e1theta1.44_epsilon"+num2str(epsilon)+"_waa"+num2str(waa_set(cc))+"_Np8_flx_0_0.mat";
    load("Ns24_HF/"+filenm);
    Eg=min(sortcell(En_set));
    E_all=sortcell(En_set);
    E_global(:,cc)=E_all(1:size(E_global,1));


    Eg_Ka=[Eg_Ka,En_set{Kg_set(1)}(1:8)];
    Eg_Kb=[Eg_Kb,En_set{Kg_set(2)}(1:8)];
    Eg_Kc=[Eg_Kc,En_set{Kg_set(3)}(1:8)];

    E_outside=[];
    for dd=1:length(En_set)
        if length(find(Kg_set==dd))==0
            
            E_outside=[E_outside;En_set{dd}(1:8)];
        end
    end
    E_outside=sort(E_outside);
    Es_outside=[Es_outside,E_outside];

end


Eg=[Eg_Ka(1:1,:);Eg_Kb(1:1,:);Eg_Kc(1:1,:)];
Eg=min(Eg,[],1);

Eexcitation=min([Eg_Ka(2,:);Eg_Kb(2,:);Eg_Kc(2,:);Es_outside(1,:)],[],1);


patch_x1=0;
patch_x2=0.74;
patch_y1=-2;
patch_y2=10;
fn=patch([patch_x1 patch_x2,patch_x2 patch_x1], [patch_y1,patch_y1,patch_y2,patch_y2], [0.8 0.8 0.8],'LineStyle','none');hold on;


f1=plot(waa_set,Eg_Ka(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kb(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f3=plot(waa_set,Eg_Kc(1,:)-Eg,'ro-','MarkerSize',6);hold on;

f9=plot(waa_set,Eexcitation-Eg,'k--','LineWidth',1.5);hold on;


    

set(gca,'fontsize',12);
xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);

title(['(f) $\nu_{\rm{total}}=3+1/3$'],'interpreter','latex','fontsize',15)
axis([0,0.8,-0,3])
xticks([0 0.2,0.4,0.6,0.8])
text(0.2,0.35,'FCI',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0,pos_y0+0.04,length_x,length_y];
subplot('Position',pos);
fig=gcf;

pathnm="D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27_HF_2valley\scan_waa_epsilon\";
load(pathnm+"Egap_waa_epsilon_Np18.mat");

Egap=Eexcitation_set-Eg_max_set;
%Egap=Eexcitation_set-Eg_set;

value_range=[0,5];

my_colormap1=hot(150);
my_colormap1=my_colormap1(end:-1:1,:);
my_colormap2=abyss(200);
my_colormap2=my_colormap2(1:end,:);
my_colormap=[my_colormap1;my_colormap2];

% pos=[pos_x0,pos_y0,length_x,length_y];
% subplot('Position',pos);

Egap=treat_negative(Egap);
imagesc('XData',waa_set,'YData',epsilon_set,'CData',permute(Egap,[2,1]));
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,2,8]);
% set(gca,'xticklabel',{[]})
% set(gca,'yticklabel',{[]})
%text(-0.11,1.5,'$\nu=$',  'Interpreter', 'latex', 'fontsize', 17);
%text(-0.4,0.7,'QHC, $\nu=2/3$',  'Interpreter', 'latex', 'fontsize', 12);
set(gca,'fontsize',12);
xticks([0 0.2,0.4,0.6,0.8])

title(['(g) $\nu_{\rm{total}}=3+2/3$'],'interpreter','latex','fontsize',15)

set(gca,'TickDir','out');
set(gca,'TickLength',[0.03, 0.01])

xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);
ylabel('$\epsilon$','Interpreter','latex','FontSize',18);


h = axes(fig,'visible','off'); 
% h.Title.Visible = 'on';
% h.XLabel.Visible = 'on';
% h.YLabel.Visible = 'on';
% ylabel(h,'yaxis','FontWeight','bold');
% xlabel(h,'xaxis','FontWeight','bold');
%title(h,'title');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[1/3+pos_x0,pos_y0+0.04,length_x,length_y];
subplot('Position',pos);
fig=gcf;

pathnm="D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\scan_waa_epsilon\";
load(pathnm+"Egap_waa_epsilon_Np14.mat");


Egap=Eexcitation_set-Eg_max_set;
%Egap=Eexcitation_set-Eg_set;

value_range=[0,5];

my_colormap1=hot(150);
my_colormap1=my_colormap1(end:-1:1,:);
my_colormap2=abyss(200);
my_colormap2=my_colormap2(1:end,:);
my_colormap=[my_colormap1;my_colormap2];

% pos=[pos_x0,pos_y0,length_x,length_y];
% subplot('Position',pos);

Egap=treat_negative(Egap);
imagesc('XData',waa_set,'YData',epsilon_set,'CData',permute(Egap,[2,1]));
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,2,10]);
% set(gca,'xticklabel',{[]})
% set(gca,'yticklabel',{[]})
%text(-0.11,1.5,'$\nu=$',  'Interpreter', 'latex', 'fontsize', 17);
%text(-0.4,0.7,'QHC, $\nu=2/3$',  'Interpreter', 'latex', 'fontsize', 12);

set(gca,'fontsize',12);
xticks([0 0.2,0.4,0.6,0.8])
yticks([2,6,10])

set(gca,'TickDir','out');
set(gca,'TickLength',[0.03, 0.01])

title(['(h) $\nu_{\rm{total}}=3+1/2$'],'interpreter','latex','fontsize',15)

xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);
%ylabel('$\epsilon$','Interpreter','latex','FontSize',18);


h = axes(fig,'visible','off'); 
% h.Title.Visible = 'on';
% h.XLabel.Visible = 'on';
% h.YLabel.Visible = 'on';
% ylabel(h,'yaxis','FontWeight','bold');
% xlabel(h,'xaxis','FontWeight','bold');
%title(h,'title');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[2/3+pos_x0,pos_y0+0.04,length_x,length_y];
subplot('Position',pos);
fig=gcf;

pathnm="D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley\scan_waa_epsilon\";
load(pathnm+"Egap_waa_epsilon_Np8.mat");


Egap=Eexcitation_set-Eg_max_set;
%Egap=Eexcitation_set-Eg_set;

value_range=[0,5];

my_colormap1=hot(150);
my_colormap1=my_colormap1(end:-1:1,:);
my_colormap2=abyss(200);
my_colormap2=my_colormap2(1:end,:);
my_colormap=[my_colormap1;my_colormap2];

% pos=[pos_x0,pos_y0,length_x,length_y];
% subplot('Position',pos);

Egap=treat_negative(Egap);
imagesc('XData',waa_set,'YData',epsilon_set,'CData',permute(Egap,[2,1]));
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,2,8]);
% set(gca,'xticklabel',{[]})
% set(gca,'yticklabel',{[]})
%text(-0.11,1.5,'$\nu=$',  'Interpreter', 'latex', 'fontsize', 17);
%text(-0.4,0.7,'QHC, $\nu=2/3$',  'Interpreter', 'latex', 'fontsize', 12);

set(gca,'fontsize',12);
xticks([0 0.2,0.4,0.6,0.8])

title(['(i) $\nu_{\rm{total}}=3+1/3$'],'interpreter','latex','fontsize',15)
text(0.3,-1.8,'Gap (meV)',  'Interpreter', 'latex', 'fontsize', 15);

set(gca,'TickDir','out');
set(gca,'TickLength',[0.03, 0.01])

xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);
%ylabel('$\epsilon$','Interpreter','latex','FontSize',18);


h = axes(fig,'visible','off'); 

cb=colorbar(h,'Position',[0.87 0.5 0.03 0.4]);
colormap(cb,my_colormap)
caxis(h,value_range); 
%ylabel(cb,'Gap (meV)','Rotation',270,'interpreter','latex','fontsize',13)
%xlabel(cb,'Gap (meV)','Rotation',0,'interpreter','latex','fontsize',13)

cb.Location= 'southoutside';
cb.Position=[0.1 0.04 0.7 0.02];
%cb.Label.HorizontalAlignment = 'right';
cb.FontSize=12;




set(gcf,'Position',[100 100 670 600])



exportgraphics(gcf,'gap_data.eps','ContentType','vector')







function ys=treat_negative(ys)
for cc=1:length(ys)
    if ys(cc)<0
        ys(cc)=0;
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