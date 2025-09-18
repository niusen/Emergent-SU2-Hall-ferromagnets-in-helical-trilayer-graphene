clc;clear;close all;
restoredefaultpath;
pos_x0=0.08;
pos_y0=0.13;
length_x=0.22;
length_y=0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,0.5+pos_y0,length_x,length_y];
subplot('Position',pos);

addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley')

N=36;
En_set=cell(N);
filenm0="_e1theta1.44_epsilon4_waa0.6_Np24_flx_0_0.mat";
for cc=1:N
filenm="Kind"+num2str(cc)+filenm0;
load(filenm);
En_set{cc}=eu;
end


V1=[6,0];
V2=[0,6];
[Pset,N,g1,g3,kk1,kk2,V1_,V2_]=sort_K(V1,V2,k_total_set);
assert(N==36); 
N1=6;
N2=6;
pp_set=zeros(N,2);
pp_ind=zeros(N,2);
for c1=0:N1-1
    for c2=0:N2-1
        pp_set(c1+N1*c2+1,:)=(c1*V2_+c2*V1_);
        pp_ind(c1+N1*c2+1,:)=[c1,c2];
    end
end

pp_mod=mod(pp_set,1);
for cc=1:size(Pset,1)
    if abs(pp_mod(cc,1)-1)<1e-10
        pp_mod(cc,1)=0;
    end
    if abs(pp_mod(cc,2)-1)<1e-10
        pp_mod(cc,2)=0;
    end
end



order=zeros(size(Pset,1),1);
for cc=1:size(Pset,1)
    kk=k_total_set(cc,:);
    cond1=abs(pp_mod(:,1)-kk(1))<1e-10;
    cond2=abs(pp_mod(:,2)-kk(2))<1e-10;
    pos=find(cond1.*cond2);
    
    order(cc)=pos;
end
order=order-1;

assert( length(unique(order))==N);



Eg=min(sortcell(En_set));
E_all=sortcell(En_set);
for c1=1:size(En_set,1)
    En=En_set{c1};
    Es=En-Eg;
    %Es=En;
    plot(order(c1),Es(1:6),'k.','markersize',12);hold on;
end
axis([-1,N,-1,5]);

xlabel('$k_1+N_1 k_2$','interpreter','latex','fontsize',13)
ylabel('$E_n-E_0$ (meV)','interpreter','latex','fontsize',13)

% text(5,0.3,'$\times$ $2$ each ($6$ total) ',  'Interpreter', 'latex', 'fontsize', 12);
text(4,0.75,['$N_{g}=1+1+1$'],  'Interpreter', 'latex', 'fontsize', 12,'color','r');
% text(1,13,'(a) $N_s=32$','interpreter','latex','fontsize',13)
text(-8,5.5,'(a)','interpreter','latex','fontsize',14)

rectangle('position',[0-0.4,-0.2,30,0.5],'edgecolor','r')

set(gca,'fontsize',12)
xticks([0 12  24  36 ])
xticklabels({'0','12','24','36'})
% title("Hall crystal at $\nu=2/3$",'interpreter','latex')
title("$\nu_{\rm{total}}=3+\frac{2}{3}$, QHC",'interpreter','latex','fontsize',12,'Color','r')
set(gca,'TickLength',[0.03, 0.01])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pos=[pos_x0+1/3,0.5+pos_y0,length_x,length_y];
subplot('Position',pos);


addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley\Chern')
load Chern_theta1.44_waa0.6_epsilon4_Np24_Lflx_6_K1.mat


flx_set;
F_sumx=-sum(imag(F_grid),1)/2/pi;


F_integral=zeros(1,length(flx_set));

for cc=2:length(flx_set)
    F_integral(cc)=sum(F_sumx(1:cc-1));
end


plot(flx_set,F_integral,'ko--');hold on;

% for cc=1:length(flx_set)
%     cc
%     plot(flx_set(cc),F_integral(cc),'ko');hold on;
% end


text(-0.3,0.05,'(b)','interpreter','latex','fontsize',14)

axis([0,1,-1,0]);



set(gca,'fontsize',12);
xlabel('Flux $\Theta_2/2\pi$','Interpreter','latex','FontSize',12);
% ylabel('$\int_{0}^{t}d\theta_2 \int_{0}^{2\pi}d\theta_1 F(\theta_1,\theta_2)$','Interpreter','latex','FontSize',18);
ylabel('$C(\Theta_2)$','Interpreter','latex','FontSize',12);
%text(0.3,-0.15,"$\mathbf{\Gamma}$ sector",'interpreter','latex','FontSize',12,'Rotation',-35)
title("Chern number",'interpreter','latex','FontSize',12)
text(0.2,-0.05,['1 state, $C_{\rm mean}=-1$'],  'Interpreter', 'latex', 'fontsize', 12,'color','r','Rotation',-35);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





pos=[pos_x0+2/3-0.04,pos_y0+0.5,length_x*1.2,length_y];
subplot('Position',pos);



load("D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\arXiv2405_08887\data\flatter\mymap.mat")

addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley\paircorrel')

load Sq_36_24_waa0.6.mat;





for cc=1:6
    p1=corners(cc,:);
    p2=corners(cc+1,:);
    plot([p1(1),p2(1)],[p1(2),p2(2)],'k-','LineWidth',6);hold on;
    
end



value_range=[min(real(Sq_values(:))),max(real(Sq_values(:)))];
value_range=[0,0.08];
for cq=1:length(Sq_values)
    
    cq
    Sq_values(cq)
    coord=Sq_points(cq,1)*g1+Sq_points(cq,2)*g3;
     
    Color=determine_color(value_range,real(Sq_values(cq)));
    plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',25);hold on;%colorbar;

end
text(-0.75,0.3,"$\mathbf{K}$",'interpreter','latex','FontSize',12)
text(-0.75,-0.3,"$\mathbf{K}'$",'interpreter','latex','FontSize',12)
text(-0.5,0.6,'$\mathbf{M}_1$','interpreter','latex','fontsize',12)
text(-0.83,0,'$\mathbf{M}_2$','interpreter','latex','fontsize',12)
text(-0.5,-0.6,'$\mathbf{M}_3$','interpreter','latex','fontsize',12)

text(0.04,-0.03,'$\mathbf{\Gamma}$','interpreter','latex','fontsize',12)
axis([-0.6,0.6,-0.6,0.6]);
title("$S(q)$",'interpreter','latex','FontSize',12)
%ylabel('Structure factor $S(\rm q)$','interpreter','latex','fontsize',13)
caxis(value_range);
ax=gca;
colormap(gca,mycmap);
cb=colorbar;
%ylabel(cb,"Structure factor $S(\rm{q})$",'Rotation',270,'interpreter','latex','fontsize',13)
% title("$S(q)$",'interpreter','latex')
text(-1,0.65,'(c)','interpreter','latex','fontsize',14)
%title("Structure factor $S(\rm{q})$",'interpreter','latex','FontSize',12)


ax.FontSize=12;
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);


addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley')

N=36;
En_set=cell(N);
filenm0="_e1theta1.44_epsilon4_waa0.6_Np12_flx_0_0.mat";
for cc=1:N
filenm="Kind"+num2str(cc)+filenm0;
load(filenm);
En_set{cc}=eu;
end


V1=[6,0];
V2=[0,6];
[Pset,N,g1,g3,kk1,kk2,V1_,V2_]=sort_K(V1,V2,k_total_set);
assert(N==36); 
N1=6;
N2=6;
pp_set=zeros(N,2);
pp_ind=zeros(N,2);
for c1=0:N1-1
    for c2=0:N2-1
        pp_set(c1+N1*c2+1,:)=(c1*V2_+c2*V1_);
        pp_ind(c1+N1*c2+1,:)=[c1,c2];
    end
end

pp_mod=mod(pp_set,1);
for cc=1:size(Pset,1)
    if abs(pp_mod(cc,1)-1)<1e-10
        pp_mod(cc,1)=0;
    end
    if abs(pp_mod(cc,2)-1)<1e-10
        pp_mod(cc,2)=0;
    end
end



order=zeros(size(Pset,1),1);
for cc=1:size(Pset,1)
    kk=k_total_set(cc,:);
    cond1=abs(pp_mod(:,1)-kk(1))<1e-10;
    cond2=abs(pp_mod(:,2)-kk(2))<1e-10;
    pos=find(cond1.*cond2);
    
    order(cc)=pos;
end
order=order-1;

assert( length(unique(order))==N);



Eg=min(sortcell(En_set));
E_all=sortcell(En_set);
for c1=1:size(En_set,1)
    En=En_set{c1};
    Es=En-Eg;
    %Es=En;
    plot(order(c1),Es(1:6),'k.','markersize',12);hold on;
end
axis([-1,N,-0.3,3]);

xlabel('$k_1+N_1 k_2$','interpreter','latex','fontsize',13)
ylabel('$E_n-E_0$ (meV)','interpreter','latex','fontsize',13)

% text(5,0.3,'$\times$ $2$ each ($6$ total) ',  'Interpreter', 'latex', 'fontsize', 12);
text(4,0.2,['$N_{g}=3$'],  'Interpreter', 'latex', 'fontsize', 12,'color','r');
% text(1,13,'(a) $N_s=32$','interpreter','latex','fontsize',13)
text(-8,3.5,'(d)','interpreter','latex','fontsize',14)

rectangle('position',[0-0.4,-0.1,2,0.5],'edgecolor','r')

set(gca,'fontsize',12)
xticks([0 12  24  36 ])
xticklabels({'0','12','24','36'})
% title("FQAH at $\nu=1/3$",'interpreter','latex')
title("$\nu_{\rm{total}}=3+\frac{1}{3}$, FCI",'interpreter','latex','fontsize',12,'Color','r')
set(gca,'TickLength',[0.03, 0.01])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+1/3,pos_y0,length_x,length_y];
subplot('Position',pos);

addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley\Chern')
load Chern_theta1.44_waa0.6_epsilon4_Np12_Lflx_6_K1.mat


flx_set;
F_sumx=sum(imag(-F_grid),1)/2/pi;


F_integral=zeros(1,length(flx_set));

for cc=2:length(flx_set)
    F_integral(cc)=sum(F_sumx(1:cc-1));
end


plot(flx_set,F_integral,'ko--');hold on;

% for cc=1:length(flx_set)
%     cc
%     plot(flx_set(cc),F_integral(cc),'ko');hold on;
% end

text(-0.3,0.068*2,'(e)','interpreter','latex','fontsize',14)

axis([0,1,-2,0]);



set(gca,'fontsize',12);
xlabel('Flux $\Theta_2/2\pi$','Interpreter','latex','FontSize',12);
% ylabel('$\int_{0}^{t}d\theta_2 \int_{0}^{2\pi}d\theta_1 F(\theta_1,\theta_2)$','Interpreter','latex','FontSize',18);
ylabel('$C(\Theta_2)$','Interpreter','latex','FontSize',12);
% title("$k_1+N_1k_2=0$",'interpreter','latex','FontSize',12)
title("Chern number",'interpreter','latex','FontSize',12)
%text(0.3,-0.15,"$\mathbf{\Gamma}$ sector",'interpreter','latex','FontSize',12,'Rotation',-35)
text(0.2,-0.05,['3 states, $C_{\rm mean}=-\frac{2}{3}$'],  'Interpreter', 'latex', 'fontsize', 12,'color','r','Rotation',-35);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+2/3-0.04,pos_y0,length_x*1.2,length_y];
subplot('Position',pos);

load("D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\arXiv2405_08887\data\flatter\mymap.mat")

addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N36_HF_2valley\paircorrel')

load Sq_36_12_waa0.6.mat;





for cc=1:6
    p1=corners(cc,:);
    p2=corners(cc+1,:);
    plot([p1(1),p2(1)],[p1(2),p2(2)],'k-','LineWidth',6);hold on;
    
end



value_range=[min(real(Sq_values(:))),max(real(Sq_values(:)))];
for cq=1:length(Sq_values)
    
    cq
    Sq_values(cq)
    coord=Sq_points(cq,1)*g1+Sq_points(cq,2)*g3;
     
    Color=determine_color(value_range,real(Sq_values(cq)));
    plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',25);hold on;%colorbar;

end
axis([-0.6,0.6,-0.6,0.6]);
title("$S(q)$",'interpreter','latex','FontSize',12)
caxis(value_range);
ax=gca;
colormap(gca,mycmap);
cb=colorbar;
%ylabel(cb,'Structure factor $S(\rm{q})$','Rotation',270,'interpreter','latex','fontsize',13)
% title("$S(q)$",'interpreter','latex')
text(-1,0.65,'(f)','interpreter','latex','fontsize',14)
%title("Structure factor $S(\rm{q})$",'interpreter','latex','FontSize',12)
%ylabel(cb,"Structure factor $S(\rm{q})$",'Rotation',270,'interpreter','latex','fontsize',13)

%text(0.05,-0.03,'$\mathbf{\Gamma}$','interpreter','latex','fontsize',14)

ax.FontSize=12;
axis off

%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,'Position',[100 100 650 380])


exportgraphics(gcf,'nu_2_3_and_1_3.eps','ContentType','vector')

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