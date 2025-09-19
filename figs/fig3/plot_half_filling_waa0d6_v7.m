clc;clear;close all;
restoredefaultpath;
pos_x0=0.1;
pos_y0=0.12;
length_x=0.36;
length_y=0.32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,0.5+pos_y0,length_x,length_y];
subplot('Position',pos);

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley')

N=28;
En_set=cell(N);
filenm0="_e1theta1.44_epsilon4_waa0.6_Np14_flx_0_0.mat";
for cc=1:N
filenm="Kind"+num2str(cc)+filenm0;
load(filenm);
En_set{cc}=eu;
end



V1=[-2,6];
V2=[4,2];
[Pset,N,g1,g3,kk1,kk2,V1_,V2_]=sort_K(V1,V2,k_total_set);
assert(N==28); 
N1=2;
N2=14;
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
axis([-2,N,-1,14]);

xlabel('$k_1+N_1 k_2$','interpreter','latex','fontsize',13)
ylabel('$E_n-E_0$ (meV)','interpreter','latex','fontsize',13)

% text(5,0.3,'$\times$ $2$ each ($6$ total) ',  'Interpreter', 'latex', 'fontsize', 12);
text(-1.5,2,['$N_{g}=3+4+4+4=N_p+1$'],  'Interpreter', 'latex', 'fontsize', 11.5,'color','r');
% text(1,13,'(a) $N_s=32$','interpreter','latex','fontsize',13)
text(-7,15,'(a)','interpreter','latex','fontsize',14)

rectangle('position',[0-0.4,-0.3,20,0.6],'edgecolor','r')

set(gca,'fontsize',12)
xticks([0 7 14 21 28 ])
xticklabels({'0','7', '14','21', '28'})
% title("Hall crystal at $\nu=1/2$",'interpreter','latex')
%title("Spectrum",'interpreter','latex')
title("$\nu_{\rm{total}}=3+\frac{1}{2}$, QHF",'interpreter','latex','fontsize',12,'Color','r')
set(gca,'TickLength',[0.03, 0.01])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,pos_y0,length_x*0.9,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\order')
load order_theta1.44_waa0.6_epsilon4_Np14;



sec=1;
N=size(k_set,1);
pos=sum(n_keep_all(1:sec-1))+1:sum(n_keep_all(1:sec));
nks=zeros(N,n_keep_all(sec));
for cc=1:length(pos)
    nks(:,cc)=diag(ob_all(:,:,pos(cc),pos(cc)))';
end


plot(sort(nks(:,1)),'b.','MarkerSize',10,'LineWidth',0.5);hold on;





load 3_4.mat;
N=size(k_set,1);
eus=eig(rho);
plot(sort(eus),'r.','MarkerSize',10);hold on;
axis([1,N,0,1]);
title('$\mathbf{Q}=\mathbf{M}_2$','interpreter','latex','FontSize',12)

text(10,0.7,'$\mathbf{\Phi}_{k}$','interpreter','latex','fontsize',14,'Color','b')

text(7,0.15,'$a \mathbf{\Phi}_{k} + b \mathbf{\Phi}_{k+\mathbf{Q}}$','interpreter','latex','fontsize',14,'color','r')








axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)

title('Orbital occupation','interpreter','latex','FontSize',13)
text(-6,1.2,"(c)",'interpreter','latex','FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+0.5,pos_y0+0.5,length_x,length_y*0.95];
subplot('Position',pos);



%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N16_HF_2valley')
load Es_16_8_waa0.6.mat;
Ns=16;
%3,2,2,2
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(9)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:9
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(10)-E_all(1),'bo','MarkerSize',7);hold on;


%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N20B_HF_2valley')
load Es_20_10_waa0.6.mat;
Ns=20;
%3,2,2,2
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(9)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:11
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(12)-E_all(1),'bo','MarkerSize',7);hold on;



%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley')
load Es_24A_12_waa0.6.mat;
Ns=24;
%3,3,3,4
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(13)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:13
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(14)-E_all(1),'bo','MarkerSize',7);hold on;




%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley')
load Es_28_14_waa0.6.mat;
Ns=28;
%3,4,4,4
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(15)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:15
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(16)-E_all(1),'bo','MarkerSize',7);hold on;



%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N32_HF_2valley')
load Es_32_16_waa0.6.mat;
Ns=32;
%5,4,4,4
% plot(1/Ns,E_all(1)-E_all(1),'ro','MarkerSize',7);hold on;
% plot(1/Ns,E_all(17)-E_all(1),'ro','MarkerSize',7);hold on;
for cc=1:17
    plot(1/Ns,E_all(cc)-E_all(1),'ro','MarkerSize',7);hold on;
end
plot(1/Ns,E_all(18)-E_all(1),'bo','MarkerSize',7);hold on;


%axis([0.0,1/16,0,8]);
axis([0.01,1/16,0,8.2]);
set(gca,'fontsize',12)
% xticks([0 1/36 1/30 1/27  1/24 ])
% xticklabels({'0','1/36','1/30', '1/27','1/24'})
yticks([0 4 8])
xticks([0 1/32  1/28  1/24 1/20 1/16])
xticklabels({'0','32','28','24','20','16'})
ax = gca; 
ax.XAxis.FontSize=9;
% xticks([0 1/32   1/16])
% xticklabels({'0','1/32','1/16'})
%set(gca,'TickLength',[0.03, 0.01])

xlabel('$1/N_s$','Interpreter','latex','FontSize',12);
ylabel('$E_n-E_0$ (meV)','interpreter','latex','fontsize',12)
%title("Spectral gap",'interpreter','latex','fontsize',12)

text(0.01,9.2,'(b)','interpreter','latex','fontsize',14)

title("Spectra",'interpreter','latex','FontSize',12)

text(0.017,2.5,"$N_g=17\,\,15 \,\,13\,\,\,\,\,\,\,11 \,\,\,\,\,\,\,\,\,\,\,9$",'interpreter','latex','fontsize',12,'color','r')

text(0.017,-0.9,"$N_s=$",'interpreter','latex','fontsize',12,'color','k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+0.5,pos_y0,length_x,length_y];
subplot('Position',pos);


Ns=[16,20,24,28,32];

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N16_HF_2valley\order')
load("Ns16\Sz_waa0.6_option2.mat");
Sz_16=real(eig(Op_z));

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N20B_HF_2valley\order')
load("Ns20\Sz_waa0.6_option2.mat");
Sz_20=real(eig(Op_z));

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley\order_Np12')
load("Ns24\Sz_waa0.6_option2.mat");
Sz_24=real(eig(Op_z));

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\order')
load("Ns28\Sz_waa0.6_option2.mat");
Sz_28=real(eig(Op_z));

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N32_HF_2valley\order')
load("Ns32\Sz_waa0.6_option2.mat");
Sz_32=real(eig(Op_z));



markerlength=30;
Markerwidth=0.6;
LW=2;

N=16;
Sz=Sz_16;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end

N=20;
Sz=Sz_20;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end

N=24;
Sz=Sz_24;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end

N=28;
Sz=Sz_28;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end

N=32;
Sz=Sz_32;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end

set(gca,'fontsize',12)
axis([15,35,-8.5,8.5])
%yticks([-8,-4,0,4,8])


yticks([-8,-6,-4,-2,0,2,4,6,8])
yticklabels({'-8','','-4','','0','','4','','8'})

xticks([16,20,24,28,32])
xlabel('$N_s$','interpreter','latex','FontSize',12)
%ylabel('$S^z$','interpreter','latex','FontSize',12)

text(10,8.8,'(d)','interpreter','latex','fontsize',14)

title("$\hat{S}^z$ eigenvalues",'interpreter','latex','FontSize',14)
set(gca,'TickLength',[0.03, 0.01])
%%%%%%%%%%%%%%%%%%%%







set(gcf,'Position',[100 100 550 380])


exportgraphics(gcf,'QHF_1_2.eps','ContentType','vector')



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