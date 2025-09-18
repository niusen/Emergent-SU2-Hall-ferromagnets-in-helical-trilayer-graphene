clc;clear;close all;
restoredefaultpath;
pos_x0=0.12;
pos_y0=0.12;
length_x=0.35;
length_y=0.34;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,pos_y0+0.48,length_x,length_y];
subplot('Position',pos);
nK=2;
addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band')
stack="ABA";%"ABA","AAA"
a0 = 2.46e-10;
% theta_ = 1.1322795155315042;
if stack=="ABA"
    theta_ = 1.44;
elseif stack=="AAA"
    theta_ = 0.68;
end
theta=theta_/180*pi;
kd = 4*pi/(3*a0);
ktheta = 2*kd*sin(theta/2);

hbar = 1.054571817e-34;
J_to_meV = 6.24150636e21; %meV

kappa=0.6;
wab = 110;%mev
waa = wab*kappa;
% v0 = 5.581706932471464;
% v0 = 0.8481*10^6;%m/s;
v0 = 1*10^6;%m/s
v=v0*hbar*J_to_meV;%meV


%ktheta = wab/(alpha*v0)
am = a0/(2*sin(theta/2));

parameters.valley=-1;
parameters.theta_=theta_;


parameters.ktheta=ktheta;
parameters.v=v;
parameters.waa=waa;
parameters.wab=wab;%110.7;%mev

alpha = wab/(2*parameters.v*kd*sin(theta/2));



a1=2*pi/ktheta*[sqrt(3)/2,1/2]*(4/3);
a2=2*pi/ktheta*[-sqrt(3)/2,1/2]*(4/3);

if stack=="ABA"
    Phi=[0,1,-1]*(2*pi/3);
elseif stack=="AAA"
    Phi=[0,0,0]*(2*pi/3);
end

% AB_potential1=[40,0];
% AB_potential2=[40,0];
% AB_potential3=[40,0];
AB_potential1=[-4,4];
AB_potential2=[-6,6];
AB_potential3=[-8,8];
parameters.AB_potential=[AB_potential1;AB_potential2;AB_potential3];


L=40;
%nK=4;%number of g vectors in each side
E_keep=1;

[H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,[0,0],nK);

b1m=sqrt(3)*ktheta*[0.5, -sqrt(3)/2];
b2m=sqrt(3)*ktheta*[0.5, sqrt(3)/2];


k_gamma1=[0,0];
k_K=-q1;
k_Kp=-2*q1;
k_gamma2=k_Kp+q3;

PA=q2+q1;%K'
PB=q2;%K
PC=q2-q1;%Gamma1
PD=q2-q3;%Gamma2

M=(q2-q3)/2;

% 
% path1=k_line(PA,PB,L);
% path2=k_line(PB,PC,L);
% path3=k_line(PC,PD,2*L);
% path4=k_line(PD,PA,L);


path1=k_line(PB,PA,L);
path2=k_line(PA,M,L/2);
path3=k_line(M,PD,L);
path4=k_line(PD,PB,L);

path=[path1;path2(2:end,:);path3(2:end,:);path4(2:end,:)];


Es_set=zeros(E_keep*2,length(path));
Ea=[];
Eb=[];
for cc=1:length(path)
    [H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,path(cc,:),nK);
    dimm=length(eu);
    %Es_set(:,cc)=eu(dimm/2-E_keep+1:dimm/2+E_keep);
    Es_set(:,cc)=eu(dimm/2+E_keep:dimm/2+E_keep);
    plot(cc, Es_set(:,cc),'k.');hold on;
    siz=length(eu);
    Ea=[Ea,eu(siz/2)];
    Eb=[Eb,eu(siz/2+1)];
end
ylabel('$E$ (meV)','interpreter','latex')
xticks([1 L 1.5*L-1,2.5*L-2,3.5*L-3])
xticklabels({'K','K^{\prime}','M','\Gamma','K'})


%%%%%%
addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\HF_energy_2valley');
load e1theta1.44_epsilon4_waa0.6_L24.mat;

L0=sqrt(length(E_band1));
assert(abs(mod(L0,3))<1e-10);
E_total=E_band1+E_Hatree1+E_Hatree2*2;

PA=q2+q1;%K'
PB=q2;%K
PC=q2-q1;%Gamma1
PD=q2-q3;%Gamma2

M=(q2-q3)/2;


path1=k_line(PB,PA,L0/3+1);
path2=k_line(PA,M,L0/6+1);
path3=k_line(M,PD,L0/2+1);
path4=k_line(PD,PB,L0/3+1);

% path1_=linspace(1,L0/3+1,L0/3+1);
% path2_=linspace(length(path1_),length(path1_)+L0/6,L0/6+1);
% path3_=linspace(length(path1_)+length(path2_)-1,length(path1_)+length(path2_)-1+L0/2,L0/2+1);
% path4_=linspace(length(path1_)+length(path2_)+length(path3_)-2,length(path1_)+length(path2_)+length(path3_)-2+L0/3,L0/3+1);

L=40;
path1_=linspace(1,L,L0/3+1);
path2_=linspace(L, 1.5*L-1,L0/6+1);
path3_=linspace(1.5*L-1,2.5*L-2,L0/2+1);
path4_=linspace(2.5*L-2,3.5*L-3,L0/3+1);


k_set=k_set(:,1)*b1m+k_set(:,2)*b2m;

Es1=zeros(size(path1,1),1);
for cc=1:length(path1)
    kk=path1(cc,:);
    pos=find_pos(kk,k_set,b1m,b2m);
    Es1(cc)=E_total(pos);
end

Es2=zeros(size(path2,1),1);
for cc=1:length(path2)
    kk=path2(cc,:);
    pos=find_pos(kk,k_set,b1m,b2m);
    Es2(cc)=E_total(pos);
end

Es3=zeros(size(path3,1),1);
for cc=1:length(path3)
    kk=path3(cc,:);
    pos=find_pos(kk,k_set,b1m,b2m);
    Es3(cc)=E_total(pos);
end

Es4=zeros(size(path4,1),1);
for cc=1:length(path4)
    kk=path4(cc,:);
    pos=find_pos(kk,k_set,b1m,b2m);
    Es4(cc)=E_total(pos);
end


plot(path1_,Es1,'r.--');hold on;
plot(path2_,Es2,'r.--');hold on;
plot(path3_,Es3,'r.--');hold on;
plot(path4_,Es4,'r.--');hold on;
plot(linspace(0,1000,10),linspace(0,0,10),'k--');hold on;
%%%%%%


set(gca,'fontsize',12)
% xticks([0 12  24  36 ])
% xticklabels({'0','12','24','36'})
% title("FQAH at $\nu=1/3$",'interpreter','latex')

set(gca,'TickLength',[0.03, 0.01])
%rectangle('position',[1,2,3.5*L-2,24],'edgecolor','r','LineStyle','--');

axis([0,3.5*L,0,40])
%text(84,20,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 12);
%text(84,-20,'$C=1$',  'Interpreter', 'latex', 'fontsize', 12);
ar=annotation('doublearrow',[0.365 0.365],[0.795 0.915],'Color','r');
ar.Head1Width=6;
ar.Head1Length=5;
ar.Head2Width=ar.Head1Width;
ar.Head2Length=ar.Head1Length;


ar=annotation('doublearrow',[0.285 0.285],[0.63 0.8],'Color','k');
ar.Head1Width=6;
ar.Head1Length=5;
ar.Head2Width=ar.Head1Width;
ar.Head2Length=ar.Head1Length;

%text(-60,17,'renormalize',  'Interpreter', 'latex', 'fontsize', 12);

%title(['$w_{AA}/w_{AB}=0.6,$ ', '{\boldmath$\beta$}','$=(\sigma,K)$'],'interpreter','latex','fontsize',14)
%title(['$w_{AA}/w_{AB}=0.6,$ ','$(\sigma,K)$'],'interpreter','latex','fontsize',14)
title(['Single flavor ','$(\sigma,K)$'],'interpreter','latex','fontsize',14)

text(-40,45,'(a)',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
%set(gca, 'YTick', []);



% text(13,25,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 11,'Color','r');
% text(75,11,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 11);
% text(50,-15,'$C=1$',  'Interpreter', 'latex', 'fontsize', 11);
text(3.5*L-50,14,'$0<\nu<1$',  'Interpreter', 'latex', 'fontsize', 11);
text(5,37,'$\nu_{\rm{total}}=3+\nu$',  'Interpreter', 'latex', 'fontsize', 11,'Color','r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+0.5,pos_y0+0.48,length_x,length_y];
subplot('Position',pos);




plot(linspace(0,2,10),linspace(7,7,10),'Color','k');hold on;
plot(linspace(2.5,4.5,10),linspace(7,7,10),'Color','k');hold on;
plot(linspace(5,7,10),linspace(7,7,10),'Color','k');hold on;
plot(linspace(7.5,9.5,10),linspace(7,7,10),'Color','k');hold on;


plot(linspace(0,2,10),linspace(3,3,10),'Color','k');hold on;
plot(linspace(2.5,4.5,10),linspace(3,3,10),'Color','k');hold on;
plot(linspace(5,7,10),linspace(3,3,10),'Color','k');hold on;
plot(linspace(7.5,9.5,10),linspace(3,3,10),'Color','k');hold on;


text(0,8.5,'$0<\nu<1$',  'Interpreter', 'latex', 'fontsize', 12,'Color','k');





text(-0.5,1,'$(\sigma,K),$',  'Interpreter', 'latex', 'fontsize', 12,'Color','k');
text(2.2,1,'$(\bar{\sigma},K),$',  'Interpreter', 'latex', 'fontsize', 12,'Color','k');
text(4.9,1,'$(\sigma,K^{\prime}),$',  'Interpreter', 'latex', 'fontsize', 12,'Color','k');
text(7.7,1,'$(\bar{\sigma},K^{\prime})$',  'Interpreter', 'latex', 'fontsize', 12,'Color','k');

text(0,4.5,'$\nu_{\rm{total}}=3+\nu$',  'Interpreter', 'latex', 'fontsize', 12,'Color','r');

plot(0.4,7+0.2,'marker','.','Color','k','MarkerSize',20);hold on;
plot(0.4,3+0.2,'marker','.','Color','k','MarkerSize',20);hold on;

for x0=2.5:2.5:7.5
plot(x0+0.4,3+0.2,'marker','.','Color','k','MarkerSize',20);hold on;
plot(x0+0.4+0.6,3+0.2,'marker','.','Color','k','MarkerSize',20);hold on;
plot(x0+0.4+0.6*2,3+0.2,'marker','.','Color','k','MarkerSize',20);hold on;
end



text(-2.5,11,'(b)',  'Interpreter', 'latex', 'fontsize', 15);

axis([0,10,0,10]);
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################

clc;clear;

restoredefaultpath;
pos_x0=0.30;
pos_y0=0.45;
length_x=0.55;
length_y=0.05;
y_distance=0.065;
fig=gcf;
%%%%%%%%%%%%%%%%%%%
waa_set=0:0.05:0.8;


value_range=[0,5];

my_colormap1=hot(150);
my_colormap1=my_colormap1(end:-1:1,:);
my_colormap2=abyss(200);
my_colormap2=my_colormap2(1:end,:);
my_colormap=[my_colormap1;my_colormap2];

pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);
addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27')
load("Egap_Np18.mat");
Egap=treat_negative(Egap);
imagesc('XData',waa_set,'CData',Egap);
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,0.5,1]);
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
%text(-0.11,1.5,'$\nu=$',  'Interpreter', 'latex', 'fontsize', 17);
text(-0.4,0.7,'QHC, $\nu=2/3$',  'Interpreter', 'latex', 'fontsize', 12);

text(-0.42,1.2,'(c)',  'Interpreter', 'latex', 'fontsize', 15);



pos=[pos_x0,pos_y0-y_distance,length_x,length_y];
subplot('Position',pos);
addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N27_HF_2valley')
load("Egap_Np18_HF.mat");
Egap=treat_negative(Egap);
imagesc('XData',waa_set,'CData',Egap);
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,0.5,1]);
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
text(-0.4,0.7,'QHC, $\nu_{\rm{total}}=3+2/3$',  'Interpreter', 'latex', 'fontsize', 12);



pos=[pos_x0,pos_y0-y_distance*2,length_x,length_y];
subplot('Position',pos);
addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28')
load("Egap_Np14.mat");
Egap=treat_negative(Egap);
imagesc('XData',waa_set,'CData',Egap);
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,0.5,1]);
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
text(-0.4,0.7,'QHF, $\nu=1/2$',  'Interpreter', 'latex', 'fontsize', 12);


pos=[pos_x0,pos_y0-y_distance*3,length_x,length_y];
subplot('Position',pos);
addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley')
load("Egap_Np14_HF.mat");
Egap=treat_negative(Egap);
imagesc('XData',waa_set,'CData',Egap);
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,0.5,1]);
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
text(-0.4,0.7,'QHF, $\nu_{\rm{total}}=3+1/2$',  'Interpreter', 'latex', 'fontsize', 12);



pos=[pos_x0,pos_y0-y_distance*4,length_x,length_y];
subplot('Position',pos);
addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A')
load("Egap_Np8.mat");
Egap=treat_negative(Egap);
imagesc('XData',waa_set,'CData',Egap);
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,0.5,1]);
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
text(-0.4,0.7,'FCI, $\nu=1/3$',  'Interpreter', 'latex', 'fontsize', 12);



pos=[pos_x0,pos_y0-y_distance*5,length_x,length_y];
subplot('Position',pos);
addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N24A_HF_2valley')
load("Egap_Np8_HF.mat");
Egap=treat_negative(Egap);
imagesc('XData',waa_set,'CData',Egap);
caxis(value_range);
ax=gca;
colormap(gca,my_colormap);
axis([0,0.8,0.5,1]);
%set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8])
set(gca,'fontsize',12);
xlabel('$w_{AA}/w_{AB}$','interpreter','latex','fontsize',15)
text(-0.4,0.7,'FCI, $\nu_{\rm{total}}=3+1/3$',  'Interpreter', 'latex', 'fontsize', 12);



h = axes(fig,'visible','off'); 
% h.Title.Visible = 'on';
% h.XLabel.Visible = 'on';
% h.YLabel.Visible = 'on';
% ylabel(h,'yaxis','FontWeight','bold');
% xlabel(h,'xaxis','FontWeight','bold');
%title(h,'title');
cb=colorbar(h,'Position',[0.87 0.1 0.03 0.4]);
colormap(cb,my_colormap)
caxis(h,value_range); 
ylabel(cb,'Gap (meV)','Rotation',270,'interpreter','latex','fontsize',13)







set(gcf,'Position',[100 100 550 420])




exportgraphics(gcf,'band_property_new.eps','ContentType','vector')





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



function Egap=get_gap(epsilon,waa_set)
Kg_set=[6+1];
Ng_set=[1];
E_global=zeros(10,length(waa_set));

Eg_Ka=[];

Es_outside=[];

for cc=1:length(waa_set)
    filenm="two_band_theta1.44_waa"+num2str(waa_set(cc))+"_epsilon"+num2str(epsilon)+"_Np12.mat";
    load(filenm);
    Eg=min(sortcell(En_set));
    E_all=sortcell(En_set);
    E_global(:,cc)=E_all(1:size(E_global,1));


    Eg_Ka=[Eg_Ka,En_set{Kg_set(1)}(1:8)];

    E_outside=[];
    for dd=1:length(En_set)
        if length(find(Kg_set==dd))==0
            
            E_outside=[E_outside;En_set{dd}(1:8)];
        end
    end
    E_outside=sort(E_outside);
    Es_outside=[Es_outside,E_outside];

end


Eg=[Eg_Ka(1:1,:);];
Eg=min(Eg,[],1);

Eexcitation=min([Eg_Ka(2,:);Es_outside(1,:)],[],1);

Egap=Eexcitation-Eg;


end


function Ewidth=get_HF_bandwidth(epsilon,waa)
filenm=['e1theta1.44_epsilon',num2str(epsilon),'_waa',num2str(waa),'.mat'];
load(filenm);
L=sqrt(length(E_band1));


E_band1=reshape(E_band1,L,L);
E_Hatree1=reshape(E_Hatree1,L,L);
E_Hatree2=reshape(E_Hatree2,L,L);
E_=E_band1+E_Hatree1+E_Hatree2*2;
Ewidth=max(max(E_))-min(min(E_));

end


function pos=find_pos(kk0,k_set,b1m,b2m)
    for c1=-1:1
        for c2=-1:1
    
            kk=kk0+c1*b1m+c2*b2m;
            
            
            dis=(kk(1)-k_set(:,1)).^2+(kk(2)-k_set(:,2)).^2;
            pos=find(dis<1e-8);
            if ~isempty(pos)
                break;
            end
        end
        if ~isempty(pos)
            break;
        end
    end
    
    assert(length(pos)==1);
end


function ys=treat_negative(ys)
for cc=1:length(ys)
    if ys(cc)<0
        ys(cc)=0;
    end
end
end