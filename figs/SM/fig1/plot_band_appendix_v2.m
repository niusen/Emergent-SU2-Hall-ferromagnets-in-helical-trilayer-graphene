clc;clear;close all;
restoredefaultpath;
pos_x0=0.08;
pos_y0=0.12;
length_x=0.24;
length_y=0.34;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nK=2;%number of g vectors in each side
L=40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,pos_y0+0.48,length_x,length_y];
subplot('Position',pos);


%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band')
stack="ABA";%"ABA","AAA"
a0 = 2.46e-10;
% theta_ = 1.1322795155315042;
if stack=="ABA"
    theta_ = 1.5;
elseif stack=="AAA"
    theta_ = 0.68;
end
theta=theta_/180*pi;
kd = 4*pi/(3*a0);
ktheta = 2*kd*sin(theta/2);

hbar = 1.054571817e-34;
J_to_meV = 6.24150636e21; %meV

kappa=0;
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
AB_potential1=[-8,8];
AB_potential2=[-8,8];
AB_potential3=[-8,8];
parameters.AB_potential=[AB_potential1;AB_potential2;AB_potential3];




E_keep=2;

[H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,[0,0],nK);




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
path2=k_line(PA,M,L);
path3=k_line(M,PD,L);
path4=k_line(PD,PB,L);

path=[path1;path2(2:end,:);path3(2:end,:);path4(2:end,:)];


Es_set=zeros(E_keep*2,length(path));
Ea=[];
Eb=[];
for cc=1:length(path)
    [H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,path(cc,:),nK);
    dimm=length(eu);
    Es_set(:,cc)=eu(dimm/2-E_keep+1:dimm/2+E_keep);
    plot(cc, Es_set(:,cc),'k.');hold on;
    siz=length(eu);
    Ea=[Ea,eu(siz/2)];
    Eb=[Eb,eu(siz/2+1)];
end
ylabel('$E$ (meV)','interpreter','latex')
xticks([1 L-1 2*L-2,3*L-3,4*L-4])
xticklabels({'K','K^{\prime}','M','\Gamma','K'})

set(gca,'fontsize',12)
% xticks([0 12  24  36 ])
% xticklabels({'0','12','24','36'})
% title("FQAH at $\nu=1/3$",'interpreter','latex')
title("Spectrum",'interpreter','latex')
set(gca,'TickLength',[0.03, 0.01])
rectangle('position',[1,2,4*L-2,10],'edgecolor','r','LineStyle','--');

axis([0,4*L,-30,30])
text(50,20,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 12);
text(50,-20,'$C=1$',  'Interpreter', 'latex', 'fontsize', 12);

title(['$\theta=1.5, \kappa=0,$ \boldmath $\mu$ $=$ \boldmath $\mu_1$'],'interpreter','latex','fontsize',14)

text(-32,35,'(a)',  'Interpreter', 'latex', 'fontsize', 15,'color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+1/3,pos_y0+0.48,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band')
stack="ABA";%"ABA","AAA"
a0 = 2.46e-10;
% theta_ = 1.1322795155315042;
if stack=="ABA"
    theta_ = 1.5;
elseif stack=="AAA"
    theta_ = 0.68;
end
theta=theta_/180*pi;
kd = 4*pi/(3*a0);
ktheta = 2*kd*sin(theta/2);

hbar = 1.054571817e-34;
J_to_meV = 6.24150636e21; %meV

kappa=0;
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
AB_potential1=[0,0];
AB_potential2=[0,0];
AB_potential3=[-20,20];
parameters.AB_potential=[AB_potential1;AB_potential2;AB_potential3];



%nK=4;%number of g vectors in each side
E_keep=2;

[H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,[0,0],nK);




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
path2=k_line(PA,M,L);
path3=k_line(M,PD,L);
path4=k_line(PD,PB,L);

path=[path1;path2(2:end,:);path3(2:end,:);path4(2:end,:)];


Es_set=zeros(E_keep*2,length(path));
Ea=[];
Eb=[];
for cc=1:length(path)
    [H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,path(cc,:),nK);
    dimm=length(eu);
    Es_set(:,cc)=eu(dimm/2-E_keep+1:dimm/2+E_keep);
    plot(cc, Es_set(:,cc),'k.');hold on;
    siz=length(eu);
    Ea=[Ea,eu(siz/2)];
    Eb=[Eb,eu(siz/2+1)];
end
ylabel('$E$ (meV)','interpreter','latex')
xticks([1 L-1 2*L-2,3*L-3,4*L-4])
xticklabels({'K','K^{\prime}','M','\Gamma','K'})

set(gca,'fontsize',12)
% xticks([0 12  24  36 ])
% xticklabels({'0','12','24','36'})
% title("FQAH at $\nu=1/3$",'interpreter','latex')
title("Spectrum",'interpreter','latex')
set(gca,'TickLength',[0.03, 0.01])
rectangle('position',[1,2,4*L-2,10],'edgecolor','r','LineStyle','--');

axis([0,4*L,-30,30])
text(105,20,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 12);
text(105,-20,'$C=1$',  'Interpreter', 'latex', 'fontsize', 12);

title(['$\theta=1.5, \kappa=0,$ \boldmath $\mu$ $=$ \boldmath $\mu_2$'],'interpreter','latex','fontsize',14)

text(-32,35,'(b)',  'Interpreter', 'latex', 'fontsize', 15,'color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+2/3,pos_y0+0.48,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band')
stack="ABA";%"ABA","AAA"
a0 = 2.46e-10;
% theta_ = 1.1322795155315042;
if stack=="ABA"
    theta_ = 1.5;
elseif stack=="AAA"
    theta_ = 0.68;
end
theta=theta_/180*pi;
kd = 4*pi/(3*a0);
ktheta = 2*kd*sin(theta/2);

hbar = 1.054571817e-34;
J_to_meV = 6.24150636e21; %meV

kappa=0;
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
AB_potential1=[-8,8];
AB_potential2=[-0,0];
AB_potential3=[-8,8];
parameters.AB_potential=[AB_potential1;AB_potential2;AB_potential3];



%nK=4;%number of g vectors in each side
E_keep=2;

[H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,[0,0],nK);




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
path2=k_line(PA,M,L);
path3=k_line(M,PD,L);
path4=k_line(PD,PB,L);

path=[path1;path2(2:end,:);path3(2:end,:);path4(2:end,:)];


Es_set=zeros(E_keep*2,length(path));
Ea=[];
Eb=[];
for cc=1:length(path)
    [H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,path(cc,:),nK);
    dimm=length(eu);
    Es_set(:,cc)=eu(dimm/2-E_keep+1:dimm/2+E_keep);
    plot(cc, Es_set(:,cc),'k.');hold on;
    siz=length(eu);
    Ea=[Ea,eu(siz/2)];
    Eb=[Eb,eu(siz/2+1)];
end
ylabel('$E$ (meV)','interpreter','latex')
xticks([1 L-1 2*L-2,3*L-3,4*L-4])
xticklabels({'K','K^{\prime}','M','\Gamma','K'})

set(gca,'fontsize',12)
% xticks([0 12  24  36 ])
% xticklabels({'0','12','24','36'})
% title("FQAH at $\nu=1/3$",'interpreter','latex')
title("Spectrum",'interpreter','latex')
set(gca,'TickLength',[0.03, 0.01])
rectangle('position',[1,2,4*L-2,10],'edgecolor','r','LineStyle','--');

axis([0,4*L,-30,30])
text(105,20,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 12);
text(105,-20,'$C=1$',  'Interpreter', 'latex', 'fontsize', 12);

title(['$\theta=1.5, \kappa=0,$ \boldmath $\mu$ $=$ \boldmath $\mu_3$'],'interpreter','latex','fontsize',14)

text(-32,35,'(c)',  'Interpreter', 'latex', 'fontsize', 15,'color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band')
stack="ABA";%"ABA","AAA"
a0 = 2.46e-10;
% theta_ = 1.1322795155315042;
if stack=="ABA"
    theta_ = 1.5;
elseif stack=="AAA"
    theta_ = 0.68;
end
theta=theta_/180*pi;
kd = 4*pi/(3*a0);
ktheta = 2*kd*sin(theta/2);

hbar = 1.054571817e-34;
J_to_meV = 6.24150636e21; %meV

kappa=0;
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



%nK=4;%number of g vectors in each side
E_keep=2;

[H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,[0,0],nK);




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
path2=k_line(PA,M,L);
path3=k_line(M,PD,L);
path4=k_line(PD,PB,L);

path=[path1;path2(2:end,:);path3(2:end,:);path4(2:end,:)];


Es_set=zeros(E_keep*2,length(path));
Ea=[];
Eb=[];
for cc=1:length(path)
    [H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,path(cc,:),nK);
    dimm=length(eu);
    Es_set(:,cc)=eu(dimm/2-E_keep+1:dimm/2+E_keep);
    plot(cc, Es_set(:,cc),'k.');hold on;
    siz=length(eu);
    Ea=[Ea,eu(siz/2)];
    Eb=[Eb,eu(siz/2+1)];
end
ylabel('$E$ (meV)','interpreter','latex')
xticks([1 L-1 2*L-2,3*L-3,4*L-4])
xticklabels({'K','K^{\prime}','M','\Gamma','K'})

set(gca,'fontsize',12)
% xticks([0 12  24  36 ])
% xticklabels({'0','12','24','36'})
% title("FQAH at $\nu=1/3$",'interpreter','latex')
title("Spectrum",'interpreter','latex')
set(gca,'TickLength',[0.03, 0.01])
rectangle('position',[1,2,4*L-2,10],'edgecolor','r','LineStyle','--');

axis([0,4*L,-30,30])
text(105,20,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 12);
text(105,-20,'$C=1$',  'Interpreter', 'latex', 'fontsize', 12);

title(['$\theta=1.5, \kappa=0,$ \boldmath $\mu$ $=$ \boldmath $\mu_4$'],'interpreter','latex','fontsize',14)

text(-32,35,'(d)',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+1/3,pos_y0,length_x,length_y];
subplot('Position',pos);


%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band')
stack="ABA";%"ABA","AAA"
a0 = 2.46e-10;
% theta_ = 1.1322795155315042;
if stack=="ABA"
    theta_ = 1.5;
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




E_keep=2;

[H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,[0,0],nK);




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
path2=k_line(PA,M,L);
path3=k_line(M,PD,L);
path4=k_line(PD,PB,L);

path=[path1;path2(2:end,:);path3(2:end,:);path4(2:end,:)];


Es_set=zeros(E_keep*2,length(path));
Ea=[];
Eb=[];
for cc=1:length(path)
    [H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,path(cc,:),nK);
    dimm=length(eu);
    Es_set(:,cc)=eu(dimm/2-E_keep+1:dimm/2+E_keep);
    plot(cc, Es_set(:,cc),'k.');hold on;
    siz=length(eu);
    Ea=[Ea,eu(siz/2)];
    Eb=[Eb,eu(siz/2+1)];
end
max(Es_set(3,:))-min(Es_set(3,:))
ylabel('$E$ (meV)','interpreter','latex')
xticks([1 L-1 2*L-2,3*L-3,4*L-4])
xticklabels({'K','K^{\prime}','M','\Gamma','K'})

set(gca,'fontsize',12)
% xticks([0 12  24  36 ])
% xticklabels({'0','12','24','36'})
% title("FQAH at $\nu=1/3$",'interpreter','latex')
title("Spectrum",'interpreter','latex')
set(gca,'TickLength',[0.03, 0.01])
rectangle('position',[1,2,4*L-2,24],'edgecolor','r','LineStyle','--');

axis([0,4*L,-30,30])
text(105,20,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 12);
text(105,-20,'$C=1$',  'Interpreter', 'latex', 'fontsize', 12);

title(['$\theta=1.5, \kappa=0.6,$ \boldmath $\mu$ $=$ \boldmath $\mu_4$'],'interpreter','latex','fontsize',14)

text(-32,35,'(e)',  'Interpreter', 'latex', 'fontsize', 15,'color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+2/3,pos_y0,length_x,length_y];
subplot('Position',pos);


%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band')
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



%nK=4;%number of g vectors in each side
E_keep=2;

[H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,[0,0],nK);




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
path2=k_line(PA,M,L);
path3=k_line(M,PD,L);
path4=k_line(PD,PB,L);

path=[path1;path2(2:end,:);path3(2:end,:);path4(2:end,:)];


Es_set=zeros(E_keep*2,length(path));
Ea=[];
Eb=[];
for cc=1:length(path)
    [H,eu,ev,q1,q2,q3,basis]=band_single_k(parameters,Phi,path(cc,:),nK);
    dimm=length(eu);
    Es_set(:,cc)=eu(dimm/2-E_keep+1:dimm/2+E_keep);
    plot(cc, Es_set(:,cc),'k.');hold on;
    siz=length(eu);
    Ea=[Ea,eu(siz/2)];
    Eb=[Eb,eu(siz/2+1)];
end
max(Es_set(3,:))-min(Es_set(3,:))
ylabel('$E$ (meV)','interpreter','latex')
xticks([1 L-1 2*L-2,3*L-3,4*L-4])
xticklabels({'K','K^{\prime}','M','\Gamma','K'})

set(gca,'fontsize',12)
% xticks([0 12  24  36 ])
% xticklabels({'0','12','24','36'})
% title("FQAH at $\nu=1/3$",'interpreter','latex')
title("Spectrum",'interpreter','latex')
set(gca,'TickLength',[0.03, 0.01])
rectangle('position',[1,2,4*L-2,24],'edgecolor','r','LineStyle','--');

axis([0,4*L,-30,30])
text(105,20,'$C=-2$',  'Interpreter', 'latex', 'fontsize', 12);
text(105,-20,'$C=1$',  'Interpreter', 'latex', 'fontsize', 12);

title(['$\theta=1.44, \kappa=0.6,$ \boldmath $\mu$ $=$ \boldmath $\mu_4$'],'interpreter','latex','fontsize',14)

text(-32,35,'(f)',  'Interpreter', 'latex', 'fontsize', 15,'color','k');






set(gcf,'Position',[100 100 850 420])




exportgraphics(gcf,'band_appendix.eps','ContentType','vector')





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