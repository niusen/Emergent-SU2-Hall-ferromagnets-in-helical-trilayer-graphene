clc;clear;close all;
restoredefaultpath;
pos_x0=0.06;
pos_y0=0.13;
length_x=0.18;
length_y=0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,pos_y0+0.5,length_x,length_y];
subplot('Position',pos);


Ns=[16,20,24,28,32];




%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N16_HF_2valley\order')
load("data/Ns16/Sxyz_waa0.6_option2.mat");
Sx_16=real(eig(Opx))*2;

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N20B_HF_2valley\order')
load("data/Ns20/Sxyz_waa0.6_option2.mat");
Sx_20=real(eig(Opx))*2;

load("data/Ns24/Sxyz_waa0.6_option2.mat");
Sx_24=real(eig(Opx))*2;

load("data/Ns28/Sxyz_waa0.6_option2.mat");
Sx_28=real(eig(Opx))*2;

load("data/Ns32/Sxyz_waa0.6_option2.mat");
Sx_32=real(eig(Opx))*2;



markerlength=30;
Markerwidth=0.6;
LW=2;


component_max=[];


N=16;
Sx=Sx_16;
for cc=1:length(Sx)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sx(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sx))];

N=20;
Sx=Sx_20;
for cc=1:length(Sx)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sx(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sx))];

N=24;
Sx=Sx_24;
for cc=1:length(Sx)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sx(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sx))];

N=28;
Sx=Sx_28;
for cc=1:length(Sx)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sx(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sx))];

N=32;
Sx=Sx_32;
for cc=1:length(Sx)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sx(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sx))];

set(gca,'fontsize',12)
axis([15,35,-8.5,8.5])
%yticks([-8,-4,0,4,8])


yticks([-8,-6,-4,-2,0,2,4,6,8])
yticklabels({'-8','','-4','','0','','4','','8'})

xticks([16,20,24,28,32])
xlabel('$N_s$','interpreter','latex','FontSize',12)
%ylabel('$S^z$','interpreter','latex','FontSize',12)

text(10,11,'(a)','interpreter','latex','fontsize',14)

title("$\hat{S}^x$ eigenvalues",'interpreter','latex','FontSize',14)
set(gca,'TickLength',[0.03, 0.01])

% range=2:4;
% xs=1./Nset(range);
% ys=Ees(range);
% p = polyfit(xs,ys,1);
% xx=0:0.001:1/15;
% yy=polyval(p,xx);
% plot(xx,yy,'b--');hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);

plot(1./Ns,(Ns/2-component_max)./(Ns/2),'ro'); hold on;

xs=1./Ns;
ys=(Ns/2-component_max)./(Ns/2);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'r--');hold on;

axis([0,1/15,0,0.04])

xlabel('$1/N_s$','interpreter','latex','FontSize',12)
set(gca,'TickLength',[0.03, 0.01])
set(gca,'fontsize',12)
xticks([0.02,0.04,0.06])
yticks([0,0.01,0.02,0.03,0.04])
title("Relative error",'interpreter','latex','FontSize',14)
text(-0.02,0.048,'(e)','interpreter','latex','fontsize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+1/4,pos_y0+0.5,length_x,length_y];
subplot('Position',pos);

Ns=[16,20,24,28,32];

load("data/Ns16/Sxyz_waa0.6_option2.mat");
Sy_16=real(eig(Opy))*2;

load("data/Ns20/Sxyz_waa0.6_option2.mat");
Sy_20=real(eig(Opy))*2;

load("data/Ns24/Sxyz_waa0.6_option2.mat");
Sy_24=real(eig(Opy))*2;

load("data/Ns28/Sxyz_waa0.6_option2.mat");
Sy_28=real(eig(Opy))*2;

load("data/Ns32/Sxyz_waa0.6_option2.mat");
Sy_32=real(eig(Opy))*2;




markerlength=30;
Markerwidth=0.6;
LW=2;

component_max=[];

N=16;
Sy=Sy_16;
for cc=1:length(Sy)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sy(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sy))];

N=20;
Sy=Sy_20;
for cc=1:length(Sy)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sy(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sy))];

N=24;
Sy=Sy_24;
for cc=1:length(Sy)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sy(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sy))];

N=28;
Sy=Sy_28;
for cc=1:length(Sy)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sy(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sy))];

N=32;
Sy=Sy_32;
for cc=1:length(Sy)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sy(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sy))];

set(gca,'fontsize',12)
axis([15,35,-8.5,8.5])
%yticks([-8,-4,0,4,8])


yticks([-8,-6,-4,-2,0,2,4,6,8])
yticklabels({'-8','','-4','','0','','4','','8'})

xticks([16,20,24,28,32])
xlabel('$N_s$','interpreter','latex','FontSize',12)
%ylabel('$S^z$','interpreter','latex','FontSize',12)

text(10,11,'(b)','interpreter','latex','fontsize',14)

title("$\hat{S}^y$ eigenvalues",'interpreter','latex','FontSize',14)
set(gca,'TickLength',[0.03, 0.01])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+1/4,pos_y0,length_x,length_y];
subplot('Position',pos);

plot(1./Ns,(Ns/2-component_max)./(Ns/2),'ro'); hold on;



xs=1./Ns;
ys=(Ns/2-component_max)./(Ns/2);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'r--');hold on;

axis([0,1/15,0,0.04])

xlabel('$1/N_s$','interpreter','latex','FontSize',12)
set(gca,'TickLength',[0.03, 0.01])
set(gca,'fontsize',12)
xticks([0.02,0.04,0.06])
yticks([0,0.01,0.02,0.03,0.04])
title("Relative error",'interpreter','latex','FontSize',14)
text(-0.02,0.048,'(f)','interpreter','latex','fontsize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+2/4,pos_y0+0.5,length_x,length_y];
subplot('Position',pos);

Ns=[16,20,24,28,32];

load("data/Ns16/Sxyz_waa0.6_option2.mat");
Sz_16=real(eig(Opz))*2;

load("data/Ns20/Sxyz_waa0.6_option2.mat");
Sz_20=real(eig(Opz))*2;

load("data/Ns24/Sxyz_waa0.6_option2.mat");
Sz_24=real(eig(Opz))*2;

load("data/Ns28/Sxyz_waa0.6_option2.mat");
Sz_28=real(eig(Opz))*2;

load("data/Ns32/Sxyz_waa0.6_option2.mat");
Sz_32=real(eig(Opz))*2;



markerlength=30;
Markerwidth=0.6;
LW=2;

component_max=[];

N=16;
Sz=Sz_16;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sz))];

N=20;
Sz=Sz_20;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sz))];

N=24;
Sz=Sz_24;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sz))];

N=28;
Sz=Sz_28;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sz))];

N=32;
Sz=Sz_32;
for cc=1:length(Sz)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Sz(cc))/2*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Sz))];

set(gca,'fontsize',12)
axis([15,35,-8.5,8.5])
%yticks([-8,-4,0,4,8])


yticks([-8,-6,-4,-2,0,2,4,6,8])
yticklabels({'-8','','-4','','0','','4','','8'})

xticks([16,20,24,28,32])
xlabel('$N_s$','interpreter','latex','FontSize',12)
%ylabel('$S^z$','interpreter','latex','FontSize',12)

text(10,11,'(c)','interpreter','latex','fontsize',14)

title("$\hat{S}^z$ eigenvalues",'interpreter','latex','FontSize',14)
set(gca,'TickLength',[0.03, 0.01])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+2/4,pos_y0,length_x,length_y];
subplot('Position',pos);

plot(1./Ns,(Ns/2-component_max)./(Ns/2),'ro'); hold on;

xs=1./Ns;
ys=(Ns/2-component_max)./(Ns/2);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'r--');hold on;

axis([0,1/15,0,0.04])

xlabel('$1/N_s$','interpreter','latex','FontSize',12)
set(gca,'TickLength',[0.03, 0.01])
set(gca,'fontsize',12)
xticks([0.02,0.04,0.06])
yticks([0,0.01,0.02,0.03,0.04])
title("Relative error",'interpreter','latex','FontSize',14)
text(-0.02,0.048,'(g)','interpreter','latex','fontsize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+3/4,pos_y0+0.5,length_x,length_y];
subplot('Position',pos);

Ns=[16,20,24,28,32];

load("data/Ns16/Sxyz_waa0.6_option2.mat");
Stotal_16=real(eig(Opz*Opz+Opx*Opx+Opy*Opy));

load("data/Ns20/Sxyz_waa0.6_option2.mat");
Stotal_20=real(eig(Opz*Opz+Opx*Opx+Opy*Opy));

load("data/Ns24/Sxyz_waa0.6_option2.mat");
Stotal_24=real(eig(Opz*Opz+Opx*Opx+Opy*Opy));

load("data/Ns28/Sxyz_waa0.6_option2.mat");
Stotal_28=real(eig(Opz*Opz+Opx*Opx+Opy*Opy));

load("data/Ns32/Sxyz_waa0.6_option2.mat");
Stotal_32=real(eig(Opz*Opz+Opx*Opx+Opy*Opy));



markerlength=30;
Markerwidth=0.6;
LW=2;

component_max=[];

N=16;
Stotal=Stotal_16;
Stotal=(sqrt(1+Stotal*4)-1)/2;
for cc=1:length(Stotal)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Stotal(cc))*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Stotal))];

N=20;
Stotal=Stotal_20;
Stotal=(sqrt(1+Stotal*4)-1)/2;
for cc=1:length(Stotal)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Stotal(cc))*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Stotal))];

N=24;
Stotal=Stotal_24;
Stotal=(sqrt(1+Stotal*4)-1)/2;
for cc=1:length(Stotal)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Stotal(cc))*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Stotal))];

N=28;
Stotal=Stotal_28;
Stotal=(sqrt(1+Stotal*4)-1)/2;
for cc=1:length(Stotal)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Stotal(cc))*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Stotal))];

N=32;
Stotal=Stotal_32;
Stotal=(sqrt(1+Stotal*4)-1)/2;
for cc=1:length(Stotal)
    f1=plot(linspace(N-Markerwidth,N+Markerwidth,markerlength),real(Stotal(cc))*ones(1,markerlength),'r-','linewidth',LW);hold on;
end
component_max=[component_max,max(real(Stotal))];

set(gca,'fontsize',12)
axis([15,35,0,8.5])
%yticks([-8,-4,0,4,8])


yticks([-8,-6,-4,-2,0,2,4,6,8])
yticklabels({'-8','','-4','','0','','4','','8'})

xticks([16,20,24,28,32])
xlabel('$N_s$','interpreter','latex','FontSize',12)
%ylabel('$S^z$','interpreter','latex','FontSize',12)

text(10,9.6,'(d)','interpreter','latex','fontsize',14)

title("Total $S$",'interpreter','latex','FontSize',14)
set(gca,'TickLength',[0.03, 0.01])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=[pos_x0+3/4,pos_y0,length_x,length_y];
subplot('Position',pos);

plot(1./Ns,(Ns/4-component_max)./(Ns/4),'ro'); hold on;

xs=1./Ns;
ys=(Ns/4-component_max)./(Ns/4);
p = polyfit(xs,ys,1);
xx=0:0.001:1/15;
yy=polyval(p,xx);
plot(xx,yy,'r--');hold on;

axis([0,1/15,0,0.04])

xlabel('$1/N_s$','interpreter','latex','FontSize',12)
set(gca,'TickLength',[0.03, 0.01])
set(gca,'fontsize',12)
xticks([0.02,0.04,0.06])
yticks([0,0.01,0.02,0.03,0.04])
title("Relative error",'interpreter','latex','FontSize',14)
text(-0.02,0.048,'(h)','interpreter','latex','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%
set(gcf,'Position',[100 100 850 380])




exportgraphics(gcf,'scaling_quantum_number.eps','ContentType','vector')



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