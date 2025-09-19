clc;clear;close all;
restoredefaultpath
subplot(1,2,1)
N=28;
En_set=cell(N);
filenm0="_e1theta1.44_epsilon4_waa0.1_Np14_flx_0_0.mat";
for cc=1:N
filenm="Kind"+num2str(cc)+filenm0;
load(filenm);
En_set{cc}=eu;
end
Eg=min(sortcell(En_set));
E_all=sortcell(En_set);
for c1=1:size(En_set,1)
    En=En_set{c1};
    Es=En-Eg;
    %Es=En;
    plot(c1-1,Es(1:6),'k.','markersize',10);hold on;
end

