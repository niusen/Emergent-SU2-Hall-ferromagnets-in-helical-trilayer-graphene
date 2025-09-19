clc;clear;close all;
%pathnm="D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\HF_energy_2valley\";
waa_set=0:0.05:0.8;

for cw=1:length(waa_set)
    waa=waa_set(cw);
    filenm="data\"+"e1theta1.44_epsilon4_waa"+num2str(waa);
    %filenm="e1theta1.44_epsilon4_waa"+num2str(waa);
    load(filenm);
    E_bare=E_band1;
    E_HF=E_band1+E_Hatree1+E_Hatree2*2;

    width_bare(cw)=max(real(E_bare(:)))-min(real(E_bare(:)));
    width_HF(cw)=max(real(E_HF(:)))-min(real(E_HF(:)));

end

f1=plot(waa_set,width_bare,'ro');hold on;
f2=plot(waa_set,width_HF,'bo');hold on;

plot(0.6*ones(1,20),linspace(0,30,20),'k--');hold on;


    % h=legend([f1;f2],'Bare band','$\nu_{\rm{tot}}=\nu+3$ with Hartree correction');
    h=legend([f1;f2],'Non-interacting $C_{band}=-2$ band','With Hartree correction');
    set(h,'Interpreter','latex','Location','northwest')
    set(h,'fontsize',13)
    set(h,'units','normalized');
    %set(h,'NumColumns',2);
    h.ItemTokenSize = [15,10];
    %h.NumColumns=2;
set(gca,'TickLength',[0.03, 0.01])
set(gca,'fontsize',12)

ylabel('E (meV)','Interpreter','latex','FontSize',15);
xlabel('$w_{AA}/w_{AB}$','interpreter','latex','FontSize',14)
title('Band width $W$','interpreter','latex','fontsize',15)



set(gcf,'Position',[30 100 500 400])

exportgraphics(gcf,'band_width.eps','ContentType','vector')
