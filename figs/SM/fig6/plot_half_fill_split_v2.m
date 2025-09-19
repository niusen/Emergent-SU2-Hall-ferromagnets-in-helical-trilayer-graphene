clc;clear;close all;
restoredefaultpath
pos_x0=0.12;
pos_y0=0.23;
length_x=0.35;
length_y=0.6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);

%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\HF_energy_2valley');





waa_set=0:0.05:0.9;
for cw=1:length(waa_set)
    waa=waa_set(cw);
    filenm="data\e1theta1.44_epsilon4_waa"+num2str(waa)+"_N28";
    %filenm="e1theta1.44_epsilon4_waa"+num2str(waa);
    load(filenm);
    E_bare=E_band1;
    E_HF=E_band1+E_Hatree1+E_Hatree2*2;

    gaps_M=energy_connected_by_Mpoints(k_set,E_HF);
    width_HF_M(cw)=max(gaps_M);

    width_bare(cw)=max(real(E_bare(:)))-min(real(E_bare(:)));
    width_HF(cw)=max(real(E_HF(:)))-min(real(E_HF(:)));

end

f1=plot(waa_set,width_HF,'r.-');hold on;
%f12=plot(waa_set,width_HF_M,'m.-');hold on;
f2=plot(waa_set,V_typical*ones(1,length(waa_set)),'b.-');hold on;

set(gca,'fontsize',12);
axis([0,0.8,-0.3,40])
xticks([0 0.2,0.4,0.6,0.9])
xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);

title(['(a)'],'interpreter','latex','fontsize',15)

ylabel('E (meV)','Interpreter','latex','FontSize',15);

text(0.1,30,'Interaction $U$',  'Interpreter', 'latex', 'fontsize', 12,'color','b');
text(0.2,3,'Band width $W$',  'Interpreter', 'latex', 'fontsize', 12,'color','r');
%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0+1/2,pos_y0,length_x,length_y];
subplot('Position',pos);




%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\order')
waa_set=[0,0.2,0.4,0.6,0.8,0.9];
for cw=1:length(waa_set)
    waa=waa_set(cw);
    filenm="kinetics_theta1.44_waa"+num2str(waa)+"_epsilon4_Np14";
    %filenm="e1theta1.44_epsilon4_waa"+num2str(waa);
    load("data/"+filenm);
    eks=[];
    for ck=1:length(kinetics_all)
        kinetics=kinetics_all{ck};
        for cs=1:length(kinetics)
            eks=[eks;kinetics{cs}];
        end
    end
    eks=real(eks);
    kinetic_diff(cw)=max(eks)-min(eks);
end

f1=plot(waa_set,kinetic_diff,'k.-','LineWidth',1.5);hold on;



%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\spectrum_scan_waa\');
Kg_set=[0+1,1+1,11+1,15+1];
Ng_set=[3,4,4,4];
E_global=zeros(10,length(waa_set));
waa_set=0:0.05:0.9;
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
    filenm="spectrum_scan_waa\"+"Kind"+num2str(dd)+filenm0;
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
Eg_max=max([Eg_Ka(3,:);Eg_Kb(4,:);Eg_Kc(4,:);Eg_Kd(4,:);],[],1);


f1=plot(waa_set,Eg_max-Eg,'r.-','LineWidth',1.5);hold on;
f2=plot(waa_set,Eexcitation-Eg_max,'b.-','LineWidth',1.5);hold on;



ylabel('E (meV)','Interpreter','latex','FontSize',15);
    
set(gca,'fontsize',12);
axis([0,0.9,-0,8])
xticks([0 0.2,0.4,0.6,0.8])
xlabel('$w_{AA}/w_{AB}$','Interpreter','latex','FontSize',15);

title(['(b)'],'interpreter','latex','fontsize',15)
%text(0.2,2,'QHF',  'Interpreter', 'latex', 'fontsize', 15,'color','k');
text(0.5,6.5,'Gap $\Delta_{gap}$',  'Interpreter', 'latex', 'fontsize', 12,'color','b');
text(0.05,1,'$\Delta_{E_g}$, ',  'Interpreter', 'latex', 'fontsize', 14,'color','r');
text(0.25,1,'$\Delta_{E_k}$',  'Interpreter', 'latex', 'fontsize', 14,'color','k');
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


set(gcf,'Position',[100 100 520 250])



exportgraphics(gcf,'split_half_fill.eps','ContentType','vector')





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

function partners=find_partner(k_set,kk,Q)



        kk=kk+Q;
        dk1=kk(1)-k_set(:,1);
        dk2=kk(2)-k_set(:,2);

        cond1=(abs(round(dk1)-dk1)<1e-10);
        cond2=(abs(round(dk2)-dk2)<1e-10);

        pos=find((cond1.*cond2)==1);
        assert(length(pos)==1);
        partners=pos;
   
end

function dEs=energy_connected_by_Mpoints(kset,Es)
   Q1=[0,0.5];
   Q2=[0.5,0];
   Q3=[0.5,0.5];
   dEs=zeros(size(kset,1),1);
   for cc=1:size(kset,1)
       kk=kset(cc,:);
    p1=find_partner(kset,kk,Q1);
    p2=find_partner(kset,kk,Q2);
    p3=find_partner(kset,kk,Q3);
    eset=[Es(cc),Es(p1),Es(p2),Es(p3)];
    dEs(cc)=max(eset)-min(eset);
   end
end
