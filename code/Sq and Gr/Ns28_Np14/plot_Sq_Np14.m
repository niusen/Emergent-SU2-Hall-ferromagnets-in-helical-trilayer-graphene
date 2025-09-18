clc;clear;close all;
load("D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\arXiv2405_08887\data\flatter\mymap.mat")
K_id_set=[0+1,1+1,11+1,15+1];%
n_keep_set=[3,4,4,4];%






N=28;
Np=14;

waa=0.6;

filenm="FBZ_paircorrel_theta1.44_waa"+num2str(waa,'%.1f')+"_epsilon4_Np14_sec";
% filenm="FBZ_paircorrel_theta1.44_waa"+"0.0"+"_epsilon4_Np14_sec";
load(filenm+num2str(K_id_set(1))+".mat");
Ng=sum(n_keep_set);



Sq_all=zeros(length(q_FBZ),Ng,2*nK+1,2*nK+1);
Sq_all_subtract=Sq_all;
step=0;
for ck=1:length(K_id_set)
    load(filenm+num2str(K_id_set(ck))+".mat");
    %%%%%%%%%%%%
    
    pos=find((q_FBZ(:,1).^2+q_FBZ(:,2).^2)<1e-10);
    [l1,l3]=size(rho_q1);
    rho_q1=reshape(rho_q1,[2*nK+1,2*nK+1,l3]);
    rho_q2=reshape(rho_q2,[2*nK+1,2*nK+1,l3]);
    rho_q1=permute(rho_q1,[3,1,2]);
    rho_q2=permute(rho_q2,[3,1,2]);
    rhoqrhoq=reshape(rho_q1.*rho_q2,[1,l3,2*nK+1,2*nK+1]);

    %Sq_set(pos,:,nK+1,nK+1)=Sq_set(pos,:,nK+1,nK+1)-rhoqrhoq(1,:,nK+1,nK+1);
    %assert(norm(Sq_set(:)-abs(Sq_set(:)))<1e-4);
    %Sq_set_subtract=abs(Sq_set_subtract);
    %%%%%%%%%%%%
    for cs=1:n_keep_set(ck)
        step=step+1;



        Sq_all(:,step,:,:)=Sq_set(:,cs,:,:);
        Sq_all_subtract(:,step,:,:)=Sq_set(:,cs,:,:);
        Sq_all_subtract(pos,step,nK+1,nK+1)=Sq_all_subtract(pos,step,nK+1,nK+1)-rhoqrhoq(1,cs,nK+1,nK+1);


    end
end




nK=double(nK);
Sq_average=real(sum(Sq_all,2))/Ng;
Sq_average_subtract=real(sum(Sq_all_subtract,2))/Ng;

[corners,in_FBZ,qx_set,qy_set]=determine_in_FBZ(k_set,q_FBZ);


for cc=1:6
    p1=corners(cc,:);
    p2=corners(cc+1,:);
    plot([p1(1),p2(1)],[p1(2),p2(2)],'k-','LineWidth',6);hold on;
    
end

%define momentum space and real space lattice vectors
g1=[1,0,0];
g3=[-1/2,sqrt(3)/2,0];
z_3d=[0,0,1];

a1=2*pi*cross(z_3d,g3)/(dot(g1, cross(z_3d,g3)));
a2=2*pi*cross(z_3d,g1)/(dot(g3, cross(z_3d,g1)));
a1=a1(1:2);
a2=a2(1:2);


K_set=zeros(size(q_FBZ,1),2,2*nK+1,2*nK+1);
R_set=zeros(size(q_FBZ,1),2,2*nK+1,2*nK+1);
gr_set=zeros(size(q_FBZ,1),1,2*nK+1,2*nK+1);
for ca=-nK:nK
    for cb=-nK:nK
        for c1=1:size(q_FBZ,1)
            kk=q_FBZ(c1,1)*g1+q_FBZ(c1,2)*g3+ca*g1+cb*g3;
            K_set(c1,:,ca+nK+1,cb+nK+1)=kk(1:2);
            Ri=2*pi*cross(z_3d,kk)/(dot(g3, cross(z_3d,g1)));
            R_set(c1,:,ca+nK+1,cb+nK+1)=Ri(1:2);
        
        end
    end
end







value_range=[min(real(Sq_all_subtract(:))),max(real(Sq_all_subtract(:)))];


for ca=-2:2
    for cb=-2:2
        for cq=1:size(q_FBZ,1)
            
            cq
            Sq_average_subtract(cq,1,ca+nK+1,cb+nK+1)
            coord=K_set(cq,:,ca+nK+1,cb+nK+1);
            if norm(coord)<1e-6
                Color=determine_color(value_range,0);
            else
                Color=determine_color(value_range,real(Sq_average_subtract(cq,1,ca+nK+1,cb+nK+1)));
            end
            plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',25);hold on;%colorbar;
        
        end
    end
end

%%%%%%%%%%%%%%%%
figure
for cc=1:6
    p1=corners(cc,:);
    p2=corners(cc+1,:);
    plot([p1(1),p2(1)],[p1(2),p2(2)],'k-','LineWidth',6);hold on;
    
end

Sq_points=[];
Sq_values=[];
for ca=-1:1
    for cb=-1:1
        for c1=1:size(q_FBZ,1)
            kk=q_FBZ(c1,:)+[ca,cb];

            in_FBZ=determine_in_FBZ_single(k_set,kk);
            if in_FBZ==1
                Sq_points=[Sq_points;kk];
                Sq_values=[Sq_values,Sq_average_subtract(c1,1,ca+nK+1,cb+nK+1)];
            end
            
        end
    end
end


value_range=[min(real(Sq_values(:))),max(real(Sq_values(:)))];
for cq=1:length(Sq_values)
    
    cq
    Sq_values(cq)
    coord=Sq_points(cq,1)*g1+Sq_points(cq,2)*g3;

    Color=determine_color(value_range,real(Sq_values(cq)));
    
    plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',25);hold on;%colorbar;

end

caxis(value_range);
ax=gca;
colormap(gca,mycmap);
cb=colorbar;

save("Sq_28_14"+"_waa"+num2str(waa)+".mat","corners","Sq_values","Sq_points","g1","g3");


%%%%%%%%%%%%%%%%%%

for ca=-nK:nK
    for cb=-nK:nK
        for c1=1:size(q_FBZ,1)
            coord=R_set(c1,:,ca+nK+1,cb+nK+1);
            gr=get_gr(Sq_average,K_set,coord)*N/Np/Np;
            
            gr_set(c1,1,ca+nK+1,cb+nK+1)=gr;


        end
    end
end
% 
% figure
% gr_all=real(gr_set(:));
% 
% value_range=[min(gr_all),max(gr_all)];
% 
% for ca=-nK:nK
%     for cb=-nK:nK
%         for cq=1:size(q_FBZ,1)
% 
%             cq
%             gr_set(cq,1,ca+nK+1,cb+nK+1)
%             coord=R_set(cq,:,ca+nK+1,cb+nK+1);
% 
%             Color=determine_color(value_range,real(gr_set(cq,1,ca+nK+1,cb+nK+1)));
%             plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',10);hold on;%colorbar;
% 
%         end
%     end
% end
% Rx=R_set(:,1,:,:);
% Ry=R_set(:,2,:,:);
% 
% 
% caxis(value_range);
% ax=gca;
% colormap(gca,mycmap);
% cb=colorbar;
% 
% set(gcf,'Position',[30 100 600 900])



%%%%%%%%%%%%%%%%
% figure
% %region of gr
% ay=[0,4*pi/sqrt(3)];%lattice constant
% ax=[2*pi,0];
% 
% 
% coord0=ax-2*ay;
% Dx=4;
% Dy=4;
% LL=100;
% 
% gr_points=[];
% gr_values=[];
% for cx=0:LL
%     for cy=0:LL
%         coord=coord0+cx*Dx/LL*ax+cy*Dy/LL*ay;
%         gr_points=[gr_points;[coord]];
%         gr=get_gr(Sq_average,K_set,coord)*N/Np/Np;
%         gr_values=[gr_values,gr];
%     end
% end
% 
% 
% value_range=[min(real(gr_values)),max(real(gr_values))];
% for cq=1:length(gr_values)
% 
%     cq
%     gr_values(cq)
%     coord=gr_points(cq,:);
% 
%     Color=determine_color(value_range,real(gr_values(cq)));
%     plot(coord(1),coord(2),'marker','.','color',Color,'MarkerSize',10);hold on;%colorbar;
% 
% end
% 
% caxis(value_range);
% ax=gca;
% colormap(gca,mycmap);
% cb=colorbar;

%save("gr_28_14"+"_waa"+num2str(waa)+".mat","gr_values","gr_points","ax","ay","coord0","Dx","Dy");






set(gcf,'Position',[30 100 600 550])



function y=get_gr(Sq,K_set,R)
y=sum(sum(sum(Sq.*exp(i*(K_set(:,1,:,:)*R(1)+K_set(:,2,:,:)*R(2))))));
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



function [corners,in_FBZ,qx_set,qy_set]=determine_in_FBZ(k_set,q_list)



    g1=[1,0];
    g3=[-1/2,sqrt(3)/2];
    
    theta=pi/3;
    R_60=[cos(theta),-sin(theta);sin(theta),cos(theta)];
    corner1=[0.5,0.5/sqrt(3)];
    corners=zeros(7,2);
    for cc=1:7
        corners(cc,:)=R_60^(cc-1)*corner1';
        %plot(corners(cc,1),corners(cc,2),'ro');hold on;
    end
    
    
    g_set=zeros(6,2);
    
    for cc=1:size(g_set,1)
        %g_set(cc,:)=(4*pi/sqrt(3)/aM)*[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)];
        g_set(cc,:)=R_60^(cc-1)*g1';
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    klabel_set=k_set;
    
    nK=5;
    K_shift=-nK:nK;
    ca_cb=zeros((2*nK+1)*(2*nK+1),2);
    step=1;
    for ca=1:2*nK+1
        for cb=1:2*nK+1
            ca_cb(step,:)=[ca,cb];
            step=step+1;
        end
    end
    
    qx_set=zeros(size(q_list,1),1);
    qy_set=zeros(size(q_list,1),1);
    
    in_FBZ=zeros(size(q_list,1),1);
    
    for cc=1:size(q_list,1)
    
                k0=(q_list(cc,1))*g1+(q_list(cc,2))*g3;
                distance=norm(k0);
                distance0=norm(g1)/2;
    
                ds0=norm(k0);
                dss=zeros(6,1);
                for dd=1:length(dss)
                    dss(dd)=sqrt((k0(1)-g_set(dd,1))^2+(k0(2)-g_set(dd,2))^2);
                end
                qx_set(cc)=k0(1);
                qy_set(cc)=k0(2);
                if ds0<=min(dss)+0.01
                    %plot(k0(1),k0(2),'ko','MarkerSize',10,'LineWidth',4);hold on;
    
                    in_FBZ(cc)=1;
                    
                end
    
                % if (distance/distance0)<=1+1e-8
                %     plot(k0(1),k0(2),'ko','MarkerSize',10,'LineWidth',4);hold on;
                %     break;
                % end
            
        
    end
end


function in_FBZ=determine_in_FBZ_single(k_set,kk)



    g1=[1,0];
    g3=[-1/2,sqrt(3)/2];
    
    theta=pi/3;
    R_60=[cos(theta),-sin(theta);sin(theta),cos(theta)];
    corner1=[0.5,0.5/sqrt(3)];
    corners=zeros(7,2);
    for cc=1:7
        corners(cc,:)=R_60^(cc-1)*corner1';
        %plot(corners(cc,1),corners(cc,2),'ro');hold on;
    end
    
    
    g_set=zeros(6,2);
    
    for cc=1:size(g_set,1)
        %g_set(cc,:)=(4*pi/sqrt(3)/aM)*[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)];
        g_set(cc,:)=R_60^(cc-1)*g1';
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    klabel_set=k_set;
    
    nK=5;
    K_shift=-nK:nK;
    ca_cb=zeros((2*nK+1)*(2*nK+1),2);
    step=1;
    for ca=1:2*nK+1
        for cb=1:2*nK+1
            ca_cb(step,:)=[ca,cb];
            step=step+1;
        end
    end
    
    qx_set=0;
    qy_set=0;
    
    in_FBZ=0;
    
    
    
                k0=(kk(1))*g1+(kk(2))*g3;
                distance=norm(k0);
                distance0=norm(g1)/2;
    
                ds0=norm(k0);
                dss=zeros(6,1);
                for dd=1:length(dss)
                    dss(dd)=sqrt((k0(1)-g_set(dd,1))^2+(k0(2)-g_set(dd,2))^2);
                end
                qx_set=k0(1);
                qy_set=k0(2);
                if ds0<=min(dss)+0.01
                    %plot(k0(1),k0(2),'ko','MarkerSize',10,'LineWidth',4);hold on;
    
                    in_FBZ=1;
                    
                end
    
                % if (distance/distance0)<=1+1e-8
                %     plot(k0(1),k0(2),'ko','MarkerSize',10,'LineWidth',4);hold on;
                %     break;
                % end
            
        
    
end


function Sq=Sq_single_K(K_id_set,n_keep_set,Sq_all,id)
Sq=Sq_all(:,1)*0;
step=1;
    for cc=1:length(K_id_set)
        for dd=1:n_keep_set(cc)
            if cc==id
                Sq=Sq+Sq_all(:,step);
            end
            step=step+1;
        end
    end
Sq=Sq/n_keep_set(id);
end

function [kx_set,ky_set]=shift_kset(k_set)



    g1=[1,0];
    g3=[-1/2,sqrt(3)/2];
    
    theta=pi/3;
    R_60=[cos(theta),-sin(theta);sin(theta),cos(theta)];
    corner1=[0.5,0.5/sqrt(3)];
    corners=zeros(7,2);
    for cc=1:7
        corners(cc,:)=R_60^(cc-1)*corner1';
        %plot(corners(cc,1),corners(cc,2),'ro');hold on;
    end
    
    
    g_set=zeros(6,2);
    
    for cc=1:size(g_set,1)
        %g_set(cc,:)=(4*pi/sqrt(3)/aM)*[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)];
        g_set(cc,:)=R_60^(cc-1)*g1';
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    klabel_set=k_set;
    
    nK=5;
    K_shift=-nK:nK;
    ca_cb=zeros((2*nK+1)*(2*nK+1),2);
    step=1;
    for ca=1:2*nK+1
        for cb=1:2*nK+1
            ca_cb(step,:)=[ca,cb];
            step=step+1;
        end
    end
    
    kx_set=zeros(size(k_set,1),1);
    ky_set=zeros(size(k_set,1),1);
    
    for cc=1:size(k_set,1)
        kk=k_set(cc,:);
        for c1=-3:3
            for c2=-3:3
                k0=(kk(1)+c1)*g1+(kk(2)+c2)*g3;
                distance=norm(k0);
                distance0=norm(g1)/2;
    
                ds0=norm(k0);
                dss=zeros(6,1);
                for dd=1:length(dss)
                    dss(dd)=sqrt((k0(1)-g_set(dd,1))^2+(k0(2)-g_set(dd,2))^2);
                end

                if ds0<=min(dss)
                    %plot(k0(1),k0(2),'ko','MarkerSize',10,'LineWidth',4);hold on;
                    kx_set(cc)=k0(1);
                    ky_set(cc)=k0(2);

                    
                end

            end
        end
            
        
    end
end