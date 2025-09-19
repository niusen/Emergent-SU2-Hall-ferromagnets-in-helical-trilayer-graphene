clc;clear;close all;
Q1=-[0,0.5];
Q2=[0.5,0];
Q3=[0.5,0.5];

Q=Q1;

subplot(3,4,1);
s_direction="x";
plot_texture(s_direction,Q1,Q2,Q3,Q);

text(-1.8,0,'$\tilde{\psi}_{\mathbf{M}_2,\mathbf{M}_3}, \mathbf{Q}=\mathbf{M}_1$','interpreter','latex','fontsize',14)


subplot(3,4,2);
s_direction="y";
plot_texture(s_direction,Q1,Q2,Q3,Q);

subplot(3,4,3);
s_direction="z";
plot_texture(s_direction,Q1,Q2,Q3,Q);

subplot(3,4,4);
s_direction="xy";
plot_texture(s_direction,Q1,Q2,Q3,Q);



%%%%%%%%%%%%%%%%%%%%
Q=Q2;

subplot(3,4,5);
s_direction="x";
plot_texture(s_direction,Q1,Q2,Q3,Q);

text(-1.8,0,'$\tilde{\psi}_{\mathbf{M}_1,\mathbf{M}_3}, \mathbf{Q}=\mathbf{M}_2$','interpreter','latex','fontsize',14)

subplot(3,4,6);
s_direction="y";
plot_texture(s_direction,Q1,Q2,Q3,Q);

subplot(3,4,7);
s_direction="z";
plot_texture(s_direction,Q1,Q2,Q3,Q);

subplot(3,4,8);
s_direction="xy";
plot_texture(s_direction,Q1,Q2,Q3,Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%
Q=Q3;

subplot(3,4,9);
s_direction="x";
plot_texture(s_direction,Q1,Q2,Q3,Q);

text(-1.8,0,'$\tilde{\psi}_{\mathbf{M}_2,\mathbf{M}_1}, \mathbf{Q}=\mathbf{M}_3$','interpreter','latex','fontsize',14)

subplot(3,4,10);
s_direction="y";
plot_texture(s_direction,Q1,Q2,Q3,Q);

subplot(3,4,11);
s_direction="z";
plot_texture(s_direction,Q1,Q2,Q3,Q);

subplot(3,4,12);
s_direction="xy";
plot_texture(s_direction,Q1,Q2,Q3,Q);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



set(gcf,'Position',[30 100 1200 600])







function plot_texture(s_direction,Q1,Q2,Q3,Q)
load("mymap.mat")

load order_theta1.44_waa0.6_epsilon4_Np14;





%load correlation matrix
if prod(Q==Q1)
    load 2_4.mat;
elseif prod(Q==Q2)
    load 3_4.mat;
elseif prod(Q==Q3)
    load 2_3.mat;
end



partners=find_partner(k_set,Q);

[k_set_FBZ,kx_set,ky_set]=shift_kset(k_set);

%shift special k point by G
if prod(Q==Q3)
    k_set_FBZ(12,1)=k_set_FBZ(12,1)-1;
    kx_set(12)=kx_set(12)-1;
end


[corners,in_FBZ]=determine_in_FBZ(k_set);

nk=zeros(size(k_set,1),1);

if prod(Q==Q1)
    for cc=1:3
        p1=corners(cc,:);
        p2=corners(cc+1,:);
        plot([p1(1),p2(1)],[p1(2),p2(2)],'k-','LineWidth',6);hold on;
    end
    p1=corners(1,:);
    p2=corners(3+1,:);
    plot(linspace(p1(1),p2(1),20),linspace(p1(2),p2(2),20),'k--','LineWidth',3);hold on;
elseif prod(Q==Q2)
    for cc=2:4
        p1=corners(cc,:);
        p2=corners(cc+1,:);
        plot([p1(1),p2(1)],[p1(2),p2(2)],'k-','LineWidth',6);hold on;
    end
    p1=corners(2,:);
    p2=corners(4+1,:);
    plot(linspace(p1(1),p2(1),20),linspace(p1(2),p2(2),20),'k--','LineWidth',3);hold on;
elseif prod(Q==Q3)
    for cc=3:5
        p1=corners(cc,:);
        p2=corners(cc+1,:);
        plot([p1(1),p2(1)],[p1(2),p2(2)],'k-','LineWidth',6);hold on;
    end
    p1=corners(3,:);
    p2=corners(5+1,:);
    plot(linspace(p1(1),p2(1),20),linspace(p1(2),p2(2),20),'k--','LineWidth',3);hold on;
end


N=size(k_set,1);
value_range=[-0.5,0.5];

%half FBZ
if prod(Q==Q1)
    extra_inds=[17,20];
elseif prod(Q==Q2)
    extra_inds=[8,22];
    
elseif prod(Q==Q3)
    extra_inds=[7,12,18];
end


sx=[0,1;1,0]/2;
sy=[0,-i;i,0]/2;
sz=[1,0;0,-1]/2;

for cq=1:N
    posa=cq;
    posb=partners(cq);

    % if sum(cq==extra_inds)>0
    %     plot(kx_set(cq),ky_set(cq),'marker','.','color','r','MarkerSize',38);hold on;%colorbar;
    % end
    % if (abs(k_set_FBZ(posa,1)-k_set_FBZ(posb,1)+Q(1))<1e-10)&&(abs(k_set_FBZ(posa,2)-k_set_FBZ(posb,2)+Q(2))<1e-10)
    % 
    %     %Color=determine_color(value_range,real(nk(cq)));
    %     plot(kx_set(cq),ky_set(cq),'marker','.','color','r','MarkerSize',38);hold on;%colorbar;
    % % elseif (abs(k_set_FBZ(posb,1)-k_set_FBZ(posa,1)+Q(1))<1e-10)&&(abs(k_set_FBZ(posb,2)-k_set_FBZ(posa,2)+Q(2))<1e-10)
    % %     %Color=determine_color(value_range,real(nk(cq)));
    % %     plot(kx_set(cq),ky_set(cq),'marker','.','color','b','MarkerSize',38);hold on;%colorbar;
    % else
    %     cq
    %     k_set_FBZ(posa,:)
    %     k_set_FBZ(posb,:)
    % end

    rho_sub=rho([posa,posb],[posa,posb]);
    [ev,eu]=eig(rho_sub);
    assert(eu(2,2)>0.8);
    assert(eu(1,1)<0.2);
    psi=ev(:,2);

    if s_direction=="x"
        ob=psi'*sx*psi;
    elseif s_direction=="y"
        ob=psi'*sy*psi;
    elseif s_direction=="z"
        ob=psi'*sz*psi;
    elseif s_direction=="xy"
        ob1=psi'*sx*psi;
        ob2=psi'*sy*psi;
        ob=sqrt(ob1^2+ob2^2);
    end

    Color=determine_color(value_range,real(ob));
    if sum(cq==extra_inds)>0
        plot(kx_set(cq),ky_set(cq),'marker','.','color',Color,'MarkerSize',38);hold on;%colorbar;
    end
    if (abs(k_set_FBZ(posa,1)-k_set_FBZ(posb,1)+Q(1))<1e-10)&&(abs(k_set_FBZ(posa,2)-k_set_FBZ(posb,2)+Q(2))<1e-10)
     
        
        plot(kx_set(cq),ky_set(cq),'marker','.','color',Color,'MarkerSize',38);hold on;%colorbar;
    end
    
end
axis([-0.6,0.6,-0.6,0.6]);
caxis(value_range);
ax=gca;
colormap(gca,mycmap);
cb=colorbar;
%ylabel(cb,'n(k)','Rotation',270)


if s_direction=="x"
    title("$\langle \tau_{x} \rangle$",'interpreter','latex')
elseif s_direction=="y"
    title("$\langle \tau_{y} \rangle$",'interpreter','latex')
elseif s_direction=="z"
    title("$\langle \tau_{z} \rangle$",'interpreter','latex')
elseif s_direction=="xy"
title("$\sqrt{\langle \tau_{x} \rangle^2+\langle \tau_{y} \rangle^2}$",'interpreter','latex')
end


ax.FontSize=15;
axis off


end





function partners=find_partner(k_set,Q)
    siz=size(k_set,1);
    partners=zeros(siz,1);
    for cc=1:siz
        kk=k_set(cc,:);
        
        kk=kk+Q;
        dk1=kk(1)-k_set(:,1);
        dk2=kk(2)-k_set(:,2);

        cond1=(abs(round(dk1)-dk1)<1e-10);
        cond2=(abs(round(dk2)-dk2)<1e-10);

        pos=find((cond1.*cond2)==1);
        assert(length(pos)==1);
        partners(cc)=pos;
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



function [corners,in_FBZ,qx_set,qy_set]=determine_in_FBZ(k_set)



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
    
    qx_set=zeros(size(k_set,1),1);
    qy_set=zeros(size(k_set,1),1);
    
    in_FBZ=zeros(size(k_set,1),1);
    
    for cc=1:size(k_set,1)
    
                k0=(k_set(cc,1))*g1+(k_set(cc,2))*g3;
                distance=norm(k0);
                distance0=norm(g1)/2;
    
                ds0=norm(k0);
                dss=zeros(6,1);
                for dd=1:length(dss)
                    dss(dd)=sqrt((k0(1)-g_set(dd,1))^2+(k0(2)-g_set(dd,2))^2);
                end
                qx_set(cc)=k0(1);
                qy_set(cc)=k0(2);
                if ds0<=min(dss)+1e-10
                    %plot(k0(1),k0(2),'ko','MarkerSize',10,'LineWidth',4);hold on;
    
                    in_FBZ(cc)=1;
                    
                end
    
                % if (distance/distance0)<=1+1e-8
                %     plot(k0(1),k0(2),'ko','MarkerSize',10,'LineWidth',4);hold on;
                %     break;
                % end
            
        
    end
end



function [k_set_FBZ,kx_set,ky_set]=shift_kset(k_set)



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
    k_set_FBZ=k_set;
    
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

                if ds0<=min(dss)+1e-10
                    %plot(k0(1),k0(2),'ko','MarkerSize',10,'LineWidth',4);hold on;
                    kx_set(cc)=k0(1);
                    ky_set(cc)=k0(2);
                    k_set_FBZ(cc,1:2)=[(kk(1)+c1),+(kk(2)+c2)];

                    
                end

            end
        end
            
        
    end
end