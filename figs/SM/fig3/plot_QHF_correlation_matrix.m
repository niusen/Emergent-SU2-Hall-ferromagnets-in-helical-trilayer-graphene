clc;clear;close all;
restoredefaultpath;
pos_x0=0.08;
pos_y0=0.12;
length_x=0.15;
length_y=0.26;

%%addpath('D:\My Documents\Code\matlab\m\twist_bilayer\ED_TBG\Helical_trilayer\Chern_mosaic\split_band\ED_spectrum\scan_waa\theta_1.44_sublattice_4_6_8\N28_HF_2valley\order')
load order_theta1.44_waa0.6_epsilon4_Np14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0,0.5+pos_y0,length_x,length_y];
subplot('Position',pos);

sec=1;
N=size(k_set,1);
pos=sum(n_keep_all(1:sec-1))+1:sum(n_keep_all(1:sec));
nks=zeros(N,n_keep_all(sec));
for cc=1:length(pos)
    nks(:,cc)=diag(ob_all(:,:,pos(cc),pos(cc)))';
end

% plot(sort(nks(:,1)),'ro','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,2)),'bs','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,3)),'k*','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,1)),'ro','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,2)),'bs','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,3)),'k*','MarkerSize',10,'LineWidth',0.5);hold on;


axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)
ylabel('Eigenvalues','interpreter','latex','fontsize',13)
title('$\mathbf{\Phi}_{\Gamma}^{i=1,2,3}$','interpreter','latex','FontSize',16)
text(-6,1.2,"(a)",'interpreter','latex','FontSize',14)



%%%%
pos=[pos_x0+1/4,0.5+pos_y0,length_x,length_y];
subplot('Position',pos);
set(gca,'YTickLabel',[]);

sec=2;
N=size(k_set,1);
pos=sum(n_keep_all(1:sec-1))+1:sum(n_keep_all(1:sec));
nks=zeros(N,n_keep_all(sec));
for cc=1:length(pos)
    nks(:,cc)=diag(ob_all(:,:,pos(cc),pos(cc)))';
end

% plot(sort(nks(:,1)),'ro','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,2)),'bs','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,3)),'k*','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,4)),'mv','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,1)),'ro','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,2)),'bs','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,3)),'k*','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,4)),'mv','MarkerSize',10,'LineWidth',0.5);hold on;

axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)
title('$\mathbf{\Phi}_{\mathbf{M}_1}^{i=1,2,3,4}$','interpreter','latex','FontSize',16)
text(-6,1.2,"(b)",'interpreter','latex','FontSize',14)

%%%%%
pos=[pos_x0+2/4,0.5+pos_y0,length_x,length_y];
subplot('Position',pos);
set(gca,'YTickLabel',[]);

sec=3;
N=size(k_set,1);
pos=sum(n_keep_all(1:sec-1))+1:sum(n_keep_all(1:sec));
nks=zeros(N,n_keep_all(sec));
for cc=1:length(pos)
    nks(:,cc)=diag(ob_all(:,:,pos(cc),pos(cc)))';
end

% plot(sort(nks(:,1)),'ro','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,2)),'bs','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,3)),'k*','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,4)),'mv','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,1)),'ro','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,2)),'bs','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,3)),'k*','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,4)),'mv','MarkerSize',10,'LineWidth',0.5);hold on;

axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)
title('$\mathbf{\Phi}_{\mathbf{M}_2}^{i=1,2,3,4}$','interpreter','latex','FontSize',16)
text(-6,1.2,"(c)",'interpreter','latex','FontSize',14)

%%%%%

pos=[pos_x0+3/4,0.5+pos_y0,length_x,length_y];
subplot('Position',pos);
set(gca,'YTickLabel',[]);

sec=4;
N=size(k_set,1);
pos=sum(n_keep_all(1:sec-1))+1:sum(n_keep_all(1:sec));
nks=zeros(N,n_keep_all(sec));
for cc=1:length(pos)
    nks(:,cc)=diag(ob_all(:,:,pos(cc),pos(cc)))';
end

% plot(sort(nks(:,1)),'ro','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,2)),'bs','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,3)),'k*','MarkerSize',10,'LineWidth',0.5);hold on;
% plot(sort(nks(:,4)),'mv','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,1)),'ro','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,2)),'bs','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,3)),'k*','MarkerSize',10,'LineWidth',0.5);hold on;
plot((nks(:,4)),'mv','MarkerSize',10,'LineWidth',0.5);hold on;

axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)
title('$\mathbf{\Phi}_{\mathbf{M}_3}^{i=1,2,3,4}$','interpreter','latex','FontSize',16)
text(-6,1.2,"(d)",'interpreter','latex','FontSize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(0,0), (0,0.5), (0.5,0), (0.5,0.5)
%Gamma, M2, M1, M3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=[pos_x0,pos_y0,length_x,length_y];
subplot('Position',pos);


load 2_4.mat;
N=size(k_set,1);
eus=eig(rho);
plot(sort(eus),'r.','MarkerSize',10);hold on;

axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)
ylabel('Eigenvalues','interpreter','latex','fontsize',13)
%title('$\tilde{\psi}_{\mathbf{M}_1,\mathbf{M}_3}$','interpreter','latex','FontSize',16)
% title('$\tilde{\psi}_{\mathbf{M}_2,\mathbf{M}_2+\mathbf{M}_1}$','interpreter','latex','FontSize',16)
%title('$\tilde{\psi}_{\mathbf{M}_2,\mathbf{M}_3}, \mathbf{Q}=\mathbf{M}_1$','interpreter','latex','FontSize',12)
title('$\mathbf{Q}=\mathbf{M}_1$','interpreter','latex','FontSize',12)
text(-6,1.2,"(e)",'interpreter','latex','FontSize',14)



%%%%
pos=[pos_x0+1/4,pos_y0,length_x,length_y];
subplot('Position',pos);


load 3_4.mat;
N=size(k_set,1);
eus=eig(rho);
plot(sort(eus),'r.','MarkerSize',10);hold on;

axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)
ylabel('Eigenvalues','interpreter','latex','fontsize',13)
%title('$\tilde{\psi}_{\mathbf{M}_2,\mathbf{M}_3}$','interpreter','latex','FontSize',16)
title('$\mathbf{Q}=\mathbf{M}_2$','interpreter','latex','FontSize',12)
text(-6,1.2,"(f)",'interpreter','latex','FontSize',14)
%%%%
pos=[pos_x0+2/4,pos_y0,length_x,length_y];
subplot('Position',pos);

load 2_3.mat;
N=size(k_set,1);
eus=eig(rho);
plot(sort(eus),'r.','MarkerSize',10);hold on;

axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)
ylabel('Eigenvalues','interpreter','latex','fontsize',13)
%title('$\tilde{\psi}_{\mathbf{M}_1,\mathbf{M}_2}$','interpreter','latex','FontSize',16)
title('$\mathbf{Q}=\mathbf{M}_3$','interpreter','latex','FontSize',12)
text(-6,1.2,"(g)",'interpreter','latex','FontSize',14)







%%%%
pos=[pos_x0+3/4,pos_y0,length_x,length_y];
subplot('Position',pos);

load general_rho.mat;
N=size(k_set,1);
eus=eig(rho);
plot(sort(eus),'r.','MarkerSize',10);hold on;

axis([1,N,0,1]);

xticks([0 7 14  21 28 ])
set(gca,'fontsize',12)
xlabel('Index','interpreter','latex','fontsize',13)
ylabel('Eigenvalues','interpreter','latex','fontsize',13)
%title('$\tilde{\psi}_{\mathbf{\Gamma},\mathbf{M}_3}$','interpreter','latex','FontSize',16)
title("General "+"$\mathbf{Q}$",'interpreter','latex','FontSize',12)
text(-6,1.2,"(h)",'interpreter','latex','FontSize',14)














set(gcf,'Position',[100 100 800 380])
