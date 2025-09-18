clc;clear;close all;
waa=0.6;
filenm="order_theta1.44_waa"+num2str(waa)+"_epsilon4_Np16";
load(filenm);


N=size(k_set,1);
%Q2, G1/2
sector1=3;
sector2=4;
Rng=134112;
[psi1,rho1]=opt_two_component(Rng,sector1,sector2,n_keep_all,k_set,order_all,ob_all,K_ind_correct_order);



%%%%%%%%%%%%%%%%%%
Q1=-[0,0.5];
Q2=[0.5,0];
Q3=[0.5,0.5];

partners1=find_partner(k_set,1:N,Q2);%k+G1/2
partners2=find_partner(k_set,1:N,Q1);%k+G2/2
partners3=find_partner(k_set,1:N,Q3);%k+G1/2+G2/2

pos=find((1:N)-partners1'>0);
assert(length(pos)==N/2);

half_FBZ0=pos';
half_FBZ1=partners1(pos);
half_FBZ2=partners2(pos);
half_FBZ3=partners3(pos);

cc=1;
k_set([half_FBZ0(cc),half_FBZ1(cc),half_FBZ2(cc),half_FBZ3(cc)],:)

%construct pseudo-spin basis
option=2
up_indices=[half_FBZ0,half_FBZ1];
[vecs_up,eus_up]=construct_basis_pair(rho1,up_indices,"occupy",option)
dn_indices=[half_FBZ2,half_FBZ3];
[vecs_dn,eus_dn]=construct_basis_pair(rho1,dn_indices,"unoccupy",option)

value_up=Sz_from_rho(rho1,up_indices,vecs_up)
value_dn=Sz_from_rho(rho1,dn_indices,vecs_dn)

Ng=size(ob_all,3);
Op_z=zeros(Ng,Ng);
for c1=1:Ng
    for c2=1:Ng
        rho=ob_all(:,:,c1,c2);
        value_up=Sz_from_rho(rho,up_indices,vecs_up);
        value_dn=Sz_from_rho(rho,dn_indices,vecs_dn);
        Op_z(c1,c2)=value_up-value_dn;
    end
end
    

sz_values=eig(Op_z)
plot(sort(real(sz_values))/2,'ro');

save("Sz_waa"+num2str(waa)+"_option"+num2str(option)+".mat","Op_z");

function value=Sz_from_rho(rho,indices,vecs)
    value=0;
    siz=size(indices,1);
    for cc=1:siz
        ind_=indices(cc,:);
        st=vecs(cc,:);
        assert(length(st)==2);
        st=[st(1);st(2)];
        rho_sub=rho(ind_,ind_);
        value=value+st'*rho_sub*st;
    end
end

function [vecs,eus]=construct_basis_pair(rho,indices,fill,option)
    ks=indices(:,1);
    k_G1s=indices(:,2);
    assert(length(ks)==length(k_G1s));
    siz=length(ks);
    vecs=zeros(siz,2);
    eus=[];
    for cc=1:siz
        rho_sub=rho([ks(cc),k_G1s(cc)],[ks(cc),k_G1s(cc)]);
        assert(norm(rho_sub-rho_sub')<1e-10);
        rho_sub=(rho_sub+rho_sub')/2;
        [ev,eu]=eig(rho_sub);
        
        if fill=="occupy"
            vecs(cc,:)=ev(:,2);
            eus=[eus,eu(2,2)];
        elseif fill=="unoccupy"
            if option==1
                vecs(cc,:)=ev(:,1);
            elseif option==2
                %choose basis according to the observation that, the off
                %diagonal correlation matrix elements between states with opposite order parameters have opposite signs.
                %thus the partner of (alpha,beta) is (alpha,-beta), which are
                %not exactly othogonal .
                %this works only when cluster is highly symmetric (I guess need inversion symmetry)
                vecs(cc,:)=[ev(1,2);-ev(2,2)];
            end
            eus=[eus,eu(1,1)];
        end
    end
    %eus=sort(eus(:));
end


function [psi,rho]=opt_two_component(Rng,sector1,sector2,n_keep_all,k_set,order_all,ob_all,K_ind_correct_order)

sector_pos1=sum(n_keep_all(1:sector1-1))+1:sum(n_keep_all(1:sector1));
sector_pos2=sum(n_keep_all(1:sector2-1))+1:sum(n_keep_all(1:sector2));
K_ind_correct_order;

k_select=3;
pos=[order_all{sector1}(k_select),order_all{sector2}(k_select)];

Ng=sum(n_keep_all);
N=size(k_set,1);

ob_select=ob_all(pos,pos,:,:);
%ob_select=reshape(ob_select(1,3,:,:),[Ng,Ng]);



f = @(x)get_max_eu(x,ob_all,Ng,N,sector_pos1,sector_pos2);
rng(Rng)
x0=randn(2*Ng,1);
[ xsol fval]=fminsearch(f,x0,optimset('Display','off','MaxFunEvals',10000));
psi=xsol(1:Ng)+i*xsol(Ng+1:2*Ng);
psi=psi_project(psi,sector_pos1,sector_pos2);
eus=get_eu_all(psi,ob_all,Ng,N,sector_pos1,sector_pos2)

[eus,rho]=get_eu_all(psi,ob_all,Ng,N,sector_pos1,sector_pos2);



filenm=num2str(sector1)+"_"+num2str(sector2)+".mat";
%save(filenm,"rho","psi")

end

function psi=psi_project(psi,sector_pos1,sector_pos2)
psi1=psi;
psi1(sector_pos1)=psi1(sector_pos1)*0;
psi1(sector_pos2)=psi1(sector_pos2)*0;
psi=psi-psi1;
psi=psi/norm(psi);
end
function eu_res=get_max_eu(psi_all,ob_select,Ng,N,sector_pos1,sector_pos2)

    psi=psi_all(1:Ng)+psi_all(Ng+1:2*Ng)*i;
    psi=psi_project(psi,sector_pos1,sector_pos2);
    eus=get_eu(psi,ob_select,Ng,N,sector_pos1,sector_pos2);
   
    eu_res=sum(1-eus(N/2+1:N));
end


function eus=get_eu(psi,ob_all,Ng,N,sector_pos1,sector_pos2)
    psi=psi_project(psi,sector_pos1,sector_pos2);
    rho=zeros(N,N)*i;
    for c1=1:N
        for c2=1:N
            ob_=reshape(ob_all(c1,c2,:,:),[Ng,Ng]);
            rho(c1,c2)=contract(ob_,psi);
        end
    end
    assert(norm(rho-rho')<1e-10);
    rho=(rho+rho')/2;
    eus=eig(rho);

end

function [eus,rho]=get_eu_all(psi,ob_all,Ng,N,sector_pos1,sector_pos2)
    psi=psi_project(psi,sector_pos1,sector_pos2);
    rho=zeros(N,N)*i;
    for c1=1:N
        for c2=1:N
            ob_=reshape(ob_all(c1,c2,:,:),[Ng,Ng]);
            rho(c1,c2)=contract(ob_,psi);
        end
    end
    assert(norm(rho-rho')<1e-10);
    rho=(rho+rho')/2;
    eus=eig(rho);

end

function y=contract(rho,psi)
y=psi'*rho*psi;

end



function partners=find_partner(k_set0,index_list,Q)
    k_set=k_set0(index_list,:);
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