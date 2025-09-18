function [Pset,N,g1,g3,kk1,kk2,V1_,V2_]=sort_K(V1,V2,K_set)

cc=1;
g1=[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)];
cc=3;
g3=[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)];


g1_3d=[g1(:);0];
g3_3d=[g3(:);0];
z_3d=[0;0;1];

a1=2*pi*cross(z_3d,g3_3d)/(g1_3d'*cross(z_3d,g3_3d));
a2=2*pi*cross(z_3d,g1_3d)/(g3_3d'*cross(z_3d,g1_3d));

V1_3d=V1(1)*a1+V1(2)*a2;
V2_3d=V2(1)*a1+V2(2)*a2;

N_=abs(dot(cross(V1_3d,V2_3d),z_3d)/(dot(cross(a1,a2),z_3d)));
assert((abs(round(N_))-abs(N_))/abs(round(N_))<1e-14);
N=abs(round(N_));





if (V1(2)==0)&&(V2(1)==0)
    k_set=[];
    step=0;
    for c2=1:V2(2)
        for c1=1:V1(1)
    
            label=[(c1-1)/V1(1),(c2-1)/V2(2)];
                step=step+1;
                k_set=[k_set;label];
        end
    end
else
    LL=40;
    k_set=[];
    step=0;
    for c2=1:LL
        for c1=1:LL
            p1=c1-round(LL/2);
            p2=c2-round(LL/2);
            V1_=[V1(2),-V1(1)];
            V2_=[V2(2),-V2(1)];
            label=(p1*V2_+p2*V1_)/N;
            if (0<=label(1))&(label(1)<1)&(0<=label(2))&(label(2)<1)
                step=step+1;
                k_set=[k_set;label];
            end
            step=step+1;
        end
    end
end
k_set;


    V1_=[V1(2),-V1(1)]/N;
    V2_=[V2(2),-V2(1)]/N;
    %label=(p1*V2_+p2*V1_)/N;
    
    Vm=[V2(2),V1(2);-V2(1),-V1(1)];
    Pset=zeros(size(K_set,1),2);
    for cc=1:size(K_set,1)
        kk=K_set(cc,:)';
        Pset(cc,:)=inv(Vm)*kk*N;
        
    end
    
    
    
    kk1=V2_(1)*g1+V2_(2)*g3;
    kk2=V1_(1)*g1+V1_(2)*g3;
    
    
    
    
end