function [H,eu,ev,q1,q2,q3,basis,k_set,H_kinetic_derivative_x,H_kinetic_derivative_y]=band_single_k(parameters,Phi,k,nK)
valley=parameters.valley;
AB_potential=parameters.AB_potential;
v=parameters.v;
% alpha=0;
wab=parameters.wab;
waa=parameters.waa;
ktheta=parameters.ktheta;


theta=parameters.theta_/180*pi;






q1=ktheta*[0,-1];
q2=ktheta*[sqrt(3)/2,1/2];
q3=ktheta*[-sqrt(3)/2,1/2];

b1m=sqrt(3)*ktheta*[0.5, -sqrt(3)/2];
b2m=sqrt(3)*ktheta*[0.5, sqrt(3)/2];

% K1=ktheta*[-sqrt(3)/2,-1/2];
% K2=ktheta*[-sqrt(3)/2,1/2];
% K3=ktheta*[-sqrt(3)/2,1/2+1];

K1=ktheta*[0,-1];
K2=ktheta*[0,0];
K3=ktheta*[0,1];


phi=2*pi/3;

n=1;
T1=[waa wab*exp(-i*phi*(n-1));wab*exp(i*phi*(n-1)) waa];
n=2;
T3=[waa wab*exp(-i*phi*(n-1));wab*exp(i*phi*(n-1)) waa];
n=3;
T2=[waa wab*exp(-i*phi*(n-1));wab*exp(i*phi*(n-1)) waa];


% T1=[waa wab;wab waa];
% 
% T2=[waa*exp(-i*phi) wab;wab*exp(i*phi) waa*exp(-i*phi)];
% 
% T3=[waa*exp(i*phi) wab;wab*exp(-i*phi) waa*exp(i*phi)];

%[exp(i*d_vector*q1'),exp(i*d_vector*q2'),exp(i*d_vector*q3')]

T1_12=T1*exp(-i*Phi(1));
T2_12=T2*exp(-i*Phi(2));
T3_12=T3*exp(-i*Phi(3));

T1_23=T1*exp(i*Phi(1));
T2_23=T2*exp(i*Phi(2));
T3_23=T3*exp(i*Phi(3));



if valley==1
    %time reversal
    T1_12=conj(T1_12);
    T2_12=conj(T2_12);
    T3_12=conj(T3_12);

    T1_23=conj(T1_23);
    T2_23=conj(T2_23);
    T3_23=conj(T3_23);
end

if valley==-1
    vx_coe=1;
elseif valley==1
    vx_coe=-1;
end




k_set=zeros(2*nK+1,2*nK+1,2);
for c1=1:2*nK+1
    for c2=1:2*nK+1
        k_set(c1,c2,:)=k+(c1-nK-1)*b1m+(c2-nK-1)*b2m;
        %plot(k_set(c1,c2,1),k_set(c1,c2,2),'ks');hold on;
    end
end

layer_set=[1,2,3];
sz_set=[1,-1];%sublattice

%construct single particle basis
N_basis=prod([size(k_set,1),size(k_set,2),length(layer_set),length(sz_set)]);
basis=linspace(1,N_basis,N_basis);
basis=reshape(basis,[size(k_set,1),size(k_set,2),length(layer_set),length(sz_set)]);

H_kinetic=zeros(N_basis,N_basis);
H_U=zeros(N_basis,N_basis);
H_AB_potential=zeros(N_basis,N_basis);
H_kinetic_derivative_x=zeros(N_basis,N_basis);
H_kinetic_derivative_y=zeros(N_basis,N_basis);


%kinetic energy
for c1a=1:size(k_set,1)
    for c1b=1:size(k_set,2)
        k=k_set(c1a,c1b,:);
        k=k(:)';
        for c2=1:length(layer_set)
            for c3=1:length(sz_set)
                if c3==1
                    c31=1;
                    c32=2;
                elseif c3==2
                    c31=2;
                    c32=1;
                end
                ind1=basis(c1a,c1b,c2,c31);
                ind2=basis(c1a,c1b,c2,c32);

            
                if c31==2%sublattice B
                    if c2==1
                        k_=([k(1);k(2)]-valley*[K1(1);K1(2)]);
                        H_kinetic(ind1,ind2)=-v*(vx_coe*k_(1)-i*k_(2));

                        H_kinetic_derivative_x(ind1,ind2)=-vx_coe*v*(1);
                        H_kinetic_derivative_y(ind1,ind2)=-v*(-i);
                    elseif c2==2
                        k_=([k(1);k(2)]-valley*[K2(1);K2(2)]);
                        H_kinetic(ind1,ind2)=-v*(vx_coe*k_(1)-i*k_(2));

                        H_kinetic_derivative_x(ind1,ind2)=-vx_coe*v*(1);
                        H_kinetic_derivative_y(ind1,ind2)=-v*(-i);
                    elseif c2==3
                        k_=([k(1);k(2)]-valley*[K3(1);K3(2)]);
                        H_kinetic(ind1,ind2)=-v*(vx_coe*k_(1)-i*k_(2));

                        H_kinetic_derivative_x(ind1,ind2)=-vx_coe*v*(1);
                        H_kinetic_derivative_y(ind1,ind2)=-v*(-i);
                    end

                elseif c31==1%sublattice A
                    if c2==1
                        k_=([k(1);k(2)]-valley*[K1(1);K1(2)]);
                        H_kinetic(ind1,ind2)=-v*(vx_coe*k_(1)+i*k_(2));

                        H_kinetic_derivative_x(ind1,ind2)=-vx_coe*v*(1);
                        H_kinetic_derivative_y(ind1,ind2)=-v*(i);
                    elseif c2==2
                        k_=([k(1);k(2)]-valley*[K2(1);K2(2)]);
                        H_kinetic(ind1,ind2)=-v*(vx_coe*k_(1)+i*k_(2));

                        H_kinetic_derivative_x(ind1,ind2)=-vx_coe*v*(1);
                        H_kinetic_derivative_y(ind1,ind2)=-v*(i);
                    elseif c2==3
                        k_=([k(1);k(2)]-valley*[K3(1);K3(2)]);
                        H_kinetic(ind1,ind2)=-v*(vx_coe*k_(1)+i*k_(2));

                        H_kinetic_derivative_x(ind1,ind2)=-vx_coe*v*(1);
                        H_kinetic_derivative_y(ind1,ind2)=-v*(i);
                    end

                end

            end
        end
    end
end


%U term
for c1a=1:size(k_set,1)
    for c1b=1:size(k_set,2)

        for d1a=1:size(k_set,1)
            for d1b=1:size(k_set,2)
                
                for c2=1:length(layer_set)-1
                    if c2==1
                        c21=1;
                        c22=2;
                    elseif c2==2
                        c21=2;
                        c22=3;
                    end
                    

                    for c31=1:length(sz_set)%sublattice
                        for c32=1:length(sz_set)%sublattice

                            ind1=basis(c1a,c1b,c21,c31);
                            ind2=basis(d1a,d1b,c22,c32);
                            
                            if c21==1
                                kp=k_set(c1a,c1b,:);
                                kp=kp(:)'-K1*valley;
                            elseif c21==2
                                kp=k_set(c1a,c1b,:);
                                kp=kp(:)'-K2*valley;
                            elseif c21==3
                                kp=k_set(c1a,c1b,:);
                                kp=kp(:)'-K3*valley;
                            end

                            if c22==1
                                k=k_set(d1a,d1b,:);
                                k=k(:)'-K1*valley;
                            elseif c22==2
                                k=k_set(d1a,d1b,:);
                                k=k(:)'-K2*valley;
                            elseif c22==3
                                k=k_set(d1a,d1b,:);
                                k=k(:)'-K3*valley;
                            end
                            
                            if (c21==1)&&(c22==2)
                                if (sum(abs(kp-(k-valley*q1)))<(1e-10)*norm(q1)) 
                                    H_U(ind1,ind2)=T1_12(c31,c32);
                                elseif (sum(abs(kp-(k-valley*q2)))<(1e-10)*norm(q1)) 
                                    H_U(ind1,ind2)=T2_12(c31,c32);
                                elseif (sum(abs(kp-(k-valley*q3)))<(1e-10)*norm(q1)) 
                                    H_U(ind1,ind2)=T3_12(c31,c32);
                                end
                            elseif (c21==2)&&(c22==3)
                                if (sum(abs(kp-(k-valley*q1)))<(1e-10)*norm(q1)) 
                                    H_U(ind1,ind2)=T1_23(c31,c32);
                                elseif (sum(abs(kp-(k-valley*q2)))<(1e-10)*norm(q1)) 
                                    H_U(ind1,ind2)=T2_23(c31,c32);
                                elseif (sum(abs(kp-(k-valley*q3)))<(1e-10)*norm(q1)) 
                                    H_U(ind1,ind2)=T3_23(c31,c32);
                                end
                            end
                        end
                    end
                            
                end
            end
        end
    end
end

%sublattice potential
for c1a=1:size(k_set,1)
    for c1b=1:size(k_set,2)
        for c2=1:length(layer_set)
            for c3=1:length(sz_set)%sublattice
                ind=basis(c1a,c1b,c2,c3);
                H_AB_potential(ind,ind)=AB_potential(c2,c3);
            end  
        end
    end
end
H=H_kinetic+H_U+H_U'+H_AB_potential;

[ev,eu]=eig(H);
eu=diag(eu);

end

