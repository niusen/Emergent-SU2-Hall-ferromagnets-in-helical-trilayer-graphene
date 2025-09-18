using Distributed
@everywhere using Distributed 
using Revise, TensorKit
using JLD2,MAT
using KrylovKit
using JSON
using Random
using Dates
@everywhere using SparseArrays,Combinatorics
#using Arpack
using KrylovKit

cd(@__DIR__)
@everywhere include("ED_lib.jl")
include("model_lib.jl")
@show num_logical_cores = Sys.CPU_THREADS
####################
import LinearAlgebra.BLAS as BLAS
@show n_blas_thread=1;
BLAS.set_num_threads(n_blas_thread);
@show nthread=Threads.nthreads();#number of thread can only be set in job.sh file
pid=getpid();
println("pid="*string(pid));
####################
@everywhere function main(band_name,parameters,Np,V1,V2,flux_1,flux_2,simplify_interaction,nsec,n_keep,band_state_origin,Mode)
    println(parameters)
    global band_model
    
    #############################################
    theta_,theta,aM,ktheta,waa,wab,v,valley,AB_potential,Phi,kappa=get_parameters(parameters);
    #############################################
    println("theta="*string(theta_))
    println("Np="*string(Np))
    println("V1="*string(V1))
    println("V2="*string(V2))
    
    @assert 0<=flux_1<=1
    @assert 0<=flux_2<=1
 
    nK=6;#number of g vectors in each side

    H,eu,ev,g1,g3,basis,g_lattice=band_model(parameters,[0,0],nK);   
    # println(eu)
    #####################################

    #g_lattice: multiples of g vector

    g1_3d=[g1[1],g1[2], 0];
    g3_3d=[g3[1],g3[2], 0];
    z_3d=[0,0,1];

    a1=2*pi*cross(z_3d,g3_3d)/(dot(g1_3d, cross(z_3d,g3_3d)));
    a2=2*pi*cross(z_3d,g1_3d)/(dot(g3_3d, cross(z_3d,g1_3d)));


    V1_3d=V1[1]*a1+V1[2]*a2;
    V2_3d=V2[1]*a1+V2[2]*a2;

    N_=abs(cross(V1_3d,V2_3d)[3]/cross(a1,a2)[3]);
    @assert (Int(abs(round(N_)))-abs(N_))/abs(round(N_))<1e-14;
    N=Int(abs(round(N_)));
    println("cluster size N="*string(N))
    G1_3d=2*pi*cross(z_3d,V2_3d)/(dot(V1_3d, cross(z_3d,V2_3d)));
    G2_3d=2*pi*cross(z_3d,V1_3d)/(dot(V2_3d, cross(z_3d,V1_3d)));

    G1=G1_3d[1:2];
    G2=G2_3d[1:2];


    #################

    if (V1[2]==0)&&(V2[1]==0)
        k_set=zeros(N,2);
        Step=0;
        for c2=1:V2[2]
            for c1=1:V1[1]
                label=[(c1-1)/V1[1],(c2-1)/V2[2]];
                Step=Step+1;
                k_set[Step,:]=label;
            end
        end

        #determine shift of k points
        #old method
        # k1_shift=flux_1*[1/abs(V1[1]),0];
        # k2_shift=flux_2*[0,1/abs(V2[2])];

        #new method
        V1_=[V1[2],-V1[1]];
        V2_=[V2[2],-V2[1]];
        k1_shift=flux_1*V2_/N;
        k2_shift=-flux_2*V1_/N;

    else

        LL=40;
        k_set=zeros((LL*LL,2));
        Step=0;
        for c2 in range(1,LL)
            for c1 in range(1,LL)
                p1=Int(c1-round(LL/2));
                p2=Int(c2-round(LL/2));

                V1_=[V1[2],-V1[1]];
                V2_=[V2[2],-V2[1]];
                label=[(p1*V2_[1]+p2*V1_[1])//N, (p1*V2_[2]+p2*V1_[2])//N];
                
                if (0<=label[1])&(label[1]<1)&(0<=label[2])&(label[2]<1)
                    Step=Step+1;
                    k_set[Step,:]=label;
                end
            end
        end
        k_set=k_set[1:Step,:];

        #determine shift of k points
        k1_shift=flux_1*V2_/N;
        k2_shift=-flux_2*V1_/N;
    end
    println(k_set)
    @assert size(k_set,1)==N;
    klabel_set=k_set;
    # klabel_set=k_set-0.5;


    ######################
    if band_name=="h1";
        nb=Int(length(eu)/2);#index of band
    elseif band_name=="e1";
        nb=Int(length(eu)/2)+1;#index of band
    end

    
    ########################
    momentum_set,band_state,E_band=get_bandstate(parameters,band_model,N,nK,nb,g1,g3,klabel_set,k1_shift::Vector,k2_shift::Vector,"spin_independent")
    if Mode=="load"
        band_state=band_state_origin;#ensure gauge unchangd 
    end

    parameters_valley2=deepcopy(parameters);
    parameters_valley2["valley"]=1;
    momentum_set_valley2,band_state_valley2,E_band_valley2=get_bandstate(parameters_valley2,band_model,N,nK,nb,g1,g3,klabel_set,k1_shift::Vector,k2_shift::Vector,"spin_independent")


    println("k points:")
    println(momentum_set)

    @show E_band;flush(stdout);
    @show E_band_valley2;
    Fm, Fm_label=get_Fm(N,nK,band_state,momentum_set,klabel_set);

    Fm_valley2, Fm_label_valley2=get_Fm(N,nK,band_state_valley2,momentum_set_valley2,klabel_set);

    q_matrix,q_matrix_label,V_typical=get_coulomb_potential(parameters,N,nK,aM,g1,g3,klabel_set);

    @show band_width=maximum(E_band)-minimum(E_band);
    U_W_ratio=abs(V_typical)/band_width;
    println("ratio between interaction and bandwidth= "*string(U_W_ratio));
    

    ################################

    if Mode=="load"

        q_FBZ=deepcopy(k_set);
        #move k points into FBZ
        for cq=1:size(q_FBZ,1)
            for (ca,cb) in Iterators.product(-10:10, -10:10) 
                kk=[q_FBZ[cq,1]+ca,q_FBZ[cq,2]+cb];
                if determine_in_FBZ(kk)
                    q_FBZ[cq,:]=kk;
                    break;
                end
            end
            @assert determine_in_FBZ(q_FBZ[cq,:]);
        end

        #get q_mirror
        # q_mirror=zeros(Int64,size(q_FBZ,1));
        # for cq=1:size(q_FBZ,1)
        #     ka=q_FBZ[cq,:];
        #     for ci=1:size(q_FBZ,1)
        #         kb=q_FBZ[ci,:];
        #         dk=ka+kb;
        #         dk1=dk[1];
        #         dk2=dk[2];
                
        #         if ((abs.(dk1))<1e-10)&((abs.(dk2))<1e-10)
        #             q_mirror[cq]=ci;
        #             break;
        #         end
        #     end
            
        # end
        # @show q_FBZ
        # @show q_mirror

        # q_keep=ones(Int64,size(q_FBZ,1));
        # for cq=1:size(q_FBZ,1)
        #     kk=q_FBZ[cq,:];
        #     if norm(kk)<1e-10
        #         #origin point:do not keep
        #         continue;
        #     elseif q_mirror[cq]==0
        #         #no mirror point:keep
        #         continue;
        #     else
        #         p2=q_mirror[cq];
        #         if q_keep[cq]*q_keep[p2]==1
        #             q_keep[cq]=0;
        #         end
        #     end
        # end
        # @show q_keep

        #find q equals to multiple of G vectors
        println("Integer q vectors")


        rho_q_coe_1=zeros(ComplexF64,2*nK+1,2*nK+1,size(k_set,1));
        rho_q_coe_2=zeros(ComplexF64,2*nK+1,2*nK+1,size(k_set,1));
        for (ca,cb) in Iterators.product(1:2*nK+1, 1:2*nK+1) 

            Q=[ca-nK-1,cb-nK-1];
            for ck=1:size(k_set,1)
                cond=(abs.(Fm_label[ck,ck,:,:,1].-Q[1]).<1e-10).*(abs.(Fm_label[ck,ck,:,:,2].-Q[2]).<1e-10);
                pos=findall(x->x.==1,cond);
                pos=pos[1];
                @assert sum(abs.(Fm_label[ck,ck,pos[1],pos[2],:]-Q))<1e-10;
                rho_q_coe_1[ca,cb,ck]=Fm[ck,ck,pos[1],pos[2]];

                cond=(abs.(Fm_label[ck,ck,:,:,1].+Q[1]).<1e-10).*(abs.(Fm_label[ck,ck,:,:,2].+Q[2]).<1e-10);
                pos=findall(x->x.==1,cond);
                pos=pos[1];
                @assert sum(abs.(Fm_label[ck,ck,pos[1],pos[2],:]+Q))<1e-10;
                rho_q_coe_2[ca,cb,ck]=Fm[ck,ck,pos[1],pos[2]];

            end
        end
        # @show q_integer
        # @show rho_q_coe_1
        # @show rho_q_coe_2


        
        # @show sum(q_keep)

        k1234_set=Vector{Any}(undef,size(q_FBZ,1));
        V1234_set=Vector{Any}(undef,size(q_FBZ,1));

        #cdag_{k1} c_{k1+q}  cdag_{k2}  c_{k2-q}
        for cq=1:size(q_FBZ,1)
     

            k1234_all_G=Matrix{Any}(undef,2*nK+1,2*nK+1);
            V1234_set_cq=zeros(N^2,2*nK+1,2*nK+1)*im;
            for (ca,cb) in Iterators.product(1:2*nK+1, 1:2*nK+1) 
                dq=q_FBZ[cq,:]+[ca-nK-1,cb-nK-1];

                k1234=zeros(Int,(N^2,4));

                Step=0;

                for c1 in range(1,N)
                    k1_plus_q=k_set[c1,:]+dq;
                    k_label_sum=k_set.-k1_plus_q';
                    k_label_sum_int=round.(k_label_sum);
                    cond=(abs.(k_label_sum[:,1]-k_label_sum_int[:,1]).<1e-10).*(abs.(k_label_sum[:,2]-k_label_sum_int[:,2]).<1e-10);
                    pos=findall(x->x==1,cond);
                    posa=pos[1];
                    #println(pos)

                    for c2 in range(1,N)
                        
                        k2_minus_q=k_set[c2,:]-dq;
                        k_label_sum=k_set.-k2_minus_q';
                        k_label_sum_int=round.(k_label_sum);
                        cond=(abs.(k_label_sum[:,1]-k_label_sum_int[:,1]).<1e-10).*(abs.(k_label_sum[:,2]-k_label_sum_int[:,2]).<1e-10);
                        pos=findall(x->x==1,cond);
                        posb=pos[1];
                        #println(pos)

                                
                        Step=Step+1;
                        k1234[Step,:]=[c1,posa,c2,posb];

                    end
                end

                k1234_set[cq]=k1234;
                k1234_all_G[ca,cb]=k1234;
                @assert norm(k1234_all_G[ca,cb]-k1234_all_G[1,1])==0;
                V1234_new=get_V_k1234_spinup_single_q(nK,k_set,Fm,Fm_label,q_matrix_label,dq,k1234_set[cq]);
                # println(size(V1234_new))
                # println(size(V1234_set_cq))
                V1234_set_cq[:,ca,cb]=V1234_new;
            end 
            V1234_set[cq]=V1234_set_cq;
            
        end
        return nK,k_set,q_FBZ,k1234_set,V1234_set, rho_q_coe_1, rho_q_coe_2
    end

    ##################################
    #not sure if the unit is correct
    # println("not sure if the unit is correct")
    # coe=(1e10)*(1e-3);
    # q_matrix=q_matrix*coe;#to ev
    ##################################
    # println(q_matrix)




        
    #################################


    #cdag_{k1,sigma} cdag_{k2,sigma'} c_{k3,sigma'} c_{k4,sigma}
    Step=0;
    k1234=zeros(Int,(N^4,4));
    for c1 in range(1,N)
        for c2 in range(1,N)
            for c3 in range(1,N)
                for c4 in range(1,N)
                    k_label_sum=sum(klabel_set[[c1,c2],:]-klabel_set[[c3,c4],:],dims=1);
                    k_label_sum_int=round.(k_label_sum);
                    if (abs(k_label_sum[1]-k_label_sum_int[1])<1e-10)&(abs(k_label_sum[2]-k_label_sum_int[2])<1e-10)
                        Step=Step+1;
                        k1234[Step,:]=[c1,c2,c3,c4];
                    end
                end
            end
        end
    end
                        
                        
    k1234=k1234[1:Step,:];

    ########################################
    
    V_k1234_general=get_V_k1234_spinup(nK,klabel_set,Fm,Fm_label,q_matrix_label,q_matrix,k1234);
    V_k1234=V_k1234_general[:];
    if simplify_interaction
        println("simplify interaction terms")
        k1234,V_k1234=reduce_interaction(k1234,V_k1234);
    end


    ##########################################
    global HF_type
    if HF_type=="singleband_singleValley_spinful"
        E_Hatree=Hartree_Fock_correction_singleband_singleValley_spinful(N,nK,klabel_set,Fm_label,Fm,q_matrix_label,q_matrix);
        E_Hatree=E_Hatree[:];
        E_Fock=E_Hatree*0;

        @show E_Hatree;
        @show E_Fock;
    elseif HF_type=="twoValley_spinful_3+nu"
        E_Hatree1=Hartree_Fock_correction_singleband_singleValley_spinful(N,nK,klabel_set,Fm_label,Fm,q_matrix_label,q_matrix);
        E_Hatree1=E_Hatree1[:];
        E_Fock1=E_Hatree1*0;

        E_Hatree2=Hartree_Fock_correction_differentValley_spinful(N,nK,klabel_set, Fm_label,Fm, Fm_label_valley2,Fm_valley2, q_matrix_label,q_matrix);
        E_Hatree2=E_Hatree2[:];
        E_Fock2=E_Hatree2*0;

        @show E_Hatree1;
        @show E_Hatree2;
        @show E_Fock1;
        @show E_Fock2;

        E_Hatree=E_Hatree1+E_Hatree2*2;
        E_Fock=E_Fock1+E_Fock2*2;

    end


    ##########################################
    basis_up=spinless_fermion_basis_1d(N,Np);
    ########################################



    n_E=6;
    occu_set_up,particle_set_up=get_spinless_occu(N,Np,basis_up,length(basis_up));
    particle_set_up=[];

    k_total_set=klabel_set;
    # k_total_set[:,1]=k_total_set[:,1].-k1_shift*(Nup+Ndn-1);
    # k_total_set[:,2]=k_total_set[:,2].-k2_shift*(Nup+Ndn-1);
    ksector_size=zeros(Int64,size(k_total_set,1));

    @tensor Kset_up[:]:=occu_set_up[-1,1]*klabel_set[1,-2];
    occu_set_up=[];
    GC.gc(true);

    ########################
    println("start ED:");flush(stdout);

    #Qn_set=list();
    #Qn_set.append(list([Nup,Ndn]));
    En_set=Vector{Vector}(undef,N);


    # for nsec in range(1,N)
        starting_time=now()
        ksector_size,ksector_ind=get_sector_dim_spinup(N,basis_up,Kset_up,k_total_set,nsec);
        Dim=ksector_size;
        println("dim="*string(Dim));flush(stdout);
        # ksector_ind=ksector_ind[end:-1:1]
        # println(ksector_size)
        # println(ksector_total_ind)
        Now=now();
        Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
        println("time for construct basis: "*string(Time));flush(stdout);
        starting_time=now()
        #basis,H=ED_spinup(N,Nup,external_basis,H_terms)
        spin1=1;
        global break_Hamiltonian,use_HatrrFock;
        if use_HatreeFock
            Hk=kinetic_energy_1band_spinup(N,ksector_ind,E_band+E_Hatree+E_Fock);
        else
            Hk=kinetic_energy_1band_spinup(N,ksector_ind,E_band);
        end
        global use_kinetic_energy;
        println("use_kinetic_energy="*string(use_kinetic_energy));
        if use_kinetic_energy==false
            Hk=Hk*0;
        end

        if ~simplify_interaction
            # V_upup=interaction_1band_spinup(N,ksector_ind, k1234,V_k1234);
        else
            if ~break_Hamiltonian
                #k1234_upup,V_k1234_upup, k1234_dndn,V_k1234_dndn,  k1234_updn,V_k1234_updn
                # V_upup=interaction_1band_spinup(N,ksector_ind, k1234,V_k1234);
            else
                H_set=Vector{SparseMatrixCSC}(undef,0);
                #push!(H_set,Hk);
                push!(H_set,ComplexF32.(Hk));
            
                H_set=interaction_1band_spinup_break_H(H_set,N,ksector_ind, k1234,V_k1234);
            end
        end

        if Dim>2000000
            GC.gc(true);
        end
        if ~break_Hamiltonian
            H=V_upup+Hk;
        else
            #H=Hk+V_upup;
            function apply_H(x::Vector,H_set::Vector{SparseMatrixCSC})
                x1=H_set[1]*x;
                for cc=2:length(H_set)
                    x1=x1+H_set[cc]*x;
                end
                return x1
            end
            F_H(x)=apply_H(x,H_set);
        end

        Now=now();
        Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
        println("time for construct Hamiltonian: "*string(Time));flush(stdout);

        starting_time=now()
        #eu,ev=H.eigh()
        if ~break_Hamiltonian
            if Dim<10
                eu,ev=eigen(Matrix(H));
            else
                #println("compare:")
                #eu,ev=eigs(H,nev=min(n_E,Dim-2),which=:SR)
                eu,ev=eigsolve(H,min(n_E,Dim-2),:SR, ishermitian=true)
            end
        else
            GC.gc(true);
            x_init=ComplexF32.(randn(Dim)*im);x_init=x_init/norm(x_init);
            #x_init=randn(Dim)*im;x_init=x_init/norm(x_init);
            println("test speed of Mv");flush(stdout);
            @time F_H(x_init);
            @time F_H(x_init);
            eu,ev=eigsolve(F_H,x_init,min(n_E,Dim-2),:SR, ishermitian=true)
        end
        @assert norm(imag.(eu))/norm(eu)<1e-12
        eu=real.(eu);
        eu=sort(eu);
        println(eu)
        
        Now=now();
        Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
        println("time for diagonalizing: "*string(Time));flush(stdout);
        
        GC.gc(true);
        # En_set[nsec]=eu;
    # end

    ####################
    #clear memory

    GC.gc(true);

    # matwrite(filenm_, Dict(
    #     "k_total_set" => k_total_set,
    #     "En_set" => En_set
    # ); compress = false)



    #construct disctionary for basis index
    dim=length(ksector_ind);
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(N,ksector_ind[ci]);
    end
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))



    return N,eu[1:n_keep],ev[1:n_keep],k_total_set,ksector_ind,basis_dict,band_state

end
# ##########################

function determine_in_FBZ(kk)
    g1=[1 0];
    g3=[-1/2 sqrt(3)/2];
    
    theta=pi/3;
    R_60=[cos(theta) -sin(theta);sin(theta) cos(theta)];
    corner1=[0.5 0.5/sqrt(3)];
    corners=zeros(7,2);
    for cc=1:7
        corners[cc,:]=R_60^(cc-1)*corner1';
    end
    
    g_set=zeros(6,2);
    for cc=1:size(g_set,1)
        g_set[cc,:]=R_60^(cc-1)*g1';
    end
    # klabel_set=k_set;
    
    nK=5;
    K_shift=-nK:nK;
    ca_cb=zeros((2*nK+1)*(2*nK+1),2);
    step=1;
    for ca=1:2*nK+1
        for cb=1:2*nK+1
            ca_cb[step,:]=[ca cb];
            step=step+1;
        end
    end
    
    qx_set=0;
    qy_set=0;
 
    in_FBZ=false
    
    k0=(kk[1])*g1+(kk[2])*g3;
    distance=norm(k0);
    distance0=norm(g1)/2;

    ds0=norm(k0);
    dss=zeros(6,1);
    for dd=1:length(dss)
        dss[dd]=sqrt((k0[1]-g_set[dd,1])^2+(k0[2]-g_set[dd,2])^2);
    end
    qx_set=k0[1];
    qy_set=k0[2];
    if ds0<=minimum(dss)+0.001;
        in_FBZ=true
    end
    return in_FBZ
end




function main2(waa_ratio,Mode,save_K_ind::Number,K_ind)



band_name="e1";#"e1","h1"



parameters=Dict{String,Any}();
stack="ABA";
merge!(parameters,Dict{String,Any}("stack"=>stack));
if stack=="ABA"
    theta_ = 1.5;
    Phi=[0,1,-1]*(2*pi/3);
elseif stack=="AAA"
    theta_ = 0.68;
    Phi=[0,0,0]*(2*pi/3);
end

merge!(parameters,Dict{String,Any}("theta_"=>theta_));
merge!(parameters,Dict{String,Any}("Phi"=>Phi));
merge!(parameters,Dict{String,Any}("a0"=>2.46e-10));

merge!(parameters,Dict{String,Any}("wab"=>110));
merge!(parameters,Dict{String,Any}("waa"=>waa_ratio*parameters["wab"]));
merge!(parameters,Dict{String,Any}("v"=>1*(10^6)*(1.054571817e-34)*(6.24150636e21)));
merge!(parameters,Dict{String,Any}("valley"=>-1));

#sublattice potential
# AB_potential1=[20 -20];
# AB_potential2=[0 0];
# AB_potential3=[20 -20];
# AB_potential1=[-4 4];
# AB_potential2=[-6 6];
# AB_potential3=[-8 8];
AB_potential1=[-8 8];
AB_potential2=[-8 8];
AB_potential3=[-8 8];
merge!(parameters,Dict{String,Any}("AB_potential"=>[AB_potential1;AB_potential2;AB_potential3]));#ev

merge!(parameters,Dict{String,Any}("screen_type"=>"No"));#"No","two_gate","single_gate","yukawa"
merge!(parameters,Dict{String,Any}("kappa"=>0));#ev, yukawa screening of coulomb interaction
merge!(parameters,Dict{String,Any}("d_sc"=>3));#nm, screen length, gate distance

merge!(parameters,Dict{String,Any}("epsilon"=>4));

band_model=band_single_k_Helical_trilayer;#
global use_kinetic_energy,break_Hamiltonian,band_model,simplify_interaction,use_HatreeFock
use_kinetic_energy=true;
break_Hamiltonian=true;#break big Hamiltonian into sections so that construction sparse Hamiltonian is possible
simplify_interaction=true;
@show use_HatreeFock=true;

global HF_type
@show HF_type="twoValley_spinful_3+nu";#"twoValley_spinful_3+nu", "singleband_singleValley_spinful"

Np=8;


# V1=[3,0];# 15 sites
# V2=[0,5];# 15 sites

V1=[4,0];# 16 sites
V2=[0,4];# 16 sites
# V1=[3,-3];# 18 sites
# V2=[0,6];# 18 sites
# V1=[4,-4];# 20 sites
# V2=[0,5];# 20 sites
# V1=[1,4];# 24E sites
# V2=[5,-4];# 24E sites
# V1=[2,4]; #26 sites
# V2=[0,13]; #26 sites
# V1=[6,-3];#integer vector, 27
# V2=[3,-6];#integer vector, 27
# V1=[-2,6];#integer vector: 28 sites
# V2=[4,2];#integer vector: 28 sites
# V1=[3,3];#integer vector: 30D sites
# V2=[5,-5];#integer vector: 30D sites
# V1=[2,4];#integer vector: 32 sites
# V2=[6,-4];#integer vector: 32 sites
# V1=[6,0];# 36 sites
# V2=[0,6];# 36 sites

flux_1=0;
flux_2=0;

#################################
#define momentum sectors and degeneracy of ground states
# K_id_set=[24+1,25+1];#N26
# n_keep_set=[1,1];#N26

# K_id_set=[16+1,17+1,18+1,19+1];#N20
# n_keep_set=[1,2,1,2];#N20


# K_id_set=[21+1];#N27,Np=9
# n_keep_set=[3];#N27,Np=9

# K_id_set=[6+1,18+1,21+1];#N27,Np=18
# n_keep_set=[1,1,1];#N27,Np=18

# K_id_set=[4+1,14+1];#N32
# n_keep_set=[4,2];#N32

# K_id_set=[20+1,22+1];#N24A MR
# n_keep_set=[2,4];#N24A MR

# K_id_set=[2+1,17+1];#N24A CFL
# n_keep_set=[2,2];#N24A CFL

# K_id_set=[1+1,11+1,15+1];#N28 MR
# n_keep_set=[2,2,2];#N28 MR

# K_id_set=[4+1,14+1];#N32 MR
# n_keep_set=[4,2];#N32 MR

# K_id_set=[0+1,16+1,26+1];#N36, Np=24
# n_keep_set=[1,1,1];#N36, Np=24

# K_id_set=[0+1];#N36, Np=12
# n_keep_set=[3];#N36, Np=12

# K_id_set=[0+1,1+1,11+1,15+1];#N28 Np=14
# n_keep_set=[3,4,4,4];#N28 Np=14

K_id_set=[0+1,2+1,8+1,10+1];#N16 Np=8
n_keep_set=[3,2,2,2];#N16 Np=8

#define variable that stores ground states, energies, basis dictionaries


#save ground states after diagonalizing Hamiltonian

    if Mode=="save"
        ck=save_K_ind;
        @show  K_id_set[ck],n_keep_set[ck]
        N,eus,evs,k_total_set,ksector_ind,basis_dict,band_state=main(band_name,parameters,Np, V1, V2, flux_1,flux_2, simplify_interaction, K_id_set[ck],n_keep_set[ck],[],Mode);
        # eu_set[ck]=eus;
        # ev_set[ck]=evs;
        # ksector_ind_set[ck]=ksector_ind;
        # basis_dict_set[ck]=basis_dict;

        filenm_ = "theta"*string(parameters["theta_"])*"_waa"*string(waa_ratio)*"_epsilon"*string(parameters["epsilon"])*"_Np"*string(Np)*"_sec"*string(K_id_set[ck])*".jld2";
        jldsave(filenm_;N,eus,evs,k_total_set,ksector_ind,basis_dict,band_state);


    elseif Mode=="load"
        ck=K_ind;
        filenm_ = "theta"*string(parameters["theta_"])*"_waa"*string(waa_ratio)*"_epsilon"*string(parameters["epsilon"])*"_Np"*string(Np)*"_sec"*string(K_id_set[ck])*".jld2";
        data=load(filenm_);
        band_state_origin=data["band_state"];#load band_state without rediagonalization, to ensure unchanged gauge of single particle states


        nK,k_set,q_FBZ,k1234_set,V1234_set, rho_q_coe_1, rho_q_coe_2=main(band_name,parameters,Np, V1, V2, flux_1,flux_2, simplify_interaction, 0,0,band_state_origin,Mode);

        


        N=data["N"]
        eus=data["eus"];
        evs=data["evs"];
        ksector_ind=data["ksector_ind"];
        basis_dict=data["basis_dict"];
        k_total_set=data["k_total_set"];
        
        data=[];


        ###############
        nk=zeros(size(k_set,1),n_keep_set[ck])*im;
        #momentum distribution
        println("compute nk");flush(stdout);

        for cp=1:length(ind_start:size(k_set,1))
            println(cp+ind_start-1)
            ind_,value=kinetic_energy_1band_spinup_ob(N,ksector_ind, k_set,cp+ind_start-1,evs[1:n_keep_set[ck]])  
            nk[ind_,:]=value;
        end
            

        println(nk);flush(stdout);
        ###############




        Sq_set=zeros(size(q_FBZ,1),n_keep_set[ck],2*nK+1,2*nK+1)*im;
        ob_terms_set=Vector{Any}(undef,size(q_FBZ,1))
        for cq=1:size(q_FBZ,1)
            
            
            @show cq;flush(stdout);
            #v0s=evs[1:n_keep_set[ck]];


            v0s=zeros(n_keep_set[ck],length(ksector_ind))*im;
            for cs=1:n_keep_set[ck]
                v0s[cs,:]=evs[cs];
            end

            # @show k1234_set[cq]
            # @show V1234_set[cq]
            #########################
            
            @show N_terms=size(k1234_set[cq],1);

            ob_terms=zeros(N_terms,n_keep_set[ck])*im;
            for cp=1:N_terms

                ind_,value=interaction_1band_spinup_qset_distributed(N,ksector_ind, group_ind[cp,1]:group_ind[cp,2], k1234_set[cq][group_ind[cp,1]:group_ind[cp,2],:],v0s);
                ob_terms[ind_,:]=value;
                
            end
            
            @tensor Sq_[:]:=ob_terms[1,-1]*V1234_set[cq][1,-2,-3];
            Sq_set[cq,:,:,:]=Sq_/N;
            ob_terms_set[cq]=ob_terms;
            GC.gc(true);
            #########################
            # Ham_q=interaction_1band_spinup_single_q_distributed(N,ksector_ind, k1234_set[cq], V1234_set[cq],v0);
            # for cs=1:n_keep_set[ck]
            #     v0=evs[cs];
            #     # @show size(Ham_q)
            #     # @show size(v0)
            #     if size(Ham_q,1)<length(v0)
            #         Sq_set[cq,cs]=v0[1:size(Ham_q,1)]'*Ham_q*v0[1:size(Ham_q,1)]/N;
            #     else
            #         Sq_set[cq,cs]=v0'*Ham_q*v0/N;
            #     end
            # end


            
            
        end
        #println(Sq_set)

        ################################################

        
        # Ham_q1=kinetic_energy_1band_spinup(N,ksector_ind, rho_q_coe_1[cq,:]);
        # Ham_q2=kinetic_energy_1band_spinup(N,ksector_ind, rho_q_coe_2[cq,:]);
        # for cs=1:n_keep_set[ck]
        #     v0=evs[cs];
        #     rho_q1[cq,cs]=(v0'*Ham_q1*v0)/sqrt(N);
        #     rho_q2[cq,cs]=(v0'*Ham_q2*v0)/sqrt(N);

        # end


        #evaluate <rho(q)> and <rho(-q)>
        l1,l2,l3=size(rho_q_coe_1);
        rho_q_coe_1=reshape(rho_q_coe_1,(l1*l2,l3));
        rho_q_coe_2=reshape(rho_q_coe_2,(l1*l2,l3));

        rho_q1=zeros(l1*l2,n_keep_set[ck])*im;
        rho_q2=zeros(l1*l2,n_keep_set[ck])*im;



        for cp=1:length(ind_start:size(rho_q_coe_1,1))
            println(cp+ind_start-1)
            ind_,value1,value2=kinetic_energy_1band_spinup_Eband(N,ksector_ind,cp+ind_start-1,evs[1:n_keep_set[ck]], rho_q_coe_1[cp+ind_start-1,:], rho_q_coe_2[cp+ind_start-1,:])  
            rho_q1[ind_,:]=value1/sqrt(N);
            rho_q2[ind_,:]=value2/sqrt(N);
        end


        filenm_ = "FBZ_paircorrel_theta"*string(parameters["theta_"])*"_waa"*string(waa_ratio)*"_epsilon"*string(parameters["epsilon"])*"_Np"*string(Np)*"_sec"*string(K_id_set[K_ind])*".mat";
        if K_ind>1#only save these information for the first sector
            k1234_set=0;
            V1234_set=0;
        end
        matwrite(filenm_, Dict(
            "k1234_set"=>k1234_set,
            "V1234_set"=>V1234_set,
            "nk"=>nk,
            "ob_terms_set"=>ob_terms_set,
            "nK"=>nK,
            "k_set" =>k_set,
            "q_FBZ" => q_FBZ,
            "Sq_set" =>Sq_set,
            "rho_q1" =>rho_q1,
            "rho_q2"=>rho_q2,
            "eu_set"=>eus
        ); compress = false)        
    end
end


#algorithm:
#step 1: compute all structure factor at point q in FBZ, where contribution from q+G,q-G, q+2G... are included
#step 2: compute <n(q)> for q=0,G,-G,2G,-2G,.....


# println("Distributed")

# @distributed for Nsec in 1:6
#     @show Nsec;
#     main2(Nsec)
# end



waa_ratio=0.0;


@show save_K_ind=4;# from 1 to number of ground state K sectors
@show K_ind=save_K_ind;# from 1 to number of ground state K sectors

@show Mode="save";#"save", "load"
main2(waa_ratio,Mode, save_K_ind,K_ind)




