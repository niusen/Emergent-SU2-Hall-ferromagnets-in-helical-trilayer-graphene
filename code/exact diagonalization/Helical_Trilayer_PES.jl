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
####################
import LinearAlgebra.BLAS as BLAS
@show n_blas_thread=1;
BLAS.set_num_threads(n_blas_thread);
@show nthread=Threads.nthreads();#number of thread can only be set in job.sh file
Base.Sys.set_process_title("C"*string(nthread)*"_ED")
pid=getpid();
println("pid="*string(pid));
####################
@everywhere function main(band_name,parameters,Np,V1,V2,flux_1,flux_2,simplify_interaction,nsec,n_keep)
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



    return N,eu[1:n_keep],ev[1:n_keep],k_total_set,ksector_ind,basis_dict

end
# ##########################

function main2(waa_ratio,Mode,Np,NA, save_K_ind::Number, PES_Nsec_set)



band_name="e1";#"e1","h1"



parameters=Dict{String,Any}();
stack="ABA";
merge!(parameters,Dict{String,Any}("stack"=>stack));
if stack=="ABA"
    theta_ = 1.44;
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
AB_potential1=[-4 4];
AB_potential2=[-6 6];
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




# V1=[3,0];# 15 sites
# V2=[0,5];# 15 sites

# V1=[4,0];# 16 sites
# V2=[0,4];# 16 sites
# V1=[3,-3];# 18 sites
# V2=[0,6];# 18 sites
# V1=[4,-4];# 20 sites
# V2=[0,5];# 20 sites
# V1=[4,-4];# 24 sites
# V2=[0,6];# 24 sites
# V1=[2,4]; #26 sites
# V2=[0,13]; #26 sites
V1=[6,-3];#integer vector, 27
V2=[3,-6];#integer vector, 27
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

K_id_set=[6+1,18+1,21+1];#N27,Np=18
n_keep_set=[1,1,1];#N27,Np=18

# K_id_set=[4+1,14+1];#N32
# n_keep_set=[4,2];#N32

# K_id_set=[6+1,13+1,20+1];#N24A Np 8
# n_keep_set=[1,1,1];#N24A Np 8

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

#define variable that stores ground states, energies, basis dictionaries
eu_set=Vector{Any}(undef,length(K_id_set));
ev_set=Vector{Any}(undef,length(K_id_set));
ksector_ind_set=Vector{Any}(undef,length(K_id_set));
basis_dict_set=Vector{Any}(undef,length(K_id_set));
k_total_set=[];#store all K points
N=[];#size of cluster

#save ground states after diagonalizing Hamiltonian

    if Mode=="save"
        ck=save_K_ind;
        @show  K_id_set[ck],n_keep_set[ck]
        N,eus,evs,k_total_set,ksector_ind,basis_dict=main(band_name,parameters,Np, V1, V2, flux_1,flux_2, simplify_interaction, K_id_set[ck],n_keep_set[ck]);
        # eu_set[ck]=eus;
        # ev_set[ck]=evs;
        # ksector_ind_set[ck]=ksector_ind;
        # basis_dict_set[ck]=basis_dict;

        filenm_ = "theta"*string(parameters["theta_"])*"_waa"*string(waa_ratio)*"_epsilon"*string(parameters["epsilon"])*"_Np"*string(Np)*"_sec"*string(K_id_set[ck])*".jld2";
        jldsave(filenm_;N,eus,evs,k_total_set,ksector_ind,basis_dict);


    elseif Mode=="load"
        for ck=1:length(K_id_set)

            
            filenm_ = "theta"*string(parameters["theta_"])*"_waa"*string(waa_ratio)*"_epsilon"*string(parameters["epsilon"])*"_Np"*string(Np)*"_sec"*string(K_id_set[ck])*".jld2";
            data=load(filenm_);
            N=data["N"]
            eu_set[ck]=data["eus"];
            ev_set[ck]=data["evs"];
            ksector_ind_set[ck]=data["ksector_ind"];
            basis_dict_set[ck]=data["basis_dict"];
            k_total_set=data["k_total_set"];
            data=[];
        end
        GC.gc(true);

        NB=Np-NA;
        @assert NB>0;
        @assert NA>0;

        PES_set=Vector{Any}(undef,N);
        starting_time=now()

        @assert maximum(PES_Nsec_set)<=N
        # Threads.@threads for nsecA in PES_Nsec_set;#1:N,  choose momentum sector for subsystem A
        for nsecA in PES_Nsec_set;#1:N,  choose momentum sector for subsystem A
            println("nsecA= "*string(nsecA));flush(stdout);
    
            KA=k_total_set[nsecA,:];
    
    
            ksector_ind_A,basis_dict_A,occu_A,particle_set_A=get_basis(N,NA,k_total_set,nsecA);
            @show dimA=length(ksector_ind_A);flush(stdout);
            rhoA=ComplexF32.(zeros(dimA,dimA)*im);
    
            for ck=1:length(K_id_set);#iterate ground state, 1:length(K_id_set)     ,  choose momentum sector for global state
                # ksector_ind_total,basis_dict_total,occu_total,particle_set_total=get_basis(N,NA+NB,k_total_set,K_id_set[ck]);
                println("iterate ck");flush(stdout);
    
    
                Ktotal=k_total_set[K_id_set[ck],:];
                KB=Ktotal-KA;
                jg=(abs.(k_total_set[:,1].-(mod.(KB[1],1))).<1e-10).*(abs.(k_total_set[:,2].-(mod.(KB[2],1))).<1e-10)
                nsecB=findall(x->x==1,jg);#get momentum sector for subsystem B
                @assert length(nsecB)==1;
                nsecB=nsecB[1];
    
                cs=1;
                E=eu_set[ck][cs];
                psi=ev_set[ck][cs];
                ksector_ind=ksector_ind_set[ck];
                basis_dict=basis_dict_set[ck];
                occu_total,particle_set=get_occu(N,Np,ksector_ind);
    
                dim_total=length(psi);
    
                
                #@show varinfo()
                ksector_ind_B,basis_dict_B,occu_B,particle_set_B=get_basis(N,NB,k_total_set,nsecB);
                @show dimB=length(ksector_ind_B);flush(stdout);#dimension of B
    
                for cs=1:n_keep_set[ck];#1:n_keep_set[ck], iterate ground state degeneracy
                    println("iterate cs");flush(stdout);
    
    
    
                    E=eu_set[ck][cs];
                    psi=ev_set[ck][cs];
    
                    psi_matrix=zeros(ComplexF32, dimA,dimB);#reshape vector psi to a matrix
                    
                    nthread=Threads.nthreads();
    
                    @show dimA_section=Int(ceil(dimA/nthread));
                    dimA_devided_=1:dimA_section:dimA;
                    dimA_devided=zeros(Int64,length(dimA_devided_),2);
                    dimA_devided[1,1]=1;
                    dimA_devided[1,2]=dimA_devided_[2]-1;
                    for ccc=2:length(dimA_devided_)-1
                        dimA_devided[ccc,1]=dimA_devided_[ccc];
                        dimA_devided[ccc,2]=dimA_devided_[ccc+1]-1;
                    end
                    dimA_devided[end,1]=dimA_devided_[end];
                    dimA_devided[end,2]=dimA;
    
                    println("iterate dimA");flush(stdout);
    
                    Threads.@threads for ite_section=1:size(dimA_devided,1)####
                        @show ite_section;
                        @show Threads.threadid();flush(stdout);
                        store_vec=zeros(ComplexF32, dimA_devided[ite_section,2]-dimA_devided[ite_section,1]+1,    dimB);
                        for c1=dimA_devided[ite_section,1]:dimA_devided[ite_section,2] #iterate space of A
                            if mod(c1,100)==0
                                @show c1;flush(stdout);
                            end
                            idA=ksector_ind_A[c1];
                            strA=get_bit(N,idA);
                            vecA=occu_A[c1,:];
                            idA_dense=basis_dict_A[strA];
    
                            occu_left=occu_total.-vecA';
                            occu_=sum(abs.(occu_left),dims=2);
                            pos=findall(x->x==NB,occu_);
    
                            
                            for id_total_dense in pos #iterate found global configs
                                
                                vecB=occu_left[id_total_dense[1],:];
                                strB=prod(string.(vecB));
                                idB_dense=basis_dict_B[strB];
    
    
                                coe=psi[id_total_dense[1]];
    
                                #get sign
                                particle_A=particle_set_A[c1,:];
                                sign_total=0;
                                for cn=1:NA
                                    sign_total=sign_total+sum(vecB[1:particle_A[cn]-1]);
                                end
    
                                store_vec[c1-dimA_devided[ite_section,1]+1,idB_dense]=ComplexF32(coe*(-1)^sign_total);
                                #psi_matrix[idA_dense,idB_dense]=ComplexF32(coe*(-1)^sign_total);
                                    
                                
                            end
                            
                        end
                        psi_matrix[dimA_devided[ite_section,1]:dimA_devided[ite_section,2],:]=store_vec;
                        println("finish one section of dimA");flush(stdout);
                        # GC.gc(true);
    
                    end
    
    
                    GC.gc(true);
                    rhoA=rhoA+psi_matrix*psi_matrix';
                    psi_matrix=nothing;
    
                end
                GC.gc(true);
            end
            @show typeof(rhoA)
            rhoA=rhoA/binomial(Np,NA);
            rhoA=rhoA/sum(n_keep_set);#average over degenerate ground states
    
            eu,ev=eigen(rhoA);
    
            PES=eu;
    
            filenm_ = "PES_theta"*string(parameters["theta_"])*"_waa"*string(waa_ratio)*"_epsilon"*string(parameters["epsilon"])*"_Np"*string(Np)*"_NA_"*string(NA)*"_sec"*string(nsecA)*".mat";
            matwrite(filenm_, Dict(
                "PES" => PES,
                "eu_set"=>eu_set
            ); compress = false)
    
            Now=now();
            Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
            println("time consumed: "*string(Time));flush(stdout);
        end
     
    end
end





waa_ratio=0.4;


@show Mode="save";#"save", "load"
@show NA=4;
@show Np=18;
@show PES_Nsec_set=1:27;


@show save_K_ind=1;# from 1 to number of ground state K sectors
main2(waa_ratio,Mode, Np,NA, save_K_ind, PES_Nsec_set)

@show save_K_ind=2;# from 1 to number of ground state K sectors
main2(waa_ratio,Mode, Np,NA, save_K_ind, PES_Nsec_set)

@show save_K_ind=3;# from 1 to number of ground state K sectors
main2(waa_ratio,Mode, Np,NA, save_K_ind, PES_Nsec_set)



