using Revise, TensorKit
using JLD2,MAT
using KrylovKit
using JSON
using Random
using Dates
using SparseArrays,Combinatorics
#using Arpack
using KrylovKit

cd(@__DIR__)
include("ED_lib.jl")
include("model_lib.jl")

####################
import LinearAlgebra.BLAS as BLAS
n_cpu=4;
BLAS.set_num_threads(n_cpu);
println("number of cpus: "*string(BLAS.get_num_threads()))
Base.Sys.set_process_title("C"*string(n_cpu)*"_ED")
pid=getpid();
println("pid="*string(pid));
####################
function main(band_name,parameters,Np,V1,V2,flux_1,flux_2,simplify_interaction,n_levels,K_ind)
    println(parameters)
    global band_model
    
    #############################################
    theta_,theta,aM,ktheta,waa,wab,v,valley,AB_potential,Phi,kappa=get_parameters(parameters);
    #############################################
    println("theta="*string(theta_))
    println("Np="*string(Np))

    
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



    n_E=10;
    occu_set_up,particle_set_up=get_spinless_occu(N,Np,basis_up,length(basis_up));
    particle_set_up=[];
    
    k_total_set=klabel_set;
    # k_total_set[:,1]=k_total_set[:,1].-k1_shift*(Nup+Ndn-1);
    # k_total_set[:,2]=k_total_set[:,2].-k2_shift*(Nup+Ndn-1);
    ksector_size=zeros(Int64,size(k_total_set,1));

    @tensor Kset_up[:]:=occu_set_up[-1,1]*klabel_set[1,-2];
    occu_set_up=[];


    ########################
    println("start ED:");flush(stdout);

    #Qn_set=list();
    #Qn_set.append(list([Nup,Ndn]));
    En_set=Vector{Vector}(undef,N);


    nsec =K_ind;
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
        
        En_set[nsec]=eu;
    


        function get_up_occu_K(N,ksector_ind)
            particle_occu_set_K=Matrix{Int8}(undef,length(ksector_ind),Np);
            for cc=1:length(ksector_ind)
                Str=get_bit(N,ksector_ind[cc]);
                Str_up=Str[1:N];
                Str_up=split(Str_up, "");
                pos=findall(x->x.=="1",Str_up)
                particle_occu_set_K[cc,:]=pos;
            end
            return particle_occu_set_K
        end
        particle_occu_set_K=get_up_occu_K(N,ksector_ind);

        GC.gc(true);


    return klabel_set[K_ind,:],eu[1:n_levels],ev[1:n_levels], band_state, particle_occu_set_K

end
# ##########################

function mainmain(n_level,K_ind_set,L_flx,waa_ratio,filenm_end)

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



Np=18;


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
# V1=[3,3];#integer vector: 30 sites
# V2=[5,-5];#integer vector: 30 sites
# V1=[2,4];#integer vector: 32 sites
# V2=[6,-4];#integer vector: 32 sites



@show band_name;#"h1","e1"
@show V1
@show V2
#############################
# n_levels=1;

flx_set=range(0,L_flx,L_flx+1)/L_flx;
# K_ind_set=Int.([1,2,3]);#index for total momentum sector
K_set=Vector{Any}(undef,length(K_ind_set));

nE=10;
Eset=Array{Float64}(undef,L_flx+1,L_flx+1,nE,length(K_ind_set));
F_grid=Array{Float64}(undef,L_flx,L_flx,length(K_ind_set))*im;


for K_ind in eachindex(K_ind_set)

    band_state_set=Matrix{Any}(undef,L_flx+1,L_flx+1);
    particle_occu_set_K=nothing;
    
    psi_set=Array{Any}(undef,L_flx+1,L_flx+1,n_levels);
    Ux_grid=Matrix{Number}(undef,L_flx,L_flx+1);
    Uy_grid=Matrix{Number}(undef,L_flx+1,L_flx);
    for fl1 in range(1,L_flx+1), fl2 in range(1,L_flx+1)

        flux_1=flx_set[fl1];
        flux_2=flx_set[fl2];

        KK,eus,evs,band_state, particle_occu_set_K=main(band_name,parameters,Np,V1,V2, flux_1,flux_2, simplify_interaction, nE,K_ind_set[K_ind]);

        K_set[K_ind]=KK;
        psi_set[fl1,fl2,1:n_levels]=evs[1:n_levels];
        Eset[fl1,fl2,1:nE,K_ind]=eus[1:nE];


        band_state_set[fl1,fl2]=band_state;

        
    end
    

    

    ovx_band=nothing;
    ovy_band=nothing;

    for fl1 in range(1,L_flx), fl2 in range(1,L_flx+1)  
        
        band_state1=band_state_set[fl1,fl2];
        band_state2=band_state_set[fl1+1,fl2];
        ovx_band=zeros(size(band_state1,1))*im;
        for c1=1:size(band_state1,1)
            ovx_band[c1]=dot(band_state1[c1,:,:,:,:,:][:],band_state2[c1,:,:,:,:,:][:]);
        end

        band_correction=zeros(size(particle_occu_set_K,1))*im;
        for cc=1:size(particle_occu_set_K,1)
            band_correction[cc]=prod(ovx_band[particle_occu_set_K[cc,:]])
        end

        ux=zeros(n_levels,n_levels)*im;
        for e1=1:n_levels
            for e2=1:n_levels
                ux[e1,e2]=dot(psi_set[fl1,fl2,e1],band_correction.*psi_set[fl1+1,fl2,e2]);
            end
        end
        ux=det(ux);
        Ux_grid[fl1,fl2]=ux/abs(ux);
    end

    for fl1 in range(1,L_flx+1), fl2 in range(1,L_flx)   
        band_state1=band_state_set[fl1,fl2];
        band_state2=band_state_set[fl1,fl2+1];
        ovy_band=zeros(size(band_state1,1))*im;
        for c1=1:size(band_state1,1)
            ovy_band[c1]=dot(band_state1[c1,:,:,:,:,:][:],band_state2[c1,:,:,:,:,:][:]);
        end
        
        band_correction=zeros(size(particle_occu_set_K,1))*im;
        for cc=1:size(particle_occu_set_K,1)
            band_correction[cc]=prod(ovy_band[particle_occu_set_K[cc,:]])
        end
        
        uy=zeros(n_levels,n_levels)*im;
        for e1=1:n_levels
            for e2=1:n_levels
                uy[e1,e2]=dot(psi_set[fl1,fl2,e1],band_correction.*psi_set[fl1,fl2+1,e2]);
            end
        end
        uy=det(uy);
        Uy_grid[fl1,fl2]=uy/abs(uy);
    end

    for fl1 in range(1,L_flx)
        for fl2 in range(1,L_flx)
            F_grid[fl1,fl2,K_ind]=Ux_grid[fl1,fl2]*Uy_grid[fl1+1,fl2]/Ux_grid[fl1,fl2+1]/Uy_grid[fl1,fl2];
            F_grid[fl1,fl2,K_ind]=log(F_grid[fl1,fl2,K_ind]);
        end
    end
    
end

#C=sum(F_grid[:,:,:])/2/pi;

filenm_ = "Chern_theta"*string(parameters["theta_"])*"_waa"*string(waa_ratio)*"_epsilon"*string(parameters["epsilon"])*"_Np"*string(Np)*"_Lflx_"*string(L_flx)*filenm_end;
# jldsave(filenm_*".jld2";flx_set,K_ind_set,K_set,Eset,F_grid);
matwrite(filenm_*".mat", Dict(
    "Eset" => Eset,
    "F_grid" => F_grid,
    "flx_set" =>Vector(flx_set),
    "K_set"=>K_set,
    "K_ind_set"=>K_ind_set
); compress = false)


end

L_flx=6;
waa_ratio=0.4;

n_levels=1;
K_ind_set=Int.([6+1]);
mainmain(n_levels,K_ind_set,L_flx,waa_ratio,"A");

# n_levels=1;
# K_ind_set=Int.([18+1]);
# mainmain(n_levels,K_ind_set,L_flx,waa_ratio,"B");

# n_levels=1;
# K_ind_set=Int.([21+1]);
# mainmain(n_levels,K_ind_set,L_flx,waa_ratio,"C");





