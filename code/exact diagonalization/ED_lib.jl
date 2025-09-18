
using LinearAlgebra
using TensorKit




function get_bit(L::Int,x::Int)
    Str=bitstring(x)
    y=Str[end-L+1:end];
    return y
end

function spinless_fermion_basis_1d(N,Nf)
    Dim=binomial(N,Nf);
    combs=Combinatorics.combinations(1:N, Nf);
    basis=Vector{Int64}(undef,Dim);
    abc=2;
    Step=1;
    for cc in combs
        digits=(cc.-1);
        basis[Step]=sum(abc.^digits)
        Step=Step+1
    end

    return basis
end




function get_spinless_occu(N,Np,states,siz)
    particle_set=zeros(Int8,siz,Np);
    occu_set=zeros(Int8,siz,N);
    for cc in range(1,siz)
        Str=get_bit(N,states[cc]);
        #occu_set[cc,:] = np.fromstring(Str,'u1') - ord('0')
        #occu_set[cc,:] = np.array(list(Str), dtype=int)
        Str_=split(Str, "");
        occu_set[cc,:]=[parse(Int8, ss) for ss in Str_]
        pos=findall(x->x.=="1",Str_)
        particle_set[cc,:]=pos
    end

    return occu_set,particle_set
end



function find_element_spinful(N, states_up, states_dn, Kdiff_bool, ksector_total_ind, ksector_ind, sizup, sizdn)
    Step = 0
    for c1 in 1:sizup
        for c2 in 1:sizdn
            if Kdiff_bool[c1, c2]
                ksector_ind[Step+1, :] = [c1, c2]
                ksector_total_ind[Step+1] = states_up[c1] * 2^N + states_dn[c2]
                Step += 1
            end
        end
    end
    return ksector_total_ind, ksector_ind
end

function get_sector_dim_spinful(N, states_up, states_dn, Kset_up, Kset_dn, k_total_set, sector_id)
    c1 = sector_id
    k_total = k_total_set[c1, :]
    sizup = size(Kset_up, 1)
    sizdn = size(Kset_dn, 1)

    @tensor K_set_a[:]:=Kset_up[-1,-3]*ones(sizdn)[-2]
    @tensor K_set_b[:]:=ones(sizup)[-1]*Kset_dn[-2,-3]
    K_set = K_set_a+K_set_b 
    @tensor Kdiff_set[:]:=ones(sizup, sizdn)[-1,-2]*k_total[-3]
    Kdiff_set = K_set - Kdiff_set
    Kdiff_int_set = round.(Kdiff_set)
    Kdiff_bool = abs.(Kdiff_set - Kdiff_int_set) .<= 1e-10
    Kdiff_bool = Kdiff_bool[:, :, 1] .* Kdiff_bool[:, :, 2]
    ksector_size = sum(Kdiff_bool)

    ksector_total_ind = zeros(Int, ksector_size)
    ksector_ind = zeros(Int, (ksector_size, 2))
    ksector_total_ind, ksector_ind = find_element_spinful(N, states_up, states_dn, Kdiff_bool, ksector_total_ind, ksector_ind, sizup, sizdn)
    return ksector_size, ksector_total_ind, ksector_ind
end
###################


function find_element_spinup(N, states_up, Kdiff_bool, ksector_ind, sizup)
    Step = 0
    for c1 in 1:sizup
        if Kdiff_bool[c1]
            ksector_ind[Step+1] = states_up[c1]
            Step += 1
        end
    end
    return ksector_ind
end



function get_sector_dim_spinup(N, states_up, Kset_up, k_total_set, sector_id)
    c1 = sector_id
    k_total = k_total_set[c1, :]
    sizup = size(Kset_up, 1)
    K_set = Kset_up
    @tensor k_total[:]:=ones(sizup)[-1]*k_total[-2]
    Kdiff_set = K_set -k_total
    Kdiff_int_set = round.(Kdiff_set)
    Kdiff_bool = abs.(Kdiff_set .- Kdiff_int_set) .<= 1e-10
    Kdiff_bool = Kdiff_bool[:, 1] .* Kdiff_bool[:, 2]
    ksector_size = sum(Kdiff_bool)

    ksector_ind = zeros(Int, ksector_size)
    ksector_ind = find_element_spinup(N, states_up, Kdiff_bool, ksector_ind, sizup)
    return ksector_size, ksector_ind
end

##############################

function interaction_1band_spinup(N,manybody_ind, k1234,V_k1234)
    
    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind);
    @assert dim<(2^32); #use Int32 to save index 
    if length(size(V_k1234))==3
        V_k1234=V_k1234[:,spin1,spin2];
    elseif (length(size(V_k1234))==2)&&(size(V_k1234,2)==1)
        V_k1234=V_k1234;
    elseif length(size(V_k1234))==1
        V_k1234=V_k1234;
    else
        error("unknown case")
    end
    #construct disctionary for basis index
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(N,manybody_ind[ci]);
    end
    # basis_dict=Dict(zip(ind_str,manybody_ind))
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))
   
    #ind_empty_Complex=zeros(ComplexF64,dim*N);
    ind_empty_Complex=zeros(ComplexF32,dim*N);
    ind_empty_Int=zeros(Int32,dim*N);
 
    x_ind=deepcopy(ind_empty_Int);
    y_ind=deepcopy(ind_empty_Int);
    value_ind=deepcopy(ind_empty_Complex);
    N_terms=size(k1234,1);
    
    Step_=0;
    for c1=1:dim
        state1=get_bit(N,manybody_ind[c1]);
        for ct=1:N_terms
            (k1,k2,k3,k4)=k1234[ct,:];
            
            
            if ( (k1==k2)||(k3==k4))#same spin and the same momentum
                continue;
            end
            
            ele=V_k1234[ct];
            if ele==0
                continue;
            end
            
            #C_k4sigma
            state_new,fermi_sign1,state_exist=basis_update_destruction_spinless(N,state1,k4);
            if state_exist==false
                continue;
            end
            
            #C_k3sigma
            state_new,fermi_sign2,state_exist=basis_update_destruction_spinless(N,state_new,k3);
            if state_exist==false
                continue;
            end
            
            #C_k2sigma
            state_new,fermi_sign3,state_exist=basis_update_creation_spinless(N,state_new,k2);
            if state_exist==false
                continue;
            end
            
            #C_k1sigma
            state_new,fermi_sign4,state_exist=basis_update_creation_spinless(N,state_new,k1);
            if state_exist==false
                continue;
            end
            
            #prod(string.(state_new))

            # println([state1,state_new])
            ind_right=c1;
            ind_left=basis_dict[state_new];
            Step_=Step_+1;
            #println([ele,fermi_sign1,fermi_sign2,fermi_sign3,fermi_sign4])
            
            
            if Step_>length(x_ind)#whther the number of matrix elements exceed estimation
                println("extend memory for save Hamiltonian");flush(stdout);
                append!(x_ind, ind_empty_Int);
                append!(y_ind, ind_empty_Int);
                append!(value_ind, ind_empty_Complex);
            end

            x_ind[Step_]=ind_left;
            y_ind[Step_]=ind_right;
            value_ind[Step_]=ele*fermi_sign1*fermi_sign2*fermi_sign3*fermi_sign4;
            
            
        end
    end
    # println(x_ind)
    # println(y_ind)
    x_ind=x_ind[1:Step_];
    y_ind=y_ind[1:Step_];
    value_ind=value_ind[1:Step_];
    #println("Number of matrix elements: "*string(Step_))
    if dim>2000000
        GC.gc(true);#clean memory
    end
    V_upup=sparse(x_ind,y_ind,value_ind);

    return V_upup
end



function interaction_1band_spinup_break_H(H_set::Vector{SparseMatrixCSC},N,manybody_ind, k1234,V_k1234)
    
    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind);
    @assert dim<(2^32); #use Int32 to save index 
    if length(size(V_k1234))==3
        V_k1234=V_k1234[:,spin1,spin2];
    elseif (length(size(V_k1234))==2)&&(size(V_k1234,2)==1)
        V_k1234=V_k1234;
    elseif length(size(V_k1234))==1
        V_k1234=V_k1234;
    else
        error("unknown case")
    end
    #construct disctionary for basis index
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(N,manybody_ind[ci]);
    end
    # basis_dict=Dict(zip(ind_str,manybody_ind))
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))
   
    #ind_empty_Complex=zeros(ComplexF64,dim*N);
    ind_empty_Complex=zeros(ComplexF32,dim*N);
    ind_empty_Int=zeros(Int32,dim*N);
 
    x_ind=deepcopy(ind_empty_Int);
    y_ind=deepcopy(ind_empty_Int);
    value_ind=deepcopy(ind_empty_Complex);
    N_terms=size(k1234,1);
    
    Step_=0;
    for c1=1:dim
        state1=get_bit(N,manybody_ind[c1]);
        for ct=1:N_terms
            (k1,k2,k3,k4)=k1234[ct,:];
            
            
            if ( (k1==k2)||(k3==k4))#same spin and the same momentum
                continue;
            end
            
            ele=V_k1234[ct];
            if ele==0
                continue;
            end
            
            #C_k4sigma
            state_new,fermi_sign1,state_exist=basis_update_destruction_spinless(N,state1,k4);
            if state_exist==false
                continue;
            end
            
            #C_k3sigma
            state_new,fermi_sign2,state_exist=basis_update_destruction_spinless(N,state_new,k3);
            if state_exist==false
                continue;
            end
            
            #C_k2sigma
            state_new,fermi_sign3,state_exist=basis_update_creation_spinless(N,state_new,k2);
            if state_exist==false
                continue;
            end
            
            #C_k1sigma
            state_new,fermi_sign4,state_exist=basis_update_creation_spinless(N,state_new,k1);
            if state_exist==false
                continue;
            end
            
            #prod(string.(state_new))

            # println([state1,state_new])
            ind_right=c1;
            ind_left=basis_dict[state_new];
            Step_=Step_+1;
            #println([ele,fermi_sign1,fermi_sign2,fermi_sign3,fermi_sign4])
            
            


            x_ind[Step_]=ind_left;
            y_ind[Step_]=ind_right;
            value_ind[Step_]=ele*fermi_sign1*fermi_sign2*fermi_sign3*fermi_sign4;
            

            if Step_==length(x_ind)#when the number of matrix elements is full
                x_ind=x_ind[1:Step_];
                y_ind=y_ind[1:Step_];
                value_ind=value_ind[1:Step_];
                Vm_sect=sparse(x_ind,y_ind,value_ind,dim,dim);
                push!(H_set, Vm_sect);
                GC.gc(true);
                H_set_memory=Base.summarysize(H_set)/1024/1024/1024;
                println("Memory cost for H_set: "*string(H_set_memory)*" GB.");flush(stdout);



                println("extend memory for save Hamiltonian");flush(stdout);
                x_ind=deepcopy(ind_empty_Int);
                y_ind=deepcopy(ind_empty_Int);
                value_ind=deepcopy(ind_empty_Complex);
                Step_=0;#A new Hamiltonian section
            end


            
        end
    end
    # println(x_ind)
    # println(y_ind)
    x_ind=x_ind[1:Step_];
    y_ind=y_ind[1:Step_];
    value_ind=value_ind[1:Step_];
    #println("Number of matrix elements: "*string(Step_))
    if dim>2000000
        GC.gc(true);#clean memory
    end

    Vm_sect=sparse(x_ind,y_ind,value_ind,dim,dim);
    push!(H_set, Vm_sect);
    H_set_memory=Base.summarysize(H_set)/1024/1024/1024;
    println("Memory cost for H_set: "*string(H_set_memory)*" GB.");flush(stdout);

    return H_set
end


function interaction_1band_spinup_break_H_new(N,manybody_ind, k1234,V_k1234)
    global H_set
    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind);
    @assert dim<(2^32); #use Int32 to save index 
    if length(size(V_k1234))==3
        V_k1234=V_k1234[:,spin1,spin2];
    elseif (length(size(V_k1234))==2)&&(size(V_k1234,2)==1)
        V_k1234=V_k1234;
    elseif length(size(V_k1234))==1
        V_k1234=V_k1234;
    else
        error("unknown case")
    end
    #construct disctionary for basis index
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(N,manybody_ind[ci]);
    end
    # basis_dict=Dict(zip(ind_str,manybody_ind))
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))
   
    #ind_empty_Complex=zeros(ComplexF64,dim*N);
    ind_empty_Complex=zeros(ComplexF32,dim*N);
    ind_empty_Int=zeros(Int32,dim*N);
 
    x_ind=deepcopy(ind_empty_Int);
    y_ind=deepcopy(ind_empty_Int);
    value_ind=deepcopy(ind_empty_Complex);
    N_terms=size(k1234,1);
    
    Step_=0;
    for c1=1:dim
        state1=get_bit(N,manybody_ind[c1]);
        for ct=1:N_terms
            (k1,k2,k3,k4)=k1234[ct,:];
            
            
            if ( (k1==k2)||(k3==k4))#same spin and the same momentum
                continue;
            end
            
            ele=V_k1234[ct];
            if ele==0
                continue;
            end
            
            #C_k4sigma
            state_new,fermi_sign1,state_exist=basis_update_destruction_spinless(N,state1,k4);
            if state_exist==false
                continue;
            end
            
            #C_k3sigma
            state_new,fermi_sign2,state_exist=basis_update_destruction_spinless(N,state_new,k3);
            if state_exist==false
                continue;
            end
            
            #C_k2sigma
            state_new,fermi_sign3,state_exist=basis_update_creation_spinless(N,state_new,k2);
            if state_exist==false
                continue;
            end
            
            #C_k1sigma
            state_new,fermi_sign4,state_exist=basis_update_creation_spinless(N,state_new,k1);
            if state_exist==false
                continue;
            end
            
            #prod(string.(state_new))

            # println([state1,state_new])
            ind_right=c1;
            ind_left=basis_dict[state_new];
            Step_=Step_+1;
            #println([ele,fermi_sign1,fermi_sign2,fermi_sign3,fermi_sign4])
            

            x_ind[Step_]=ind_left;
            y_ind[Step_]=ind_right;
            value_ind[Step_]=ele*fermi_sign1*fermi_sign2*fermi_sign3*fermi_sign4;
            

            if Step_==length(x_ind)#when the number of matrix elements is full
                x_ind=x_ind[1:Step_];
                y_ind=y_ind[1:Step_];
                value_ind=value_ind[1:Step_];
                Vm_sect=sparse(x_ind,y_ind,value_ind,dim,dim);
                push!(H_set, Vm_sect);
                GC.gc(true);
                H_set_memory=Base.summarysize(H_set)/1024/1024/1024;
                println("Memory cost for H_set: "*string(H_set_memory)*" GB.");flush(stdout);


                println("extend memory for save Hamiltonian");flush(stdout);
                x_ind=deepcopy(ind_empty_Int);
                y_ind=deepcopy(ind_empty_Int);
                value_ind=deepcopy(ind_empty_Complex);
                Step_=0;#A new Hamiltonian section
            end
            
        end
    end
    # println(x_ind)
    # println(y_ind)
    x_ind=x_ind[1:Step_];
    y_ind=y_ind[1:Step_];
    value_ind=value_ind[1:Step_];
    #println("Number of matrix elements: "*string(Step_))
    if dim>2000000
        GC.gc(true);#clean memory
    end

    Vm_sect=sparse(x_ind,y_ind,value_ind,dim,dim);
    push!(H_set, Vm_sect);
    H_set_memory=Base.summarysize(H_set)/1024/1024/1024;
    println("Memory cost for H_set: "*string(H_set_memory)*" GB.");flush(stdout);

end



function interaction_1band_spinup_test(N,manybody_ind, k1234,V_k1234,spin1)
    spin2=spin1;
    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind);
    @assert dim<(2^32); #use Int32 to save index 
    if length(size(V_k1234))==3
        V_k1234=V_k1234[:,spin1,spin2];
    elseif (length(size(V_k1234))==2)&&(size(V_k1234,2)==1)
        V_k1234=V_k1234;
    elseif length(size(V_k1234))==1
        V_k1234=V_k1234;
    else
        error("unknown case")
    end
    #construct disctionary for basis index
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(N,manybody_ind[ci]);
    end
    # basis_dict=Dict(zip(ind_str,manybody_ind))
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))
   
    #ind_empty_Complex=zeros(ComplexF64,dim*N);
    ind_empty_Complex=zeros(ComplexF32,dim*N);
    ind_empty_Int=zeros(Int32,dim*N);
 
    x_ind=deepcopy(ind_empty_Int);
    y_ind=deepcopy(ind_empty_Int);
    value_ind=deepcopy(ind_empty_Complex);
    N_terms=size(k1234,1);


    
    Step_=0;
    for c1=1:dim
        state1=get_bit(N,manybody_ind[c1]);
        

        step_st=1;
        store_ind=ones(Int64,N*N)*(-1);
        store_ele=ones(ComplexF32,N*N)*(-1);
        for ct=1:N_terms
            (k1,k2,k3,k4)=k1234[ct,:];
            
            
            if ( (k1==k2)||(k3==k4))#same spin and the same momentum
                continue;
            end
            
            ele=V_k1234[ct];
            if ele==0
                continue;
            end
            
            #C_k4sigma
            state_new,fermi_sign1,state_exist=basis_update_destruction_spinless(N,state1,k4);
            if state_exist==false
                continue;
            end
            
            #C_k3sigma
            state_new,fermi_sign2,state_exist=basis_update_destruction_spinless(N,state_new,k3);
            if state_exist==false
                continue;
            end
            
            #C_k2sigma
            state_new,fermi_sign3,state_exist=basis_update_creation_spinless(N,state_new,k2);
            if state_exist==false
                continue;
            end
            
            #C_k1sigma
            state_new,fermi_sign4,state_exist=basis_update_creation_spinless(N,state_new,k1);
            if state_exist==false
                continue;
            end
            
            #prod(string.(state_new))

            # println([state1,state_new])
            ind_right=c1;
            ind_left=basis_dict[state_new];
            Step_=Step_+1;
            #println([ele,fermi_sign1,fermi_sign2,fermi_sign3,fermi_sign4])
            
            
            if Step_>length(x_ind)#whther the number of matrix elements exceed estimation
                Vm_tem=sparse(x_ind,y_ind,value_ind);
                println("memory for section of H:")

                Vm_memory=Base.summarysize(Vm_tem)/1024/1024;
                println("Memory cost for Vm: "*string(Vm_memory)*" Mb.");flush(stdout);

                x_ind_memory=Base.summarysize(x_ind)/1024/1024;
                println("Memory cost for x_ind: "*string(x_ind_memory)*" Mb.");flush(stdout);

                value_ind_memory=Base.summarysize(value_ind)/1024/1024;
                println("Memory cost for value_ind: "*string(value_ind_memory)*" Mb.");flush(stdout);


                

                println("extend memory for save Hamiltonian");flush(stdout);
                append!(x_ind, ind_empty_Int);
                append!(y_ind, ind_empty_Int);
                append!(value_ind, ind_empty_Complex);
            end

            x_ind[Step_]=ind_left;
            y_ind[Step_]=ind_right;
            value_ind[Step_]=ele*fermi_sign1*fermi_sign2*fermi_sign3*fermi_sign4;

            store_ind[step_st]=ind_left;
            store_ele[step_st]=ele;
            step_st=step_st+1;
            
        end

        # if c1>10
        #     error(".");
        # end
        # println("lengths:");flush(stdout);
        # println(length(store_ind[1:step_st-1]))
        # println(store_ind[1:step_st-1])
        # println(length(unique(store_ind[1:step_st-1])))
        # println(unique(store_ind[1:step_st-1]))
        # println(store_ele);flush(stdout);
        
    end

    
    # println(x_ind)
    # println(y_ind)
    x_ind=x_ind[1:Step_];
    y_ind=y_ind[1:Step_];
    value_ind=value_ind[1:Step_];
    #println("Number of matrix elements: "*string(Step_))
    if dim>2000000
        GC.gc(true);#clean memory
    end
    V_upup=sparse(x_ind,y_ind,value_ind);

    return V_upup
end





function interaction_1band_spinful(N,manybody_ind, k1234,V_k1234,spin1,spin2)
    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind);
    if length(size(V_k1234))==3
        V_k1234=V_k1234[:,spin1,spin2];
    elseif (length(size(V_k1234))==2)&&(size(V_k1234,2)==1)
        V_k1234=V_k1234;
    elseif length(size(V_k1234))==1
        V_k1234=V_k1234;
    else
        error("unknown case")
    end
    #construct disctionary for basis index
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(2*N,manybody_ind[ci]);
    end
    # basis_dict=Dict(zip(ind_str,manybody_ind))
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))
   
    ind_empty_Complex=zeros(ComplexF64,dim*N);
    ind_empty_Int=zeros(Int128,dim*N);
 
    x_ind=deepcopy(ind_empty_Int);
    y_ind=deepcopy(ind_empty_Int);
    value_ind=deepcopy(ind_empty_Complex);
    N_terms=size(k1234,1);
    
    Step_=0;
    for c1=1:dim
        state1=get_bit(2*N,manybody_ind[c1]);
        for ct=1:N_terms
            (k1,k2,k3,k4)=k1234[ct,:];
            
            
            if (spin1==spin2)&&(( (k1==k2)||(k3==k4)))#same spin and the same momentum
                continue;
            end
            
            ele=V_k1234[ct];
            if ele==0
                continue;
            end
            
            #C_k4sigma
            state_new,fermi_sign1,state_exist=basis_update_destruction(N,state1,spin1,k4);
            if state_exist==false
                continue;
            end
            
            #C_k3sigma
            state_new,fermi_sign2,state_exist=basis_update_destruction(N,state_new,spin2,k3);
            if state_exist==false
                continue;
            end
            
            #C_k2sigma
            state_new,fermi_sign3,state_exist=basis_update_creation(N,state_new,spin2,k2);
            if state_exist==false
                continue;
            end
            
            #C_k1sigma
            state_new,fermi_sign4,state_exist=basis_update_creation(N,state_new,spin1,k1);
            if state_exist==false
                continue;
            end
            
            #prod(string.(state_new))

            # println([state1,state_new])
            ind_right=c1;
            ind_left=basis_dict[state_new];
            Step_=Step_+1;
            #println([ele,fermi_sign1,fermi_sign2,fermi_sign3,fermi_sign4])
            
            
            if Step_>length(x_ind)#whther the number of matrix elements exceed estimation
                println("extend memory for save Hamiltonian");flush(stdout);
                x_ind=vcat(x_ind, ind_empty_Int);
                y_ind=vcat(y_ind, ind_empty_Int);
                value_ind=vcat(value_ind, ind_empty_Complex);
            end

            x_ind[Step_]=ind_left;
            y_ind[Step_]=ind_right;
            value_ind[Step_]=ele*fermi_sign1*fermi_sign2*fermi_sign3*fermi_sign4;
            
            
        end
    end
    # println(x_ind)
    # println(y_ind)
    x_ind=x_ind[1:Step_];
    y_ind=y_ind[1:Step_];
    value_ind=value_ind[1:Step_];

    Vm=sparse(x_ind,y_ind,value_ind);
    return Vm

end




function basis_update_creation(N,state0,spin,k)
    #state0 is a 2xN matrix. The first row corresponds to spin up; the
    #second row corresponds to spin down.
    if state0[(spin-1)*N+k]==Char(0+'0')
        state_new=state0[1:N*(spin-1)+k-1]*'1'*state0[N*(spin-1)+k+1:2*N];
        state_exist=true;
        
        fermi_sign=(-1)^(count(==('1'), state_new[1:N*(spin-1)+k-1]))

    else
        state_new="2";
        state_exist=false;
        fermi_sign=0;
    end
    return state_new,fermi_sign,state_exist

end

function basis_update_destruction(N,state0,spin,k)
    #state0 is a 2xN matrix. The first row corresponds to spin up; the
    #second row corresponds to spin down.

    if state0[N*(spin-1)+k]==Char(1+'0')
        state_new=state0[1:N*(spin-1)+k-1]*'0'*state0[N*(spin-1)+k+1:2*N];
        state_exist=true;
        fermi_sign=(-1)^(count(==('1'), state_new[1:N*(spin-1)+k-1]))
    else
        state_new="2";
        state_exist=false;
        fermi_sign=0;
    end
    return state_new,fermi_sign,state_exist
end



function basis_update_creation_spinless(N,state0,k)

 
    if state0[k]==Char(0+'0')
        state_new=state0[1:k-1]*'1'*state0[k+1:N];
        state_exist=true;
        
        
        #fermi_sign=(-1)^(sum(state_new(1,1:k-1)));
        fermi_sign=(-1)^(count(==('1'), state_new[1:k-1]))

    else
        state_new="2";
        state_exist=false;
        fermi_sign=0;
    end
    return state_new,fermi_sign,state_exist

end

function basis_update_destruction_spinless(N,state0,k)

    if state0[k]==Char(1+'0')
        state_new=state0[1:k-1]*'0'*state0[k+1:N];
        state_exist=true;
        #fermi_sign=(-1)^(sum(state_new(1,1:k-1)));
        fermi_sign=(-1)^(count(==('1'), state_new[1:k-1]))
    else
        state_new="2";
        state_exist=false;
        fermi_sign=0;
    end
    return state_new,fermi_sign,state_exist

end

function kinetic_energy_1band_spinful(N,manybody_ind,E_band)
    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind);
   
    y_ind=zeros(Int128,dim);
    value_ind=zeros(ComplexF64,dim);
    
    
    for c1=1:dim
        Str=get_bit(2*N,manybody_ind[c1]);

        Str_up=Str[1:N];
        Str_up=split(Str_up, "");
        pos=findall(x->x.=="1",Str_up)
        Eup=sum(E_band[pos,1]);

        Str_dn=Str[N+1:2*N];
        Str_dn=split(Str_dn, "");
        pos=findall(x->x.=="1",Str_dn)
        Edn=sum(E_band[pos,2]);
  
        y_ind[c1]=c1;
        value_ind[c1]=Eup+Edn;
    end

    Hk=sparse(y_ind,y_ind,value_ind);
    return Hk
end


function kinetic_energy_1band_spinup(N,manybody_ind,E_band)
    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind);
   
    y_ind=zeros(Int64,dim);
    value_ind=zeros(ComplexF64,dim);
    
    
    for c1=1:dim
        Str=get_bit(N,manybody_ind[c1]);
        Str_=split(Str, "");
        pos=findall(x->x.=="1",Str_)

        E=sum(E_band[pos]);
  
        y_ind[c1]=c1;
        value_ind[c1]=E;
    end

    HK=sparse(y_ind,y_ind,value_ind);
    return HK
end





function reduce_interaction(k1234,V_k1234)
    V_k1234_upup=V_k1234[:];


    siz=size(k1234,1);


    #only keep k1>=k2, k3>=k4 for upup and dndn
    V_k1234_upupnew=deepcopy(V_k1234_upup);
    for c1 in range(1,siz)
        k1234_=k1234[c1,:];
        (k1,k2,k3,k4)=k1234_;
        if (k1<k2)&(k3>=k4)
            k1234_new=[k2,k1,k3,k4];
            Sign=-1;
        elseif (k1<k2)&(k3<k4)
            k1234_new=[k2,k1,k4,k3];
            Sign=1;
        elseif (k1>=k2)&(k3<k4)
            k1234_new=[k1,k2,k4,k3];
            Sign=-1;
        elseif (k1>=k2)&(k3>=k4)
            continue
        else
            error("unknown case")
        end

        V_k1234_upupnew[c1]=0;


        diff=abs.(k1234[:,1].-k1234_new[1]).+abs.(k1234[:,2].-k1234_new[2]).+abs.(k1234[:,3].-k1234_new[3]).+abs.(k1234[:,4].-k1234_new[4]);
        pos=findall(x->x.<1e-10,diff);
        @assert length(pos)==1
        V_k1234_upupnew[pos[1]]=V_k1234_upupnew[pos[1]]+Sign*V_k1234_upup[c1];



    end
    pos_upup=findall(x->abs.(x).>0, V_k1234_upupnew);

    
    k1234_upup=k1234[pos_upup,:];

    V_k1234_upupnew=V_k1234_upupnew[pos_upup,:];


    return k1234_upup,V_k1234_upupnew

end

function reduce_interaction_twoband(k1234,V_k1234)
    #two bands within the same spin and valley
    V_k1234_upup=V_k1234[:];


    siz=size(k1234,1);


    #only keep k1>=k2, k3>=k4 for upup and dndn
    V_k1234_upupnew=deepcopy(V_k1234_upup);
    for c1 in range(1,siz)
        k1234_=k1234[c1,:];
        (k1,k2,k3,k4)=k1234_;
        if (k1<k2)&(k3>=k4)
            k1234_new=[k2,k1,k3,k4];
            Sign=-1;
        elseif (k1<k2)&(k3<k4)
            k1234_new=[k2,k1,k4,k3];
            Sign=1;
        elseif (k1>=k2)&(k3<k4)
            k1234_new=[k1,k2,k4,k3];
            Sign=-1;
        elseif (k1>=k2)&(k3>=k4)
            continue
        else
            error("unknown case")
        end

        V_k1234_upupnew[c1]=0;


        diff=abs.(k1234[:,1].-k1234_new[1]).+abs.(k1234[:,2].-k1234_new[2]).+abs.(k1234[:,3].-k1234_new[3]).+abs.(k1234[:,4].-k1234_new[4]);
        pos=findall(x->x.<1e-10,diff);
        @assert length(pos)==1
        V_k1234_upupnew[pos[1]]=V_k1234_upupnew[pos[1]]+Sign*V_k1234_upup[c1];



    end
    pos_upup=findall(x->abs.(x).>0, V_k1234_upupnew);

    
    k1234_upup=k1234[pos_upup,:];

    V_k1234_upupnew=V_k1234_upupnew[pos_upup,:];


    return k1234_upup,V_k1234_upupnew

end


function get_basis(N,Np,klabel_set,nsec)
    basis_up=spinless_fermion_basis_1d(N,Np);

    occu_set_up,particle_set_up=get_spinless_occu(N,Np,basis_up,length(basis_up));

    k_total_set=klabel_set;
    ksector_size=zeros(Int64,size(k_total_set,1));

    @tensor Kset_up[:]:=occu_set_up[-1,1]*klabel_set[1,-2];
    
    ksector_size,ksector_ind=get_sector_dim_spinup(N,basis_up,Kset_up,k_total_set,nsec);



    #construct disctionary for basis index
    dim=length(ksector_ind);
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(N,ksector_ind[ci]);
    end
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))

    occu_set,particle_set=get_occu(N,Np,ksector_ind);
    return ksector_ind,basis_dict,occu_set,particle_set
end


function get_occu(N,Np,ksector_ind)
    dim=length(ksector_ind);
    occu_set=zeros(Int8,dim,N);
    particle_set=Matrix{Int32}(undef,dim,Np);
    for cc=1:dim
        
        str=get_bit(N,ksector_ind[cc]);
        for dd=1:N
            occu_set[cc,dd]=parse(Int8,str[dd]);
        end
        pos=findall(x->x.==1,occu_set[cc,:]);
        particle_set[cc,1:Np]=pos;

    end
    return occu_set,particle_set
end





function interaction_1band_spinup_single_q(N,manybody_ind, k1234,V1234)
    

    #cdag_{k1} c_{k1+q}  cdag_{k2}  c_{k2-q}
    dim=length(manybody_ind);
    @assert dim<(2^32); #use Int32 to save index 

    #construct disctionary for basis index
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(N,manybody_ind[ci]);
    end
    # basis_dict=Dict(zip(ind_str,manybody_ind))
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))
   
    #ind_empty_Complex=zeros(ComplexF64,dim*N);
    ind_empty_Complex=zeros(ComplexF32,dim*N);
    ind_empty_Int=zeros(Int32,dim*N);
 
    x_ind=deepcopy(ind_empty_Int);
    y_ind=deepcopy(ind_empty_Int);
    value_ind=deepcopy(ind_empty_Complex);
    N_terms=size(k1234,1);
    
    Step_=0;
    for c1=1:dim
        state1=get_bit(N,manybody_ind[c1]);
        for ct=1:N_terms
            (k1,k2,k3,k4)=k1234[ct,:];
            
            
            ele=V1234[ct];
            if ele==0
                continue;
            end
            
        
            
            #C_k4sigma
            state_new,fermi_sign1,state_exist=basis_update_destruction_spinless(N,state1,k4);
            if state_exist==false
                continue;
            end

            #C_k3sigma
            state_new,fermi_sign3,state_exist=basis_update_creation_spinless(N,state_new,k3);
            if state_exist==false
                continue;
            end
            
            #C_k2sigma
            state_new,fermi_sign2,state_exist=basis_update_destruction_spinless(N,state_new,k2);
            if state_exist==false
                continue;
            end
            

            
            #C_k1sigma
            state_new,fermi_sign4,state_exist=basis_update_creation_spinless(N,state_new,k1);
            if state_exist==false
                continue;
            end
            
            #prod(string.(state_new))

            # println([state1,state_new])
            ind_right=c1;
            ind_left=basis_dict[state_new];
            Step_=Step_+1;
            #println([ele,fermi_sign1,fermi_sign2,fermi_sign3,fermi_sign4])
            
            
            if Step_>length(x_ind)#whther the number of matrix elements exceed estimation
                println("extend memory for save Hamiltonian");flush(stdout);
                append!(x_ind, ind_empty_Int);
                append!(y_ind, ind_empty_Int);
                append!(value_ind, ind_empty_Complex);
            end

            x_ind[Step_]=ind_left;
            y_ind[Step_]=ind_right;
            value_ind[Step_]=ele*fermi_sign1*fermi_sign2*fermi_sign3*fermi_sign4;
            
            
        end
    end
    # println(x_ind)
    # println(y_ind)
    x_ind=x_ind[1:Step_];
    y_ind=y_ind[1:Step_];
    value_ind=value_ind[1:Step_];
    println("Number of matrix elements: "*string(Step_))
    if dim>2000000
        GC.gc(true);#clean memory
    end
    V_upup=sparse(x_ind,y_ind,value_ind);

    return V_upup
end






function kinetic_energy_1band_spinup_extract_order(N,manybody_ind1,manybody_ind2,dict1,dict2,k1,k2,v1,v2)
 
    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind2);

    Ob=0;
    
    for c1=1:dim
        state1=get_bit(N,manybody_ind2[c1]);

        state_new,fermi_sign2,state_exist=basis_update_destruction_spinless(N,state1,k2);
        if state_exist==false
            continue;
        end
        

        state_new,fermi_sign1,state_exist=basis_update_creation_spinless(N,state_new,k1);
        if state_exist==false
            continue;
        end
        ind_right=c1;
        ind_left=dict1[state_new];

        Ob=Ob+(v1[ind_left])'*v2[ind_right]*fermi_sign1*fermi_sign2;
    end


    GC.gc(true);
    return Ob
end






function interaction_1band_spinup_qset_distributed(N,manybody_ind, Range, k1234::Matrix,v0s::Matrix)
    v0s=deepcopy(v0s);
    n_keep=size(v0s,1);
    dim=length(manybody_ind);

    Obs=zeros(size(k1234,1),n_keep)*im;;
    
 
  

    #cdag_{k1} c_{k1+q}  cdag_{k2}  c_{k2-q}
    
    @assert dim<(2^32); #use Int32 to save index 

    #construct disctionary for basis index
    ind_str=Vector{String}(undef,dim)
    for ci in range(1,dim)
        ind_str[ci]=get_bit(N,manybody_ind[ci]);
    end
    # basis_dict=Dict(zip(ind_str,manybody_ind))
    basis_dict=Dict(zip(ind_str,Int64.(1:dim)))
   
    #ind_empty_Complex=zeros(ComplexF64,dim*N);
    # ind_empty_Complex=zeros(ComplexF32,dim*N);
    # ind_empty_Int=zeros(Int64,dim*N);
 
    # x_ind=deepcopy(ind_empty_Int);
    # y_ind=deepcopy(ind_empty_Int);
    # value_ind=deepcopy(ind_empty_Complex);
    N_terms=size(k1234,1);
    
    Step_=0;
    for c1=1:dim
        state1=get_bit(N,manybody_ind[c1]);
        for ct=1:N_terms
            (k1,k2,k3,k4)=k1234[ct,:];
            
            
            # ele=V1234[ct,:,:];
            # if norm(ele)==0
            #     continue;
            # end
            
        
            
            #C_k4sigma
            state_new,fermi_sign1,state_exist=basis_update_destruction_spinless(N,state1,k4);
            if state_exist==false
                continue;
            end

            #C_k3sigma
            state_new,fermi_sign3,state_exist=basis_update_creation_spinless(N,state_new,k3);
            if state_exist==false
                continue;
            end
            
            #C_k2sigma
            state_new,fermi_sign2,state_exist=basis_update_destruction_spinless(N,state_new,k2);
            if state_exist==false
                continue;
            end
            

            
            #C_k1sigma
            state_new,fermi_sign4,state_exist=basis_update_creation_spinless(N,state_new,k1);
            if state_exist==false
                continue;
            end
            
            #prod(string.(state_new))

            # println([state1,state_new])
            ind_right=c1;
            ind_left=basis_dict[state_new];
            Step_=Step_+1;
            #println([ele,fermi_sign1,fermi_sign2,fermi_sign3,fermi_sign4])
            
            

            # for cs=1:n_keep
            #     Obs[cs,:,:]=Obs[cs,:,:]+(ele*fermi_sign1*fermi_sign2*fermi_sign3*fermi_sign4)*(v0s[cs][ind_right])*(v0s[cs][ind_left])';
            # end

            Obs[ct,:]=Obs[ct,:]+(v0s[:,ind_right]).*(conj.(v0s[:,ind_left]))*(fermi_sign1*fermi_sign2*fermi_sign3*fermi_sign4);
 
            
        end
    end
    # println(x_ind)
    # println(y_ind)
    # x_ind=x_ind[1:Step_];
    # y_ind=y_ind[1:Step_];
    # value_ind=value_ind[1:Step_];
    # println("Number of matrix elements: "*string(Step_))
    # if dim>2000000
    manybody_ind=[];
    ind_str=[];
    basis_dict=[];
    v0s=[];
    k1234=[];
    GC.gc(true);#clean memory
    # end
    # V_upup=sparse(x_ind,y_ind,value_ind);

    return (Range,Obs)
end



function kinetic_energy_1band_spinup_Eband(N,manybody_ind,cq,v0s,E_band1,E_band2)


    #Cdag_k1sigma, Cdag_k2sigma', C_k3sigma',  C_k4sigma
    dim=length(manybody_ind);
    n_keep=length(v0s);


    Obs1=zeros(n_keep)*im;
    for c1=1:dim
        Str=get_bit(N,manybody_ind[c1]);
        Str_=split(Str, "");
        pos=findall(x->x.=="1",Str_)

        E=sum(E_band1[pos]);
        
        for cs=1:n_keep
            Obs1[cs]=Obs1[cs]+E*v0s[cs][c1]*(v0s[cs][c1]');
        end
    end

    Obs2=zeros(n_keep)*im;
    for c1=1:dim
        Str=get_bit(N,manybody_ind[c1]);
        Str_=split(Str, "");
        pos=findall(x->x.=="1",Str_)

        E=sum(E_band2[pos]);
        
        for cs=1:n_keep
            Obs2[cs]=Obs2[cs]+E*v0s[cs][c1]*(v0s[cs][c1]');
        end
    end

    manybody_ind=[];
    GC.gc(true);
    return (cq,Obs1,Obs2)
end

