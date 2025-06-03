

## Final Function!

using LinearAlgebra, SparseArrays

function exactsln(f,N,T,alph)

    # import Pkg; Pkg.add("FFTW")
    # import Pkg; Pkg.add("Plots")
    # import Pkg; Pkg.add("MatrixEquations")
    # using Pkg
    # Pkg.add("SpecialFunctions")
    # using Random
    # using SpecialFunctions

    #  This is the julia file with the functions! 

    include("functions.jl")



    ## Setting up the points ##

    # points in t direction

    t = pi* range(-1,1, length=N+1)[1:end-1];
    t = reshape(t,:,1)


    # points in r and z direction

    r = [cos((pi*k)/(N-1)) for k in (N-1):-1:0];
    z = [cos((pi*k)/(N-1)) for k in (N-1):-1:0];

    ## Setting up the grid ##

    rr = [ri for ri in r, zi in z, ti in t]
    zz = [zi for ri in r, zi in z, ti in t]
    tt = [ti for ri in r, zi in z, ti in t]

    rr = dropdims(rr; dims=4)
    tt = dropdims(tt; dims=4)
    zz = dropdims(zz; dims=4)



    ## Matrices ##

    # C01: Convert Chebyshev-T (T) coefficients to Chebyshev-U (U)
    C01 = spdiagm(0 => 0.5*ones(N), 2 => -0.5*ones(N-2))
    C01[1,1] = 1.0
    
    # C12: Convert Chebyshev-U to C^2, the continuous function space
    K_1 = 1.0 ./ (1:N)
    K_2 = 1.0 ./ (3:N)
    C12 = spdiagm(0 => K_1, 2 => -K_2)
    
    # C02: Convert Chebyshev-T to C^2 (It's C01 * C12)
    C02 = C01 * C12

    # D1: Differentiate Chebyshev-T and convert to Chebyshev-U
    D1 = spdiagm(1 => 1:N-1)
    
    # D2: Differentiate Chebyshev-T twice and convert to C²
    D2 = spdiagm(2 => 2 .* (2:N-1))
    
    # D0: Differentiate Chebyshev-T once
    D0 = C01 \ D1

    # R0: Multiply Cheb-T by x
    R0 = spdiagm(-1 => 0.5*ones(N-1), 1 => 0.5*ones(N-1))
    R0[2,1] = 1.0
    
    # R2: Multiply in C² space
    denom = 4:2:2*N+2
    K_sub  = ((1:N) ./ denom)[1:N-1]
    K_super = ((3:N+2) ./ denom[1:N])[2:N]
    R2 = spdiagm(-1 => K_sub, 1 => K_super)
    
    # R1: Multiply in Cheb-U basis
    R1 = C12 \ R2 * C12

    # A: Laplacian operator
    A = R2 * R2 * D2 + R2 * C12 * D1
    
    # coF: Laplacian multiplier
    coF = R2 * R2 * C02



    ## Evaluate Initial Condition at points at grid points ##

    CFS = zeros(size(rr))

    # Evaluate initial condition on the full (r, z, θ) grid
    for j in 1:N
        CFS[:, j, :] .= finit(rr[:,j,:], zz[:,j,:], tt[:,j,:], 0)
    end



    ## First 4 time steps

    CFS_func = V2C_cyl(CFS, "rzt")

    CFS_4D = complex(zeros(Float64, size(CFS_func)..., T))  # Create a new 4D array full of zeros of size (5,5,5,5)
    CFS_4D[:,:,:,1] = CFS_func # Copy over CFS_func into the first 4th slice

    Hist = Vector{Array{Float64,3}}(undef, 3)  # Equivalent to cell(3,1)

    # Plug in exact solution for the first four time steps
    for j in 1:3
        Hist[j] = Array{Float64,3}(undef, size(rr))  # Allocate space for each time slice
        for k in 1:size(rr, 2)  # Assuming k indexes the 2nd dimension
            Hist[j][:,k,:] = finit(rr[:,k,:], zz[:,k,:], tt[:,k,:], j * alph)
        end
        CFS_4D[:,:,:,j+1] = V2C_cyl(Hist[j], "rzt")
    end



    ## Full time step ##

    q = (12/25)*alph

    for d = 5:T
        RHS = (48/25)*CFS_4D[:,:,:,d-1]-(36/25)*CFS_4D[:,:,:,d-2]+(16/25)*CFS_4D[:,:,:,d-3]-(3/25)*CFS_4D[:,:,:,d-4]
        CFS_4D[:,:,:,d] = Heat3DSolver(RHS,q,N,C02,A,D2,coF)
    end


    
    ## Returning Coefficients back to Values

    CFS_4D[:,:,:,T] = C2V_cyl(CFS_4D[:,:,:,T],"rzt");

    return 

    CFS_4D

end 