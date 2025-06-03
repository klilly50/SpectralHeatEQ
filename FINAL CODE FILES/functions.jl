# import Pkg; Pkg.add("FFTW")
# import Pkg; Pkg.add("Plots")
# import Pkg; Pkg.add("MatrixEquations")
# using Pkg
# Pkg.add("SpecialFunctions")
# using Random
# using SpecialFunctions

## replicate chebtech2.vals2coeffs 
using FFTW
using LinearAlgebra

function chebv2c(vals)
    N = size(vals, 1)
    coeffs = zeros(size(vals))

    for col = 1:size(vals, 2)
        v = vals[:, col] # single column vector
        v_backwards = v[end-1:-1:2] # backwards
        v_ext = vcat(v, v_backwards) # mirrored but not the first at the end, so like 1 2 3 2
        c_ext = real(fft(v_ext)) # fft but only take real part
        c = c_ext[1:N] / (N-1) # normalizes the first half, N = first half size
        c[1] = c[1] / 2
        c[end] = c[end] / 2
        coeffs[:, col] = c
    end

    return coeffs

end

## replicate trigtech.vals2coeffs function

function trigv2c(vals)
    N = size(vals, 1)

    coeffs = (1/N)*fftshift(fft(vals, 1), 1)
    # display("coeffs below")
    # display(size(coeffs))
    # display(coeffs)

    if ( mod(N, 2) != 0 ) # N is odd
        even_odd_fix = (-1).^(-(N-1)/2:(N-1)/2)
    else # N is even
        even_odd_fix = (-1).^((-N/2):(N/2-1))
    end

    # display("evenoddfix size")
    # display(size(even_odd_fix))
    # display(even_odd_fix)

    coeffs = coeffs .* even_odd_fix

    return coeffs
end

## V2C_cyl

function V2C_cyl(vals, a)
    dim = length(a)
    SIZE = size(vals)

    a = collect(lowercase(a))

    for kk in 1:dim
        # first, "vectorize" to catch all columns at once.
        vals = reshape(vals, SIZE[kk], div(length(vals), SIZE[kk]))

        # apply the appropriate transform
        if a[kk] == 't'
            vals = trigv2c(vals)
        elseif a[kk] == '0'
            # do nothing
        else
            vals = chebv2c(vals)
        end

        # Reshape and permute back
        if dim != 1
            vals = reshape(vals, SIZE)
            vals = permutedims(vals, (2, 3, 1))
        end

        if a[kk] == 'r' || a[kk] == 'z'
            for j in 2:2:SIZE[3]
                vals[:, :, j] .= -vals[:, :, j]
            end
        end
    end

    return vals

end

## Heat 3D Solver

using MatrixEquations

function Heat3DSolver(cfs,Q,N,C02,A,D2,coF)

    # % Solves the screened Poisson equation, [I-Q*Del^2]u=CFS in coeff. space
    # % Homogeneous Dirichlet boundary conditions.
    
    # % Discretization size: N
    
    Fapprox = zeros(ComplexF64, N,N,N)
    
    # bc = [ (-1).^(0:N-1) > 0 ;  (-1).^(0:N-1) < 0 ]
    # bc = [ ((-1).^(0:N-1) .> 0) ; ((-1).^(0:N-1) .< 0) ]
    bc = Float64.(vcat(((-1).^(0:N-1) .> 0)', ((-1).^(0:N-1) .< 0)'))
    
    for k = -(N-1)//2 : (N-1)//2
        k_idx = Int(k + (N+1)//2)
        
        AA = coF - Q*(A - k^2*C02)
        BB = C02
        CC = -Q*coF
        DD = D2
        RHS = coF*cfs[:,:,k_idx]*C02'
        #%RHS = RHS(1:N-2, 1:N-2)
        
        (AA, RHS) = zeroDOF(AA, BB, RHS, bc, zeros(2,N))
        (BB, RHS) = zeroDOF(BB, AA, RHS', bc, zeros(2,N))
        RHS = RHS'
        (CC, RHS) = zeroDOF(CC, DD, RHS, bc, zeros(2,N))
        (DD, RHS) = zeroDOF(DD, CC, RHS', bc, zeros(2,N))
        RHS = RHS'
        
        AA = AA[1:end-2,3:end]
        BB = BB[1:end-2,3:end]
        CC = CC[1:end-2,3:end]
        DD = DD[1:end-2,3:end]
        RHS = RHS[1:end-2,1:end-2]
        
        X22 = gsylv(AA, BB, CC, DD, RHS )
        X12 = bc[1:2,1:2] \ (-bc[1:2,3:end]*X22)
        X21 = ( bc[1:2,1:2]' \ (-bc[1:2,3:end]*X22') )'
        X11 = bc[1:2,1:2] \ (-bc[1:2,3:end]*X21)
        X = [ X11 X12 ; X21 X22 ]
        Fapprox[:, :, k_idx] = X
    
    end

    return Fapprox
    
end
 
## zeroDOF

function zeroDOF(C1, C2, E, B, G)
    # %ZERODOF   Eliminate so degrees of freedom in the matrix equation can be
    # %removed.
    
    for ii = 1:size(B, 1) #% For each boundary condition, zero a column.
        for kk = 1:size(C1, 1)
            if ( abs(C1[kk,ii]) > 10*eps(Float64) )
                c = C1[kk, ii] #% Constant required to zero entry out.
                C1[kk,:] = C1[kk,:] - c*B[ii,:]
                E[kk,:] = E[kk,:] - ComplexF64(c) .* (reshape(G[ii,:], 1, :) * C2')[1,:]
            end
        end
    end

    return (C1, E)
    
end

##  chebtech2.coeffs2vals function

function chebc2v(coeffs)
    n = size(coeffs)[1]

    # Trivial case (constant or empty)
    if n <= 1
        return coeffs
    end

    # Check for symmetry
    # check if all even-indexed (AKA odd Chebyshev coefficients T1, T3,...) are zero
    isEven = [maximum(abs.(coeffs[2:2:end, :])) == 0]
    # check if all odd-indexed (AKA even Chebyshev coefficients T2, T4,...) are zero
    isOdd  = [maximum(abs.(coeffs[1:2:end, :])) == 0]
    # isEven and isOdd are arrays of booleans telling whether the column in coeffs is even or odd symmetric

    # Scale them by 1/2
    coeffs = copy(coeffs)  # avoid mutating input
    coeffs[2:n-1, :] .= coeffs[2:n-1, :] ./ 2

    # Mirror the coefficients (to fake a DCT using an FFT):
    tmp = vcat(coeffs, coeffs[n-1:-1:2, :])

    # Do FFT dependent on realness
    if isreal(coeffs)
        values = real.(fft(tmp))
    elseif isreal(1im .* coeffs)
        values = 1im .* real.(fft(imag.(tmp)))
    else
        values = fft(tmp)
    end

    # flip and truncate
    values = values[n:-1:1, :]
    # enforce symmetry
    values[:, isEven] = (values[:, isEven] .+ reverse(values[:, isEven])) ./ 2
    values[:, isOdd] = (values[:, isOdd] .- reverse(values[:, isOdd])) ./ 2
  
    return values
end

## trigtech.coeffs2vals function

function trigc2v(coeffs)
    N = size(coeffs, 1)
    values = similar(coeffs)

    # Fix sign alternation: invert the earlier transformation
    if ( mod(N, 2) != 0 )
        even_odd_fix = (-1.0).^(-(N-1)/2:(N-1)/2)
    else
        even_odd_fix = (-1.0).^((-N/2):(N/2-1))
    end
    coeffs = coeffs .* even_odd_fix  # bsxfun

    # Inverse FFT and shift
    # next line in MATLAB: values = ifft(ifftshift(coeffs, 1), [], 1) * n;
    values = ifft(ifftshift(coeffs, 1), 1) * N

    # translate to Julia

    # Enforce Hermitian/Skew symmetry for real-valued outputs
    # vals_full = vcat(values, values[1:1, :])
    # isHerm = vec(all(abs.(vals_full .- conj(reverse(vals_full, dims=1))) .< 1e-14, dims=1))
    # isSkew = vec(all(abs.(vals_full .+ conj(reverse(vals_full, dims=1))) .< 1e-14, dims=1))

    # values[:, isHerm] .= real.(values[:, isHerm])
    # values[:, isSkew] .= imag.(values[:, isSkew]) .* 1im

    return values
end

## C2V_cyl Function

function C2V_cyl(vals, a)
    dim = length(a)
    SIZE = size(vals)

    a = collect(lowercase(a))

    for kk in 1:dim
        # first, "vectorize" to catch all columns at once.
        vals = reshape(vals, SIZE[kk], div(length(vals), SIZE[kk]))

        # apply the appropriate transform
        if a[kk] == 't'
            for ii = 1:prod(SIZE)/SIZE[kk]
                ii = Int(ii)
                vals[:, ii] = trigc2v(vals[:,ii])
            end
        elseif a[kk] == '0'
            # do nothing
        else
            for ii = 1:prod(SIZE)/SIZE[kk]
                ii = Int(ii)
                vals[:, ii] = chebc2v(vals[:,ii])
            end
        end

        # Reshape and permute back
        if dim != 1
            vals = reshape(vals, SIZE)
            vals = permutedims(vals, (2, 3, 1))
        end 
    end

    vals = real(vals)
    return vals

end

## ngrid replication

function ndgrid(arrays::AbstractVector...)
    sizes = map(length, arrays)
    grids = ntuple(d -> begin
        shape = ones(Int, length(arrays))
        shape[d] = sizes[d]
        reshape(arrays[d], shape...) |> a -> repeat(a, sizes ./ shape)
    end, length(arrays))
    return grids
end

## func2grid Function

function func2grid(f, N, ord, coeff)

    # recreate trigpts and chebpts 

    t = pi* range(-1,1, length=N+1)[1:end-1];
    t = reshape(t,:,1)

    r = [cos((pi*k)/(N-1)) for k in (N-1):-1:0];

    # using julia version of linspace (range) and turning it into a column vector

    eq = range(-1, 1, N)
    collect(eq)

    # set up grid replicate cell() in MATLAB using Julia Array{}

    grid = Array{Any}(undef, length(ord), 1)

    # only thing that changed in loops was indexing () -> []

    for i in 1:length(ord)
        if ( ord[i] == 't' || ord[i] == 'T' )
            grid[i] = t;
        elseif ( ord[i] == 'e' || ord[i] == 'E' )
            grid[i] = eq;
        else
            grid[i] = r;
        end
    end

    GRID = ndgrid(grid...)
    X = f(GRID...)

    return

    X

end

## Functions for Printing

using Printf

function pretty_print_3Darray_matlab_style(A::Array; digits=4)
    println()

    iscomplexarray = eltype(A) <: Complex

    for k in 1:size(A, 3)
        slice = A[:, :, k]
        abs_slice = iscomplexarray ? abs.(slice) : abs.(real.(slice))
        maxval = maximum(abs_slice)
        
        println("(:,:, $k) =\n")

        if maxval == 0.0
            scale = 0
            scaled = slice
        else
            scale = floor(Int, log10(maxval))
            scaled = slice / (10.0^scale)
        end

        use_scale = maxval < 1e-3 || maxval > 1e4

        if use_scale
            println(@sprintf("   1.0e%+d *\n", scale))
        end

        for i in 1:size(slice, 1)
            for j in 1:size(slice, 2)
                val = use_scale ? scaled[i, j] : slice[i, j]
                if iscomplexarray
                    re = @sprintf("%7.*f", digits, real(val))
                    im = @sprintf("%7.*f", digits, imag(val))
                    sign_str = startswith(im, "-") ? "" : "+"
                    print(" $re $sign_str $im" * "i  ")
                else
                    print(@sprintf("%8.*f  ", digits, val))
                end
            end
            println()
        end
        println()
    end
end

function pretty_print_4Darray_matlab_style(A::Array; digits=4)
    println()
    
    if ndims(A) != 4
        error("This function only works with 4D arrays.")
    end

    iscomplexarray = eltype(A) <: Complex

    for l in 1:size(A, 4)
        for k in 1:size(A, 3)
            slice = A[:, :, k, l]
            abs_slice = iscomplexarray ? abs.(slice) : abs.(real.(slice))
            maxval = maximum(abs_slice)
            
            println("(:,:,$k,$l) =\n")
            
            if maxval == 0.0
                scale = 0
                scaled = slice
            else
                scale = floor(Int, log10(maxval))
                scaled = slice / (10.0^scale)
            end
            
            use_scale = maxval < 1e-3 || maxval > 1e4
            
            if use_scale
                println(@sprintf("   1.0e%+d *\n", scale))
            end
            
            for i in 1:size(slice, 1)
                for j in 1:size(slice, 2)
                    val = use_scale ? scaled[i, j] : slice[i, j]
                    if iscomplexarray
                        re = @sprintf("%7.*f", digits, real(val))
                        im = @sprintf("%7.*f", digits, imag(val))
                        sign_str = startswith(im, "-") ? "" : "+"
                        print(" $re $sign_str $im" * "i  ")
                    else
                        print(@sprintf("%8.*f  ", digits, val))
                    end
                end
                println()
            end
            println()
        end
    end
end

