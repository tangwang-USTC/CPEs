
"""
  For isotropic plasma:

  The normalzied kinetic moments for shell-less quasi-equilibrium state plasma 
    when the velocity space exhibits spherical symmetry without shell structure.
    When the zeroth-order amplitude of the normalized distribution function, 
    `f̂₀(v̂)`, is approximated by the MMM, 
    the `jᵗʰ`-order normalzied kinetic moment can be expressed as:
  
        `𝓜ⱼ(f̂₀) = 4π * ∫₀^∞(v̂ʲ⁺² * f̂₀) dv̂
                 = CMj0 * ∑ᵣ₌₁ᴺᴷ n̂ₐᵣ(v̂ₐₜₕᵣ)ʲ, j ≥ -2`

  where `f̂ₗ(v̂) = vₜₕ³/nₐ * f(v)` and the coefficient

        `CMj0 = 2 / √π * Γ((j+3)/2)`.
  
  When `j` is even, 

    `CMj0 = (j+1)!! / 2^(j/2),  j ∈ -2:2:N⁺`.

  When `j` is odd,

    'CMj0 = 2 / √π * ((j+1)/2)!,  j ∈ -1:2:N⁺'.

  If `is_renorm == true`,
    `𝓜ⱼ(f̂₀) /= CMj0`
  end
  
  See Ref. of Wang (2024) titled as "General relaxation model for a homogeneous plasma with spherically symmetric velocity space".

  Inputs:
    jvec: j ≥ -2, j ∈ ℕ⁺-3

  Outputs
    Mhst = MhsMMM!(Mhst,jvec,nai,vthi,nMod,ns;is_renorm=is_renorm)
    Mhst = MhsMMM!(Mhst,jvec,ns;is_renorm=is_renorm)

"""

# 2.5D, [nMod,njMs,ns]
function MhsMMM!(Mhst::AbstractArray{T},jvec::Vector{Int},
    nai::Vector{AbstractVector{T}},vthi::Vector{AbstractVector{T}},
    nMod::Vector{Int},ns::Int64;is_renorm::Bool=true) where{T}
    
    for isp in 1:ns
        if nMod == 1
            MhsMMM!(Mhst[:,isp],jvec;is_renorm=is_renorm)
        else
            vec = 1:nMod[isp]
            MhsMMM!(Mhst[:,isp],jvec,nai[isp][vec],vthi[isp][vec];is_renorm=is_renorm)
        end
    end
end

# 1.5D, [nMod,njMs]
function MhsMMM!(Mhst::AbstractVector{T},jvec::Vector{Int},
    nai::AbstractVector{T},vthi::AbstractVector{T};is_renorm::Bool=true) where{T}
    
    if prod(isone.(vthi))
        k = 0
        for j in jvec
            k += 1
            Mhst[k] = MhsMMM(j,T;is_renorm=is_renorm)
        end
        Mhst[:] *= sum_kbn(nai)
    else
        k = 0
        for j in jvec
            k += 1
            Mhst[k] = MhsMMM(j,T;is_renorm=is_renorm) * sum_kbn(nai .* vthi.^j)
        end
    end
end

# 0.5D, [nMod]
function MhsMMM(j::Int64,nai::AbstractVector{T},vthi::AbstractVector{T};is_renorm::Bool=true) where{T}
    
    if prod(isone.(vthi))
        return MhsMMM(j,T;is_renorm=is_renorm) * sum_kbn(nai)
    else
        return MhsMMM(j,T;is_renorm=is_renorm) * sum_kbn(nai .* vthi.^j)
    end
end

"""
"""

# nMod = 1 -> nai = 1, vthi = 1

# 2D, [njMs,ns]
function MhsMMM!(Mhst::AbstractArray{T},jvec::Vector{Int},ns::Int64;is_renorm::Bool=true) where{T}

    for isp in 1:ns
        MhsMMM!(Mhst[:,isp],jvec;is_renorm=is_renorm)
    end
end

# 1D, [njMs]
function MhsMMM!(Mhst::AbstractVector{T},jvec::Vector{Int};is_renorm::Bool=true) where{T}

    if is_renorm
        Mhst[:] .= 1.0
    else
        k = 0
        for j in jvec
            Mhst[k+1] = MhsMMM(j,T;is_renorm=is_renorm)
        end
    end
end

# 0D, []
function MhsMMM(j::Int,T;is_renorm::Bool=true)

    if is_renorm
        return 1.0 |> T 
    else
        return CMjL(T(j))
    end
end
