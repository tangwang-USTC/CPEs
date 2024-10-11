
"""
  The normalzied kinetic moments for quasi-equilibrium state plasma 
    when the velocity space exhibits spherical symmetry.
    When the zeroth-order amplitude of the normalized distribution function, 
    `f̂₀(v̂)`, is approximated by the KMM0, 
    the `jᵗʰ`-order normalzied kinetic moment can be expressed as:
  
        `𝓜ⱼ(f̂₀) = 4π * ∫₀^∞(v̂ʲ⁺² * f̂₀) dv̂
                 = CMj0 * ∑ᵣ₌₁ᴺᴷ{n̂ₐᵣ(v̂ₐₜₕᵣ)ʲ * ₁F₁[-j/2,3/2,-(ûₐᵣ/v̂ₐₜₕᵣ)²]}, j ≥ -2`

  where `f̂ₗ(v̂) = vₜₕ³/nₐ * f(v)`,
  `₁F₁(a,b,z)` represents the Kummer confluent hypergeometric1F1 function and the coefficient

        `CMj0 = 2 / √π * Γ((j+3)/2)`.
  
  When `j` is even, the normalzied kinetic moment can be expressed as:

    `𝓜ⱼ(f̂₀) = CMj0 * ∑ᵣ₌₁ᴺᴷ{n̂ₐᵣ(v̂ₐₜₕᵣ)ʲ * [1 + ∑ₖ₌₁^(j/2) Cj0k * (ûₐᵣ/v̂ₐₜₕᵣ)²ᵇ]} , j ∈ {(2jₚ - 2) | jₚ ∈ [0,N⁺]}`,

  where

    `CMj0 = (j+1)!! / 2^(j/2),  j ∈ -2:2:N⁺`.

  If `is_renorm == true`,
    `𝓜ⱼ(f̂₀) /= CMj0  , j ∈ {(2jₚ - 2) | jₚ ∈ [0,N⁺]}`
  end

  When `j` is odd,

    'CMj0 = 2 / √π * ((j+1)/2)!,  j ∈ -1:2:N⁺'.
  
  See Ref. of Wang (2024) titled as "General relaxation model for a homogeneous plasma with spherically symmetric velocity space".

  Inputs:
    jvec:

  Outputs
    Mhst = MhsKMM0!(Mhst,jvec,nai,uai,vthi,nMod,ns;is_renorm=is_renorm)
    Mhst = MhsKMM0!(Mhst,jvec,uai,ns;is_renorm=is_renorm)

"""

# 2.5D, [nMod,njMs,ns]
function MhsKMM0!(Mhst::AbstractArray{T},jvec::Vector{Int},
    nai::Vector{AbstractVector{T}},uai::Vector{AbstractVector{T}},vthi::Vector{AbstractVector{T}},
    nMod::Vector{Int},ns::Int64;is_renorm::Bool=true) where{T}
    
    for isp in 1:ns
        if nMod == 1
            MhsKMM0!(Mhst[:,isp],jvec,uai[isp][vec];is_renorm=is_renorm) 
        else
            vec = 1:nMod[isp]
            MhsKMM0!(Mhst[:,isp],jvec,nai[isp][vec],uai[isp][vec],vthi[isp][vec],nMod[isp];is_renorm=is_renorm) 
        end
    end
    return Mhst
end

# 1.5D, [nMod,njMs]
function MhsKMM0!(Mhst::AbstractVector{T},jvec::Vector{Int},
    nai::AbstractVector{T},uai::AbstractVector{T},vthi::AbstractVector{T},
    nMod::Int;is_renorm::Bool=true) where{T}
    
    if sum(abs.(uai)) ≤ eps(T)
        MhsMMM!(Mhst,jvec,nai,vthi;is_renorm=is_renorm) 
    else
        k = 0
        for j in jvec
            k += 1
            Mhst[k] = MhsKMM0(j,nai,uai,vthi,nMod;is_renorm=is_renorm) 
        end
    end
    return Mhst
end

# 0.5D, [nMod]
function MhsKMM0(j::Int64,nai::AbstractVector{T},uai::AbstractVector{T},vthi::AbstractVector{T},
    nMod::Int;is_renorm::Bool=true) where{T}
    
    if sum(abs.(uai)) ≤ eps(T)
        return MhsMMM(j,nai,vthi;is_renorm=is_renorm) 
    else
        Mh = 0.0
        for s = 1:nMod
            Mh += MhsKMM0(j,nai[s],uai[s],vthi[s];is_renorm=true) 
        end
        if is_renorm
            return Mh
        else
            return CMjL(j) * Mh
        end
    end
end

# 0.5D, []
function MhsKMM0(j::Int64,nai::T,uai::T,vthi::T;is_renorm::Bool=true) where{T}
    
    if sum(abs.(uai)) ≤ eps(T)
        if isone(vthi)
            return MhsMMM(j; is_renorm) * nai
        else
            return MhsMMM(j; is_renorm) * (nai * vthi^j)
        end
    else
        if iseven(j)
            if isone(vthi)
                a = 1.0
                for k in 1:j/2
                    a += CjLk(j,k) * uai^(2k)
                end
                if is_renorm
                    return a * nai
                else
                    return CMjL(j) * a * nai
                end
            else
                a = 1.0
                for k in 1:j/2
                    a += CjLk(j,k) * (uai/vthi)^(2k)
                end
                if is_renorm
                    return a * nai * vthi^j
                else
                    return CMjL(j) * a * nai * vthi^j
                end
            end
        else
            sdfgbhnm
        end
    end
end

"""
"""

# nMod = 1 -> nai = 1, vthi = 1

# 2D, [njMs,ns]
function MhsKMM0!(Mhst::AbstractArray{T},jvec::Vector{Int},uai::T,ns::Int64;is_renorm::Bool=true) where{T}

    for isp in 1:ns
        MhsKMM0!(Mhst[:,isp],jvec,uai;is_renorm=is_renorm) 
    end
    return Mhst
end

# 1D, [njMs]
function MhsKMM0!(Mhst::AbstractVector{T},jvec::Vector{Int},uai::T;is_renorm::Bool=true) where{T}

    k = 0
    for j in jvec
        Mhst[k+1] = MhsKMM0(j,uai;is_renorm=is_renorm)
    end
    return Mhst
end

# 0D, []   
function MhsKMM0(j::Int,uai::T;is_renorm::Bool=true) where{T}
    
    if abs(uai) ≤ eps(T)
        return MhsMMM(j)
    else
        if iseven(j)
            a = 1.0
            for k in 1:j/2
                a += CjLk(j,k) * uai^(2k)
            end
            if is_renorm
                return a
            else
                return CMjL(j) * a
            end
        else
            sdfgbhnm
        end
    end
end
