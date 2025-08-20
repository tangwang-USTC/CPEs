
"""
  For isotropic plasma:

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

    `𝓜ⱼ(f̂₀) = CMj0 * ∑ᵣ₌₁ᴺᴷ{n̂ₐᵣ(v̂ₐₜₕᵣ)ʲ(ûₐᵣ)ʲ⁻ˡ * [1 + ∑ₖ₌₁^(j/2) Cj0k * (ûₐᵣ/v̂ₐₜₕᵣ)²ᵇ]} , j ∈ {(2jₚ - 2) | jₚ ∈ [0,N⁺]}`,

  where

    `CMj0 = (j+1)!! / 2^(j/2),  j ∈ -2:2:N⁺`.

  If `is_renorm == true`,
    `𝓜ⱼ(f̂₀) /= CMj0  , j ∈ {(2jₚ - 2) | jₚ ∈ [0,N⁺]}`
  end

  When `j` is odd,

    'CMj0 = 2 / √π * ((j+1)/2)!,  j ∈ -1:2:N⁺'.
  
  See Ref. of Wang (2024) titled as "General relaxation model for a homogeneous plasma with spherically symmetric velocity space".

"""

"""
  Inputs:
    jvec:
    mathtype: [:Exact, :Taylor0, :Taylor1, :TaylorInf]

  Outputs
    Mhst = MhsKMM0!(Mhst,jvec,nai,uai,vthi,nMod,ns;is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype)
    Mhst = MhsKMM0!(Mhst,jvec,uai,ns;is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype)

"""

# 2.5D, [nMod,njMs,ns]
function MhsKMM0!(Mhst::AbstractArray{T},jvec::Vector{Int},
    nai::Vector{AbstractVector{T}},uai::Vector{AbstractVector{T}},vthi::Vector{AbstractVector{T}},
    nMod::Vector{Int},ns::Int64;is_renorm::Bool=true,rtol_OrjL::T=1e-10,mathtype::Symbol=:Exact) where{T}
    
    for isp in 1:ns 
        if nMod == 1
            MhsKMM0!(Mhst[:,isp],jvec,uai[isp][vec];is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype) 
        else
            vec = 1:nMod[isp] 
            MhsKMM0!(Mhst[:,isp],jvec,nai[isp][vec],uai[isp][vec],vthi[isp][vec],nMod[isp];
                    is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype) 
        end
    end
end

"""
  Inputs:
    jvec:
    mathtype: [:Exact, :Taylor0, :Taylor1, :TaylorInf]

  Outputs
    Mhst = MhsKMM0!(Mhst,jvec,nai,uai,vthi,nMod;
                    is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype)
    Mhst = MhsKMM0!(Mhst,jvec,uai;is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype)

"""

# 1.5D, [nMod,njMs]
function MhsKMM0!(Mhst::AbstractVector{T},jvec::Vector{Int},
    nai::AbstractVector{T},uai::AbstractVector{T},vthi::AbstractVector{T},
    nMod::Int;is_renorm::Bool=true,rtol_OrjL::T=1e-10,mathtype::Symbol=:Exact) where{T}
    
    if sum_kbn(abs.(uai)) ≤ eps(T)
        MhsMMM!(Mhst,jvec,nai,vthi;is_renorm=is_renorm) 
    else
        k = 0 
        for j in jvec 
            k += 1
            Mhst[k] = MhsKMM0(j,nai,uai,vthi,nMod;
                             is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype) 
        end
    end
end


"""
  Inputs:
    jvec:
    mathtype: [:Exact, :Taylor0, :Taylor1, :TaylorInf]

  Outputs
    Mhst = MhsKMM0(j,nai,uai,vthi,nMod;is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype)
    Mhst = MhsKMM0(j,uai;is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype)

"""

# 0.5D, [nMod]
function MhsKMM0(j::Int64,nai::AbstractVector{T},uai::AbstractVector{T},vthi::AbstractVector{T},
    nMod::Int;is_renorm::Bool=true,rtol_OrjL::T=1e-10,mathtype::Symbol=:Exact) where{T}
    
    if sum_kbn(abs.(uai)) ≤ eps(T)
        return MhsMMM(j,nai,vthi;is_renorm=is_renorm) 
    else
        Mh = 0.0  |> T
        for s = 1:nMod 
            Mh += MhsKMM0(j,nai[s],uai[s],vthi[s];is_renorm=true,rtol_OrjL=rtol_OrjL,mathtype=mathtype) 
        end
        if is_renorm
            return Mh
        else
            return CMjL(T(j)) * Mh
        end
    end
end

"""
  Inputs:
    jvec:
    mathtype: [:Exact, :Taylor0, :Taylor1, :TaylorInf]

  Outputs
    Mhst = MhsKMM0(j,nai,uai,vthi;is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype)
    Mhst = MhsKMM0(j,uai;is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype)

"""

# 0D, []
function MhsKMM0(j::Int64,nai::T,uai::T,vthi::T;
    is_renorm::Bool=true,rtol_OrjL::T=1e-10,mathtype::Symbol=:Exact) where{T}
    
    if abs.(uai) ≤ eps(T)
        if isone(vthi)
            return MhsMMM(j;is_renorm=is_renorm) * nai
        else
            return MhsMMM(j;is_renorm=is_renorm) * (nai * vthi^j)
        end
    else
        if iseven(j)
            if isone(vthi)
                if j == -2
                    if mathtype == :Exact
                    elseif mathtype == :Taylor0
                    else
                        ygtrttt
                        if mathtype == :Taylor1
                        elseif TaylorInf
                        else
                            aswdfcgv
                        end
                    end
                else
                    a = 1.0  |> T
                    for k in 1:Int(j/2)
                        a += CjLk(T(j),T(k)) * uai^(2k)
                    end
                end
                if is_renorm
                    return a * nai
                else
                    return CMjL(T(j)) * a * nai
                end
            else
                uhh = uai/vthi 
                if j == -2
                    if mathtype == :Exact
                        ddbhbhh
                    elseif mathtype == :Taylor0
                        dhbdhn
                    else
                        ygtrttt
                        if mathtype == :Taylor1
                        elseif TaylorInf
                        else
                            aswdfcgv
                        end
                    end
                else
                    # @show j
                    a = 1.0 |> T
                    for k in 1:Int(j/2)
                        a += CjLk(T(j),T(k)) * (uhh)^(2k)
                    end
                end
                if is_renorm
                    return a * nai * vthi^j               # * uhh^L 
                else
                    return CMjL(T(j)) * a * nai * vthi^j     # * uhh^L 
                end
            end
        else
            sdfgbhnm
            if mathtype == :Exact
            elseif mathtype == :Taylor0
            else
                ygtrttt
                if mathtype == :Taylor1
                elseif TaylorInf
                else
                    aswdfcgv
                end
            end
        end
    end
end

# j = L = 0
function MhsKMM0(nai::T) where{T}
    
    return nai[:]
end
"""
"""

# nMod = 1 -> nai = 1, vthi = 1

# 2D, [njMs,ns]
function MhsKMM0!(Mhst::AbstractArray{T},jvec::Vector{Int},uai::T,ns::Int64;
    is_renorm::Bool=true,rtol_OrjL::T=1e-10,mathtype::Symbol=:Exact) where{T}

    for isp in 1:ns 
        MhsKMM0!(Mhst[:,isp],jvec,uai;
                is_renorm=is_renorm,rtol_OrjL=rtol_OrjL,mathtype=mathtype) 
    end
end

# 1D, [njMs]
function MhsKMM0!(Mhst::AbstractVector{T},jvec::Vector{Int},uai::T;
    is_renorm::Bool=true,rtol_OrjL::T=1e-10,mathtype::Symbol=:Exact) where{T}

    k = 0 
    for j in jvec 
        Mhst[k+1] = MhsKMM0(j,uai;is_renorm=is_renorm,mathtype=mathtype)
    end
end

# 0D, []   
function MhsKMM0(j::Int,uai::T;is_renorm::Bool=true,rtol_OrjL::T=1e-10,mathtype::Symbol=:Exact) where{T}
    
    if abs(uai) ≤ eps(T)
        return MhsMMM(j)
    else
        if iseven(j)
            if j == -2
                if mathtype == :Exact
                elseif mathtype == :Taylor0
                else
                    ygtrttt
                    if mathtype == :Taylor1
                    elseif TaylorInf
                    else
                        aswdfcgv
                    end
                end
            else
                # @show j
                a = 1.0  |> T
                for k in 1:Int(j/2) 
                    a += CjLk(T(j),T(k)) * uai^(2k)
                end
            end
            if is_renorm
                return a
            else
                return CMjL(T(j)) * a
            end
        else
            if mathtype == :Exact
            elseif mathtype == :Taylor0
            else
                ygtrttt
                if mathtype == :Taylor1
                elseif TaylorInf
                else
                    aswdfcgv
                end
            end
        end
    end
end
