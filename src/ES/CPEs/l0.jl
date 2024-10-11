
"""
  Characteristics Parameter Equations (CPEs) for 0D-1V VFP equation when the velocity space is spherically symmetric.

    This form is for ModelingToolkit.

    In this situation, the order of harmonic of distribution is still zero which means the normalized kinetic moment will be:

    `M̂ⱼₗᵐ(r,t) = δₗ⁰δₘ⁰ × Mⱼₗᵐ(fₗᵐ)`
"""

"""
  nMod = 1

"""

# [j, ns]
function CPE0D1V!(eqs::Vector{Vector{Equation}},Mhj0::Vector{AbstractVector{T}}, 
    nh::AbstractVector{T}, vhth::AbstractVector{T}, jvec::Vector{Vector{Int}}, nj::Vector{Int}) where {T<:Real}

    # for isp in 1:ns
    #     for k in 1:nj[isp]
    #         eqs[isp][k] = CPE0D1V(Mhj0[isp][k],nh[isp],vhth[isp],jvec[isp][k])
    #     end
    # end
    for isp in 1:ns
        CPE0D1V!(eqs[isp],Mhj0[isp],nh[isp],vhth[isp],jvec[isp],nj[isp])
    end
end


# [j]
function CPE0D1V!(eqs::Vector{Equation},Mhj0::AbstractVector{T}, nh::T, vhth::T, jvec::Vector{Int}, nj::Int) where {T<:Real}

    for k in 1:nj
        eqs[k] = CPE0D1V(Mhj0[k],nh,vhth,jvec[k])
    end
end


# []
function CPE0D1V(Mhj0::T, nh::T, vhth::T, j::Int) where {T<:Real}

    if iseven(j)
        if j == -2
            error("Procedure is to be completed!")
        elseif j == 0
            eqs = [ 0 ~ nh - Mhj0]
        else
            eqs = [ 0 ~ (nh .* vhth .^ j) - Mhj0]
        end
    else
        error("Procedure is to be completed!")
    end
    return eqs
end


"""
  nMod ≥ 2

"""
# nMod ≥ 2


# [j, nMod, ns]
function CPE0D1V!(eqs::Vector{Vector{Equation}},Mhj0::Vector{AbstractVector{T}}, 
    nh::Vector{AbstractVector{T}}, vhth::Vector{AbstractVector{T}}, jvec::Vector{Vector{Int}}, nj::Vector{Int}) where {T<:Real}

    # for isp in 1:ns
    #     for k in 1:nj[isp]
    #         eqs[isp][k] = CPE0D1V(Mhj0[isp][k],nh[isp],vhth[isp],jvec[isp][k])
    #     end
    # end
    for isp in 1:ns
        CPE0D1V!(eqs[isp],Mhj0[isp],nh[isp],vhth[isp],jvec[isp],nj[isp])
    end
end


# [j,nMod]
function CPE0D1V!(eqs::Vector{Equation},Mhj0::AbstractVector{T}, nh::AbstractVector{T}, vhth::AbstractVector{T}, jvec::Vector{Int}, nj::Int) where {T<:Real}

    for k in 1:nj
        eqs[k] = CPE0D1V(Mhj0[k],nh,vhth,jvec[k])
    end
end

# [nMod]
function CPE0D1V(Mhj0::T, nh::AbstractVector{T}, vhth::AbstractVector{T}, j::Int) where {T<:Real}

    if iseven(j)
        if j == -2
            error("Procedure is to be completed!")
        elseif j == 0
            eqs = [ 0 ~ sum(nh) - Mhj0]
        else
            eqs = [ 0 ~ sum(nh .* vhth .^ j) - Mhj0]
        end
    else
        error("Procedure is to be completed when 'j' is odd!")
    end
    return eqs
end
