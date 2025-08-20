

"""
  Coefficient of the normalized kinetic moment `𝓜ⱼₗ`
  when `f̂ₗ(v̂)` is approximated by the KMM.

    `CMjL = (j+L+1)!! / (2L-1)!! / 2^((j-L)/2), j ∈ L:2:N⁺`.

    `iseven(j - L) == true`
"""

function CMjL(j::T,L::T) where {T<:Real}
    
    if L == 0
        if j == L 
            return 1 |> T
            # return CMLL()
        else
            return CMjL(j)
        end
    else
        if j == L 
            return 2L + 1
            # return CMLL(L)
        else
            return (2L+1) * gamma(1.5 + (j + L)/2) / gamma(1.5 + L)
        end
        # return prod((2L+1):2:(j+L+1)) / 2.0^((j-L)/2)
        # if iseven(j - L)
        #     return (2L+1) * gamma((3+j + L)/2) / gamma((3+L)/2)
        #     return prod((2L+1):2:(j+L+1)) / 2.0^((j-L)/2)
        # else
        #     return (2L+1) * gamma((3+j + L)/2) / gamma((3+L)/2)
        #     return prod((2L+1):2:(j+L+1)) / 2.0^((j-L)/2)
        # end
    end
end

CMLL(L::T) where {T<:Real} =  2L + 1
# CMLL() =  1.0

"""
  Coefficient of the normalized kinetic moment `𝓜ⱼ₀`
  when `f̂₀(v̂)` is approximated by the KMM0 (including MMM):

    `CMj0 = 2 / √π * Γ((j+3)/2)`.

    NOTES: `Γ(j)` is correct for any `j`! 
            However, `prod(3:2:j+1) / 2^(j/2)` will numerically oscillate for type of Float64!!! 
"""

# # CMj0: The theory values based on Maxwellian model (MM)
# Mhstt(j) = 2 / sqrtpi * gamma((3+j)/2)       # When `j > 30`, the errors may be so large that cannot be ignored.
# dMhstt(j) = - 2 / sqrtpi * 2 / (2+j) * gamma((4+j)/2)
# ddMhstt(j) = 2 / sqrtpi * 2 / (1+j) * gamma((3+j)/2)

#  CMjL(j) = 2 / sqrtpi * gamma((3+j)/2) 
function CMjL(j::T) where {T<:Real}

    if j == -2
        return 2 |> T
    elseif j == -1
        return 2 / sqrtpi |> T
    elseif j == 0
        return 1 |> T
    elseif j == 1.0
        return 2 / sqrtpi |> T
    elseif j == 2
        return 1.5 |> T
    elseif j == 3
        return 4 / sqrtpi |> T
    else
        return 2 / sqrtpi * gamma(1.5 + j/2) |> T
        #   if iseven(j)
        #       return prod(3:2:j+1) / 2^(j/2)
        #       # if j ≤ 30
        #       #     return 2 / sqrtpi * gamma((3+j)/2)
        #       # end
        #   else
        #       return 2.0 / sqrtpi * prod(2:1:(j+1)/2)
        #       # return 2 / sqrtpi * gamma((3+j)/2)
        #   end
    end
end
# f(j) = [2 / sqrtpi * gamma((3+j)/2) / (prod(3:2:j+1) / 2^(j/2))-1, 2 / sqrtpi * gamma((3+j)/2) / (2.0 / sqrtpi * prod(2:1:(j+1)/2))-1]



"""
  Coefficient of the Characteristic Parameter Equations (CPEs)
  when `f̂₀(v̂)` is approximated by the KMM:

  `iseven(j) == true` where `k=0:(j-L)/2`

  Inputs:
    j:
    L:
    k:  where `k=0:(j-L)/2`
    
  Outputs:
    cjLk = CjLk(j,L,k)
    cjLk = CjLk(j,k)
"""

# CjLk(j::Int,L::Int,k::Int) = 2^k * binomial(Int((j-L)/2),k) / prod((2L+3):2:(2(L+k)+1))  where `k=0:(j-L)/2`

function CjLk(j::T,L::T,k::AbstractVector{T}) where {T<:Real}

    if L == 0
        return CjLk(j,k)
    else
        jL2 = (j - L) / 2 + 1
        cjlk = gamma(jL2) ./ gamma.(jL2 .- k)
        jL2 = 1.5 + L
        cjlk .*= (gamma(jL2) ./ gamma.(jL2 .+ k))
        return cjlk ./ gamma.(1 .+ k)
    end
end

# For `ReverseDiff.jl` where `k=0:(j-L)/2`
function CjLk(j::T,L::T,k) where {T<:Real}

    if L == 0
        return CjLk(j,k)
    else
        jL2 = (j - L) / 2 + 1
        cjlk = gamma(jL2) ./ gamma.(jL2 .- k)
        jL2 = 1.5 + L
        cjlk .*= (gamma(jL2) ./ gamma.(jL2 .+ k))
        return cjlk ./ gamma.(1 .+ k)
    end
end

function CjLk(j::T,L::T,k::T) where {T<:Real}

    if L == 0
        return CjLk(j,k)
    else
        jL2 = (j - L) / 2 + 1
        cjlk = gamma(jL2) / gamma(jL2 - k)
        jL2 = 1.5 + L
        cjlk *= gamma(jL2) / gamma(jL2 + k)
        return cjlk / gamma(1 + k)
    end
end

# `L = 0`
# CjLk(j::Int,k::Int) = 2^k * binomial(Int(j/2),k) / prod(3:2:(2k+1)) where `k=0:j/2`
function CjLk(j::T,k::T) where {T<:Real}

    if j == 0
        # if k == 0
        #     return 1.0
        # elseif k == 1
        # elseif k == 2
        # else
        #     return (1 / gamma(1 - k)) * ((sqrtpi / 2) / gamma(1.5 + k)) / gamma(1 + k)
        # end
        return (1 / gamma(1 - k)) * ((sqrtpi / 2) / gamma(1.5 + k)) / gamma(1 + k)
    else
        j2 = j / 2 + 1
        # if k == 0
        #     return 1.0
        # elseif k == 1
        # elseif k == 2
        # else
        # end
        return (gamma(j2) / gamma(j2 - k)) * ((sqrtpi / 2) / gamma(1.5 + k)) / gamma(1 + k)
    end
end




"""
  Coefficient of the derivatives in the Characteristic Parameter Equations (CPEs)
  when `f̂₀(v̂)` is approximated by the KMM:

  `iseven(j) == true`
"""

# CjLk2k(j::Int,L::Int,k::Int) = 2k * CjLk

function CjLk2k(j::T,L::T,k::T) where {T<:Real}
dddd
    if L == 0
        return CjLk2k(j,k)
    else
        return 2k * CjLk(j,L,k)
    end
end

# `L = 0`
function CjLk2k(j::T,k::T) where {T<:Real}
ololloio
    return 2k * CjLk(j,k)
end








# Cj0k(j::Int) = CjLk.(j,1:Int(j/2))

# LLL = 11
# jjj = LLL + 50
# CjLk(j::Int) = CjLk.(j,LLL,1:Int((j-LLL)/2))
# label = string("j,L=",(jjj,LLL))
# pcjlk = plot(CjLk(jjj),label=label,yscale=:log10)
# xlabel!("k")
# ylabel!("CjLk(j,L,k)")
# display(pcjlk)
# CjLk(jjj)
