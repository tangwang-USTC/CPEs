
"""
  Coefficient of the normalized kinetic moment `𝓜ⱼ₀`
  when `f̂₀(v̂)` is approximated by the KMM0 (including MMM):

    `CMj0 = 2 / √π * Γ((j+3)/2)`.
"""

# # CMj0: The theory values based on Maxwellian model (MM)
# Mhstt(j) = 2 / sqrtpi * gamma((3+j)/2)       # When `j > 30`, the errors may be so large that cannot be ignored.
# dMhstt(j) = - 2 / sqrtpi * 2 / (2+j) * gamma((4+j)/2)
# ddMhstt(j) = 2 / sqrtpi * 2 / (1+j) * gamma((3+j)/2)

#  CMjL(j) = 2 / sqrtpi * gamma((3+j)/2) 
function CMjL(j::Int)

  if j == -2
      return 2
  elseif j == -1
      return 2.0 / sqrtpi
  elseif j == 0
      return 1
  elseif j == 1
      return 2.0 / sqrtpi
  elseif j == 2
      return 1.5
  elseif j == 3
      return 4.0 / sqrtpi
  else
      if iseven(j)
          return prod(3:2:j+1) / 2^(j/2)
          # if j ≤ 30
          #     return 2 / sqrtpi * gamma((3+j)/2)
          # end
      else
          return 2.0 / sqrtpi * prod(2:1:(j+1)/2)
          # return 2 / sqrtpi * gamma((3+j)/2)
      end
  end
end
# f(j) = [2 / sqrtpi * gamma((3+j)/2) / (prod(3:2:j+1) / 2^(j/2))-1, 2 / sqrtpi * gamma((3+j)/2) / (2.0 / sqrtpi * prod(2:1:(j+1)/2))-1]


"""
  Coefficient of the normalized kinetic moment `𝓜ⱼₗ`
  when `f̂ₗ(v̂)` is approximated by the KMM.

    `CMjL = (j+L+1)!! / (2L-1)!! / 2^((j-L)/2), j ∈ L:2:N⁺`.

    `iseven(j - L) == true`
"""

function CMjL(j::Int,L::Int)
    
    if L == 0
        return CMjL(j)
    else
        if iseven(j - L)
            return prod((2L+1):2:(j+L+1)) / 2.0^((j-L)/2)
        else
            edfgbn
        end
    end
end

"""
  Coefficient ofin the Characteristic Parameter Equations (CPEs)
  when `f̂₀(v̂)` is approximated by the KMM:

  `iseven(j) == true`
"""

CjLk(j,L,k) = 2^k * binomial((j-L)/2,k) / prod((2L+3):2:(2k+1)) 

# `L = 0`
CjLk(j,k) = 2^k * binomial(j/2,k) / prod(3:2:(2k+1)) 