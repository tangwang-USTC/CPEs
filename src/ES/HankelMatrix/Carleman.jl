
"""
  Carleman condition:
    
    `∑ⱼ₌₁^∞ (M₂ⱼ^(-1/(2j))) = + ∞  `

  Inputs:
    CarlM: = zeros(floor(Int,nj/2)), if dj = 1
           = zeros(nj), if dj = 2
    MjL: 矩序列
    nj: 矩的个数

  Output:
    carleman!(CarlM,MjL,jvec, nj;dj=dj)
    
"""

function carleman!(CarlM::AbstractVector{T},MjL::AbstractVector{T},jvec::AbstractVector{Int}, nj::Int;dj::Int=1) where {T}

    if dj == 1
        k = 1
        CarlM[k] = MjL[k] ^ (- 1 / jvec[k])
        for k in 2:floor(Int,nj/2)
            j2 = 2k - 1
            CarlM[k] = CarlM[k-1] + MjL[j2] ^ (- 1 / jvec[j2])
        end
    elseif dj == 2
        k = 1
        CarlM[k] = MjL[k] ^ (- 1 / jvec[k])
        for k in 2:nj
            CarlM[k] = CarlM[k-1] + MjL[k] ^ (- 1 / jvec[k])
        end
    else
        wegrgnh
    end 
end





