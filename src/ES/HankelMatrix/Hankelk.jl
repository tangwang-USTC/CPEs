
"""
  Inputs:
    hk:: zeros(k,k)
    MjL: 矩序列
    nj: 矩的个数
    k: Hankel矩阵的阶数
    Hk1type: ∈ [:modify, :shift]

  Output:
    hankelmatrix!(hk,MjL, nj, k)
    hankelmatrix!(hk,hk1,MjL, nj, k;Hk1type=Hk1type)
    
"""

function hankelmatrix!(hk::AbstractArray{T},MjL::AbstractVector{T}, nj::Int, k::Int) where {T}

    if 2k - 1 > nj
        ArgumentError("The number of moments 'nj' must larger than '2k-1'!")
    else
        for i in 1:k
            for j in 1:k
                hk[i,j] = MjL[i+j-1]
            end
        end
    end 

end


function hankelmatrix!(hk::AbstractArray{T},hk1::AbstractArray{T},MjL::AbstractVector{T},
    nj::Int, k::Int; Hk1type::Symbol=:modify) where {T}

    if 2k - 1 ≥ nj
        ArgumentError("The number of moments 'nj' must larger than '2k-1'!")
    else
        if Hk1type == :shift 
            for i in 1:k
                for j in 1:k
                    hk[i,j] = MjL[i+j-1]
                    hk1[i,j] = MjL[i+j]
                end
            end
        elseif Hk1type == :modify 
            for i in 1:k
                for j in 1:k-1
                    hk[i,j] = MjL[i+j-1]
                    hk1[i,j] = hk[i,j]
                end
            end
            for i in 1:k
                hk[i,k] = MjL[i+k-1]
                hk1[i,k] = MjL[i+k]
            end
        else
            pDgdfhf 
        end
    end 

end

"""
  Modified Hankel matrix

  Output:
    hankel1matrix!(hk1,MjL, nj, k)
"""

function hankel1matrix!(hk1::AbstractArray{T},MjL::AbstractVector{T}, nj::Int, k::Int) where {T}

    if 2k - 1 ≥ nj
        ArgumentError("The number of moments 'nj' must larger than '2k-1'!")
    else
        for i in 1:k
            for j in 1:k-1
                hk[i,j] = MjL[i+j-1]
                hk1[i,j] = hk[i,j]
            end
        end
        for i in 1:k
            hk[i,k] = MjL[i+k-1]
            hk1[i,k] = MjL[i+k]
        end
    end 

end

