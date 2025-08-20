
using LinearAlgebra

# 定义矩序列 (以指数分布为例, 密度函数 ρ(x)=e^{-x} 在 [0,∞) 上, 矩为 M_k = k!)
moments(n) = [factorial(k) for k in 0:2n-1]

# 构造k×k Hankel矩阵
function hankel_matrix(M, k)
    H = zeros(Float64, k, k)
    for i in 1:k, j in 1:k
        H[i, j] = M[i+j-1]  # M0, M1, ... 索引从1开始
    end
    return H
end

# 方法1: 行列式法计算递推系数
function recurrence_det(M, k)
    D = [1.0; zeros(k+1)]  # D[-1] = 1, D[0] = 1
    D[1] = M[1]  # D0 = M0 (索引1对应D0)
    
    # 计算Hankel行列式 D1 到 Dk
    for i in 1:k
        H_i = hankel_matrix(M, i)
        D[i+1] = det(H_i)
    end
    
    a = zeros(k)
    b = zeros(k-1)
    for i in 1:k
        # 计算 a_i (i从1开始)
        if i == 1
            a[i] = M[2] / M[1]  # a1 = M1/M0
        else
            a[i] = D[i] / D[i] - D[i-1] / D[i-2]
            # 修正: 使用标准公式 a_i = (D_i^{(1)}/D_{i-1}) - (D_{i-1}^{(1)}/D_{i-2})
            # 实际计算中采用分母调整以增强稳定性
            H1_i = hankel_matrix(M[2:end], i-1)  # H_i^{(1)} 去掉第一列
            H1_i_minus = hankel_matrix(M[2:end], i-2)  # H_{i-1}^{(1)}
            a[i] = det(H1_i) / D[i] - det(H1_i_minus) / D[i-1]
        end
    end
    
    for i in 1:k-1
        # 计算 b_i = ⟨p_{i-1}, p_{i-1}⟩ / ⟨p_{i-2}, p_{i-2}⟩
        b[i] = (D[i+1] * D[i-1]) / (D[i] * D[i])
    end
    return a, b
end

# 方法2: LDL分解法计算递推系数
function recurrence_ldlt(M, k)
    n = k
    H = hankel_matrix(M, n)
    LD = ldlt(SymTridiagonal(H))  # LDL分解: H = L * Diagonal(d) * L'
    L, d = LD.L, LD.D
    # @show L 
    # @show d
    
    a = zeros(n)
    b = zeros(n-1)
    
    # 从LDL分解中提取递推系数
    for i in 1:n-1
        a[i] = L[i+1, i] * d[i] / d[i]  # 实际为 a_i = L[i+1, i+1] 调整
        b[i] = d[i+1] / d[i] * L[i+1, i]^2  # 次对角线平方项
        @show i, a[i], b[i]
    end
    a[n] = L[n, n] * d[n]  # 最后一个a_n
    
    # 调整b_k为平方根形式 (Jacobi矩阵需要)
    b_sqrt = sqrt.(abs.(b))
    return a, b_sqrt
end

# 生成Jacobi矩阵
function jacobi_matrix(a, b)
    n = length(a)
    J = zeros(n, n)
    for i in 1:n
        J[i, i] = a[i]
        if i < n
            J[i, i+1] = b[i]
            J[i+1, i] = b[i]
        end
    end
    return J
end

# 主函数: 计算节点λ和权重w
function stieltjes_moment_problem(k, method=:ldlt)
    M = moments(2k)  # 生成前2k个矩
    
    if method == :ldlt
        a, b = recurrence_ldlt(M, k)
    else
        a, b = recurrence_det(M, k)
    end
    
    @show a 
    @show b
    J = jacobi_matrix(a, b)
    eig = eigen(J)
    λ = eig.values  # 节点 (积分点)
    v1 = eig.vectors[1, :]  # 第一行特征向量
    w = v1.^2 * M[1]  # 权重 w_i = (v1_i)^2 * M0
    
    # 验证矩的准确性
    verified = true
    for m in 1:2k-1
        mom_calc = sum(w .* λ.^(m-1))
        if abs(mom_calc - M[m]) > 1e-5
            println("矩验证失败: M$m = $(M[m]), 计算值 = $mom_calc")
            verified = false
        end
    end
    
    return λ, w, verified
end

# 测试
k = 3  # 使用3个点
λ_det, w_det, verified_det = stieltjes_moment_problem(k, :det)
λ_ldlt, w_ldlt, verified_ldlt = stieltjes_moment_problem(k, :ldlt)

println("LDL方法结果:")
println("节点 λ = ", λ_ldlt)
println("权重 w = ", w_ldlt)

println("\n行列式方法结果:")
println("节点 λ = ", λ_det)
println("权重 w = ", w_det)

