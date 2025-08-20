
# include(joinpath(pathroot,"Mathematics/maths.jl"))

maxIterKing = 500
# is_C = false                                             # maybe be conservative when `RDnuT` is moderate
# is_C = true                                              # maybe not be conservative when `RDnuT` is bigger enough

show_trace = false
# show_trace = true

# (factor, factor_abbr) = (LeastSquaresOptim.QR(), :QR)     # `=QR(), default`  # More stability
# factor = LeastSquaresOptim.Cholesky()
# factor = LeastSquaresOptim.LSMR()                       # 最差

# is_Jacobian = true                                       # maybe not be conservative when `RDnuT` is bigger enough
# # is_Jacobian = false
# is_Hessian = false
# is_constraint = false
# is_re_seed = false
# is_re_seed = true

# is_renorm = true
is_renorm = false


is_plot_MhjL = false
# is_plot_MhjL = true
is_logplot_MhjL = true

# 对 `(uhh^L)` 再次归一化
# 先格式预测，再非守恒优化，最后再守恒优化。
# 同时测试QR与Cholesky，同时比较is_Jacobian与否，择优而行。
# 若非理想，则减小步长，重新寻优。
# `nMod` 越大，优化的精度越低。

if is_change_datatype 
        include(joinpath(pathroot,"Mathematics/consts_datatype.jl"))
end

if is_re_seed
        naiL = rand(datatype,nModL)
        naiL /= sum_kbn(naiL)
        # naiL /= sum_kbn(naiL)
        # naiL = naiL / sum_kbn(naiL)
        # sum_kbn(naiL)-1, sum_kbn(naiL)-1
        
        
        vthiL = rand(datatype,nModL)
        vthiL /= sum_kbn(vthiL)
        vthiL *= nModL
        
        # L_limit = 45
        # if Vsymmtry == :sphericalMMM
        #         uaiL = zeros(datatype,nModL)
        #         uaiL /= (maximum(abs.(uaiL)) * 5)
        #         L_limit = 0
        # else
        #         uaiL = randn(datatype,nModL)
        #         uaiL /= (maximum(abs.(uaiL)) * 5)
        #         if Vsymmtry == :spherical
        #                 L_limit = 0
        #         end
        # end
        if is_MMM 
                uaiL = zeros(datatype,nModL) 
        else
                uaiL = rand(datatype,nModL) / 10000
                uaiL /= (maximum(abs.(uaiL)) * 5)
        end
end
# if datatype ≠ Float64
#         naiL = datatype.(naiL)
#         uaiL = datatype.(uaiL)
#         vthiL = datatype.(vthiL)
# end
if datatype ≠ Float64
        rtol_OrjL = datatype(rtol_OrjL)
end

kk = 4                  # Hankel 矩阵的阶数，满足：2kk - 1 ≤ njML
njML = njL + 10          # 第 ℓ 阶振幅的动理学矩的个数，满足：njML ≥ 2kk - 1
if njML < 2kk
    njML = 2kk
end
djL = 1
jvec = 0:djL:djL*(njML-1) |> Vector{Int}
jvecL = jvec .+ L

mathtype = :Exact       # [:Exact, :Taylor0, :Taylor1, :TaylorInf]

if is_MMM
        uhLN = 0.0 |> datatype
else
        uhLN = uhLNorm(naiL,uaiL,L)
        # uhLN = maximum(abs.(uaiL))
        
        # uhLN *= 1.3
        # uhLN = 1.0 |> datatype
         
        uhLNL = uhLN .^L
end

MhjL = zeros(datatype,njML)
MhsKMM!(MhjL,jvecL,L,naiL,uaiL,vthiL,uhLN,nModL;is_renorm=is_renorm,is_norm_uhL=is_norm_uhL,rtol_OrjL=rtol_OrjL,mathtype=mathtype)

MhjL = 1.0 * factorial.(0:njML-1)        # f(x) = e⁻ˣ, djL = 1, Mⱼ=factorial(j-1) where j ∈ 𝐍⁺
############# 计算Hanskel矩阵
# njML ≥ 2 * kk -1 || ArgumentError("The number of moments 'nj' must larger than '2k-1'!")
hk,hk1 = zeros(datatype,kk,kk), zeros(datatype,kk,kk)
Hk1type=:modify
# Hk1type=:shift
hankelmatrix!(hk,hk1,MhjL, njML, kk;Hk1type=:modify)

############# 验证Carleman条件
if djL == 1
        CarlM = zeros(datatype,floor(Int,njML/2))
elseif djL == 2
        CarlM = zeros(datatype,njML)
else
        ewgeh
end
carleman!(CarlM,MhjL,jvec, njML;dj=djL)
# display( plot(jvecL, CarlM, xlabel="j", ylabel="Carleman"))
if djL == 1
        display( plot(jvecL[3:2:end], diff(CarlM), xlabel="j", ylabel="diff(Carleman)"))
elseif djL == 2
        display( plot(jvecL[2:end], diff(CarlM), xlabel="j", ylabel="diff(Carleman)"))
else
end

# @show kk, (2kk -1), 3nModL, njML
if Hk1type == :shift 
    if isposdef(hk) ≠ true || isposdef(hk1) ≠ true
            error("Error: Hankel and its shifted matrices must be positive definite.")
    end
else
    if isposdef(hk) ≠ true
            error("Error: Hankel matrices must be positive definite.")
    end
end
@show isposdef(hk), isposdef(hk1)
if datatype == BigFloat
    @show cond(Float64.(hk)), cond(Float64.(hk1))
else
    @show cond(hk), cond(hk1)
end

# is_test_posdef = true
# if is_test_posdef
    # 计算Hankel矩阵的行列式
    is_Jk_method = :ldl
    is_Jk_method = :det
    if is_Jk_method == :det
        Dk = zeros(datatype,kk,3)
        if Hk1type == :shift 
            for k in 1:kk
                @show k, isposdef(hk[1:k,1:k]), isposdef(hk1[1:k,1:k])
                Dk[k,1] = det(hk[1:k,1:k])         # Hk的行列式       
                Dk[k,2] = det(hk1[1:k,1:k])        # Hk'的行列式
                Dk[k,3] = Dk[k,2] ./ Dk[k,1]       # 
            end
        else
            for k in 1:kk
                @show k, isposdef(hk[1:k,1:k])
                hk11 = zeros(datatype,k,k)
                hankel1matrix!(hk11,MhjL, njML, k)
                Dk[k,1] = det(hk[1:k,1:k])         # Hk的行列式       
                Dk[k,2] = det(hk11[1:k,1:k])        # Hk*的行列式
                Dk[k,3] = Dk[k,2] ./ Dk[k,1]       # 
            end
        end
        # Dk[k,3]
    
        # 计算递推公式系数，`ak`与`√bk`
        ak = zeros(datatype,kk)
        bk = zeros(datatype,kk)                # bk = √(bk)
        k = 1
        ak[k] = Dk[k,3]
        # akn1 = Dk[k,2] / Dk[k,1]
        # ak[k] = akn1
        if kk ≥ 2
            k = 2
            ak[k] = Dk[k,3] - Dk[k-1,3]
            # ak[k] = Dk[k,2] / Dk[k,1] - akn1
            # ak[k] = - akn1
            # akn1 = Dk[k,2] / Dk[k,1]
            # ak[k] += akn1
        
            bk[k] = √(Dk[k,1]) / Dk[k-1,1]
            if kk ≥ 3
                for k in 3:kk
                    ak[k] = Dk[k,3] - Dk[k-1,3]
                #     ak[k] = - akn1
                #     akn1 = Dk[k,2] / Dk[k,1]
                #     ak[k] += akn1
        
                    bk[k] = √(Dk[k,1] * Dk[k-2,1]) / Dk[k-1,1]
                end
            end
        end
        
        # ak2 = [1.0,3.0,5.0,7.0]    # e^{-x}
        # # Jacobi 矩阵
        Jk = zeros(datatype,kk,kk)
        k = 1
        Jk[k,k] = ak[k]
        if kk ≥ 2
            Jk[k,k+1] = bk[k+1]
            for k in 2:kk-1
                Jk[k,k-1] = bk[k]
                Jk[k,k] = ak[k]
                Jk[k,k+1] = bk[k+1]
            end
            k = kk
            Jk[k,k-1] = bk[k]
            Jk[k,k] = ak[k]
        end
    elseif is_Jk_method == :ldl     
        LDhk = ldlt(SymTridiagonal(hk))            # ldlt 程序实现方式不准确
        Lhk = LDhk.L
        Dhk = LDhk.D
        Shk = inv(Lhk)
        Jk = abs.(Shk * hk1 * Shk')
    else
        segdvb
    end

    if isposdef(Jk) ≠ true
        #     error("Error: Jacobi matrices must be positive definite.")
    end
    if datatype == BigFloat
        @show cond(Float64.(Jk)), det(Float64.(Jk))
    else
        @show cond(Jk), det(Jk)
    end

    # # # 计算特征值`λ`以及特征向量`V`
#     λ5 = eigvals(Float64.(Jk))                 # 特征值 (复数，BigFloat), for general non-symmetric matrices
#     V5 = eigvecs(Jk)
    # # λ, V = eigen(Jk)
    λ, V = eigen(Symmetric(Jk))
#     @show  λ - λ5


    # # # 验证特征方程
    errJk = Jk * V - V * diagm(λ)
    @show norm(errJk)

    
    # # # 计算特征值对应的权重
    w = V[1,:] .^2  * MhjL[1]
    if datatype == BigFloat
        @show Float64.(w) 
    else
        @show w
    end

    # 离散测度
#     λ = [0.5,2.0,10.0,50.0]
#     w = [0.999,1e-3,1e-6,1e-12]


    # 验证动理学矩的离散积分

    Mcal = zeros(datatype,njML)
    k = 1
    j = 0
    Mcal[k] = sum(w)
#     if kk ≥ 2
        k = 2
        j = 1
        Mcal[k] = sum(w .* λ.^j)
        for k in 3:njML
            jk = (k-1)
            Mcal[k] = sum(w .* λ.^jk)
        end
#     end
    if datatype == BigFloat
        RerrM = Float64.(Mcal ./ MhjL) .- 1
    else
        RerrM = Mcal ./ MhjL .- 1
    end
    RerrM

# end

# if 1 == 2

#         Diff_MhjL = Float64.(diff(MhjL[2:end]))
#         rate_MhjL = Float64.(maximum(abs.(MhjL)) / minimum(abs.(MhjL)))
        
#         if is_show_nuTi 
#             @show Float64.(uhLN), Float64.(uhLNL), rate_MhjL
#             nuTMatrix = [naiL uaiL vthiL uaiL./vthiL]  
#             nuTName = ["naiL", "uaiL", "vthiL", "uhhL"]
#             nuT = DataFrame(Float64.(nuTMatrix),:auto)
#             rename!(nuT,nuTName)  
#             @show nuT
#         end
        
#         if is_plot_MhjL
#                 label = string("sign,RMhjL=",fmtf2(rate_MhjL))
#                 pDiff_MhjL = plot(sign.(Diff_MhjL),label=label)
#                 ylabel!("diff(MhjL)")
        
#                 MhjL_sign = sign.(MhjL[2:end])
#                 MhjL_log = log.(abs.(MhjL[2:end]))
#                 label = string("jM,L=",(jvecL[end],L))
#                 pMhjLlog = plot(MhjL_sign .* MhjL_log,label=label)
#                 # xlabel!("j")
#                 ylabel!("log(|MhjL|)")
        
#                 pDMhjLlog = plot(diff(MhjL_log),label=label)
#                 xlabel!("j")
#                 ylabel!("D(log(|MhjL|))")
        
#                 label = string("jM,L=",(jvecL[end],L))
#                 pMhjL = plot(MhjL[2:end],label=label)
#                 xlabel!("j")
#                 ylabel!("MhjL")
#                 display(plot(pDiff_MhjL, pMhjLlog, pMhjL, pDMhjLlog,layout=(2,2)))
        
#                 @show fmtf2.(MhjL)
#                 @show fmtf2.(MhjL_log)
#                 @show fmtf2.(diff(MhjL_log))
            
#         end
#         # println()
#         # Msnnt = zeros(datatype,njML)
#         # Msnnt = MsnntL2fL(Msnnt,njML,L,naiL,uaiL,vthiL,nModL;is_renorm=is_renorm)
        
#         if is_optim 
#                 Nspan_optim_nuTi = [1.0,1.0, 1.0] |> Vector{datatype}
#                 DMh024 = [1.0,1.0, 1.0] |> Vector{datatype}
#                 rtol_OrjL = 1e-15 |> datatype
#                 atol_Mh = 1e-15 |> datatype
#                 rtol_Mh = 1e-15 |> datatype
                
#                 is_show_Dc = false
#                 # is_show_Dc = true
#                 println("........................................................................")
        
#                 include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_optim_kerl.jl"))
                
#                 # @warn("Checking keywords of optimization, `is_C`, `is_norm_uhL`, `is_constraint`, `is_bs` and `is_Jacobian`, respectively.")
#         end
# end
        
        