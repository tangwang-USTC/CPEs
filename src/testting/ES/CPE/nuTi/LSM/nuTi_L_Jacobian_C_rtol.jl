
using Plots,DataFrames,CSV,Format,OhMyREPL
using ModelingToolkit
using Optimization, OptimizationOptimJL, OptimizationManopt
using OptimizationMOI, Ipopt
using Optim
# using Enzyme
using Zygote
using ReverseDiff
# using Tracker

pathroot = "C:/Users/Administrator/.julia/pkgGithub/CPEs"

setprecision(128)  # 设置如256位精度
# datatype = BigFloat
datatype = Float64
is_change_datatype = true
# is_change_datatype = false

rtol_OrjL = 1e-14
RDnuT = 1e-8

is_re_seed = true
# is_re_seed = false
is_show_nuTi = false
is_optim =  true
Vsymmtry = :spherical           # ∈ [:sphericalMMM, :spherical, :axi, :general]

        
is_C = true                     # Whether utilizing the the model with hard constraints (rigid constraints)
is_C = false                    

# is_norm_uhL = true
is_norm_uhL = false
        
    is_Jacobian = true          # (=true, default) Whether Jacobian matrix will be used to improve the performance of the optimizations.
    is_Jacobian = false
    is_Hessian = false

    is_constraint = true        # Whether utilizing the the model with soft constraints
    is_constraint = false

    is_bs = true                # Fminbox constraints, however, it is not necessarily conducive to the convergence of the solution
    is_bs = false

    is_MTK = false
    is_simplify = true
    is_AD = true
    numMultistart = 1  |> Int
    Ncons = 3  |> Int


    numMultistart ≥ 1 || ArgumentError("numMultistart must not be less than 1")
    Ncons ≥ 1 || (is_constraint = false)

    include(joinpath(pathroot,"Mathematics/maths.jl"))
    include(joinpath(pathroot,"src/ES/ESs.jl")) 
    include(joinpath(pathroot,"test/run_collisions/algorithm/modules.jl"))
    include(joinpath(pathroot,"test/run_collisions/paras_alg_optim.jl"))
 

nModL = 1
is_MMM = true
if is_MMM
    nPara = 2
else
    nPara = 2
end
njL = nPara * nModL
is_anasys_L = false

is_L_even = 1
println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
if is_L_even == 1
    L = 0
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))


    # L = 0
    # include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    # L = 2
    # include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    # L = 4
    # include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    # L = 6
    # include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    # L = 8
    # include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
elseif is_L_even == 0
    L = 1
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 3
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 5
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 7
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 9
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
else
    L = 0
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 1
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 2
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 3
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 4
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 5
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 6
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 7
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 8
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
    L = 9
    include(joinpath(pathroot,"src/testting/ES/CPE/nuTi/LSM/nuTi_L_initial_rtol_kerl.jl"))
end

# @show NL_solve,Optlibary,optimizer,ADtype, (is_MTK,is_C,Ncons,is_constraint,is_bs)

1