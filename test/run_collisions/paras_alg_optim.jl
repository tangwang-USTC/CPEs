
const nhMax = 1e8       # nₑ = 10n₂₀  ~  nfuel = 10⁶n₂₀; nEP = 10¹³. The maximum absolute span of the submoments;
const uhMax = 100.0     # uh ≤ 3.0    ~  L_limit = 45. The maximum absolute span of the submoments;
const uhMin = 1e-4
const vhthMax = 50.0    # TD = 1 keV  ~  Tα = 3 MeV. The maximum absolute span of the submoments;
const vhthMin = 1e-6    # TD = 1 keV  ~  Tα = 3 MeV. The maximum absolute span of the submoments;
const vhthInitial = 4e-1# vthi 
const nhInitial = 2e-1  # nai

## Parameters: Optimization solver
## Algorithm `King`: Optimizations for `fvL`, `jMsmax = ℓ + dj * (3(nMod-1) + (3-1))`.
## Levels: [NL_solve, Optlibary, NL_solve_method, optimizer, linsolve, linesearch, preconditioner, ADtype]
        # 1, NL_solve: ∈ [:NonlinearSolve,:Optimization,:JuMP,:NLPModels] 
        # 2, Optlibary: ∈ [:Optimization,:Optim,:LeastSquaresOptim,:NLSolves,⋯]
        # 3, NL_solve_method: solver ∈ [Newton, Quasi-Newton, Krylov, TrustRegion,⋯]
        # 4, optimizer: ∈ [CG, LBFGS, KrylovTrustRegion, dogleg, LM,⋯]            -- nonlinear solver
        # 5, linsolve: ∈ [:qr,:cholesky,:lsmr,:ILU,:KLU,GMRES,:BCG,:CG,,⋯]        -- factor:linear solver
        # 6, linesearch: ∈ [:HagerZhang, :MoreThuente,:BackTracking,:StrongWolfe,:Static]     #  LineSearches.jl
        # 7, preconditioner: ∈ [:,⋯]
        # 8, ADtype: ∈ DifferentiationInterface.jl or ADTypes.jl
 ## Notes: Quasi-Newtion methods: LBFGS


 linsolve = nothing          # (default = nothing)
 linesearch = nothing        # HagerZhang method if good for cases when `tol ≤ 10⁻⁸`
 preconditioner = nothing

NL_solve = :Optimization
# # NL_solve = :NonlinearSolve
NL_solve = :LeastSquaresOptim
if NL_solve == :LeastSquaresOptim
    ADtype = :forward 
    ADtype = :central

    Optlibary = NL_solve
    NL_solve_method = :trust_region         # [:trust_region, :newton, :anderson]

    optimizer = LeastSquaresOptim.Dogleg
    # optimizer = LeastSquaresOptim.LevenbergMarquardt

    linsolve = LeastSquaresOptim.QR()       # `=QR(), default`  # More stability, maybe not conservative
    # linsolve = LeastSquaresOptim.Cholesky()
    # linsolve = LeastSquaresOptim.LSMR()
    
elseif NL_solve == :Optimization                # with nonlinear inequality/equality constraints
    # Automatic Differentional
    ADtypes = [AutoEnzyme(), AutoForwardDiff(), AutoModelingToolkit(), AutoSparseForwardDiff(), AutoModelingToolkit(true, true), 
               AutoReverseDiff(), AutoZygote(), AutoTracker(), AutoSparseReverseDiff(), AutoSparse(AutoZygote()),
               AutoFiniteDiff()]
    ADtype = nothing 
    # ADtype = AutoEnzyme()              # May be false stems from the type of variables 

    # ADtype = AutoZygote()              # cons_fun!(res::Zygote.Buffer{Float64, Vector{Float64}}, xx::Vector{Float64}, p::Vector{Float64}))
    # ADtype = AutoReverseDiff()
    ADtype = AutoForwardDiff()         # (default) cons_fun!(res::Vector{ForwardDiff.Dual{…}}, xx::Vector{ForwardDiff.Dual{…}}, p::Vector{Float64})

    # # ADtype = AutoTracker()
    # ADtype = AutoModelingToolkit()

    Optlibary = :Optimization
    # Optlibary = :Optim
    # Optlibary = :Polyalgorithms
    # Optlibary = :MathOptInterface
    # Optlibary = :Manopt

    NL_solve_method = :trust_region         # [:trust_region, :newton, :anderson]

    # "L/G" means "Local/Global"
    # "N/D" denotes "derivative-free/gradient-based" algorithm
    # "η=alpha" denotes learning rate, 
    # "ρ" denotes momentum
    # "β=(beta_mean,beta_var)" denotes tuple is the exponential decay of momentums for first and second moments estimates
    if Optlibary == :Optimization 
        optimizer = Optimization.LBFGS   #(J!,H!,cons,(lb,ub)), Quasi-Newtion, Augmented Lagrangian method -> a general-purpose nonlinear optimization solver
        NL_solve_method = :quasi_newton

        #### for training NN such as `Lux` or `Flux`
        # optimizer = Optimization.Sophia #(J!,H!,cons),
        # NL_solve_method = :newton
    elseif Optlibary == :Optim                  #(J!,H!,cons,lb,ub)
        # if is_constraint
            # # derivative-free
            # optimizer = Optim.NelderMead                    #((lb,ub),), direct search: Nelder-Mead
            # optimizer = Optim.SAMIN                         #(lb,ub,Global)), search: Simulated Annealing with bounds constraints
            # optimizer = Optim.ParticleSwarm                 #((lb,ub),Global), Particle Swarm
            # optimizer = Optim.SimulatedAnnealing            #(), Simulated Annealing
    
            # # # gradient-based (J)           # best solution
            # # Adam(alpha=0.0001,beta_mean=0.9,beta_var=0.999) # Adam optimizer
            # # AdaMax(alpha=0.002,beta_mean=0.9,beta_var=0.999)# AdaMax optimizer
            # optimizer = Optim.ConjugateGradient             #(J!,(lb,ub),linesearch), Conjugate Gradient (CG) Descent
            # optimizer = Optim.GradientDescent               #(J!,(lb,ub),linesearch), quasi-Newton: identity Hessian matrix
            # optimizer = Optim.BFGS                          #(J!,(lb,ub),linesearch, manifold), quasi-Newton: approximated Hessian matrix, Broyden-Fletcher-Goldfarb-Shanno
            optimizer = Optim.LBFGS                         #(J!,(lb,ub),linesearch, manifold), quasi-Newton: Limited-memory approximated Hessian matrix
            
            # # # Hessian-free Second Order (Hv)
            # optimizer = Optim.KrylovTrustRegion             #(J!,Hv!,lb,ub), Newton-Krylov method with Trust Region dogleg method
            
            # # # Hessian-Based Second Order (H)
            # optimizer = Optim.Newton                        #(J!,H!,linesearch), first order conditions (∇f(x)=0), stable for globally well-approximated by a quadratic
            # optimizer = Optim.NewtonTrustRegion             #(J!,H!), Newton Trust Region dogleg method with Cholesky decomposition under assumption of locally quadratic.
    
            # # # gradient-based (J), Acceleration methods
            # optimizer = Optim.OACCEL                        #(J!,linesearch, manifold), accelerates based on a minimization of a n approximation to the objective
            # optimizer = Optim.NGMRES                        #(J!,linesearch, manifold), nonlinear preconditioner on a ℓ₂ norm subspace
        # else
            # # Hessian-Based Second Order (H)
            # optimizer = Optim.IPNewton                      #(J!,H!,cons,lb,ub,linesearch), Interior-point Newton with local constraint
            # # optimizer = Optim.IPNewton(Optim.Options(allow_f_increases = true, successive_f_tol = 2)) # local constraint
        # end
        # optimizer = Optim.BFGS
    elseif Optlibary == :NLopt                  #(J!,H!,cons,lb,ub)
        # Local optimizer
        # derivative-free
        optimizer = NLopt.LN_COBYLA()                 # (Constrained Optimization BY Linear Approximations) nonlinear inequality/equality constraints, constructing successive linear approximations of the objective function and constraints
        optimizer = NLopt.LN_PRAXIS()
        optimizer = NLopt.LN_NEWUOA()
        optimizer = NLopt.LN_NEWUOA_BOUND()
        optimizer = NLopt.LN_NELDERMEAD()
        optimizer = NLopt.LN_SBPLX()
        optimizer = NLopt.LN_AUGLAG()
        optimizer = NLopt.LN_AUGLAG_EQ()
        optimizer = NLopt.LN_BOBYQA()

        # gradient-based
        optimizer = NLopt.LD_SLSQP()                   # (Sequential quadratic programming) nonlinear inequality/equality constraints
        optimizer = NLopt.LD_MMA()                     # (Method of Moving Asymptotes), improved CCSAQ for nonlinear inequality constraints
        optimizer = NLopt.LD_CCSAQ()                   # (conservative convex separable approximation)
        optimizer = NLopt.LD_LBFGS_NOCEDAL()
        optimizer = NLopt.LD_LBFGS()
        optimizer = NLopt.LD_VAR1()
        optimizer = NLopt.LD_VAR2()
        optimizer = NLopt.LD_TNEWTON()
        optimizer = NLopt.LD_TNEWTON_RESTART()
        optimizer = NLopt.LD_TNEWTON_PRECOND()
        optimizer = NLopt.LD_TNEWTON_PRECOND_RESTART()
        ######### Augmented Lagrangian algorithm, combining the objective function and the nonlinear inequality/equality constraints
        optimizer = NLopt.LD_AUGLAG()
        optimizer = NLopt.LD_AUGLAG_EQ()             # Only uses penalty functions for equality constraints, while inequality constraints are passed through to the subsidiary algorithm to be handled directly

        # Global optimizer
        # with constraint using the MOI wrapper
        optimizer = NLopt.GN_ISRES()              # (Improved Stochastic Ranking Evolution Strategy) nonlinear inequality/equality constraints
        optimizer = NLopt.GN_ORIG_DIRECT()        # nonlinear inequality constraints
        optimizer = NLopt.GN_AGS()                # nonlinear inequality constraints, Lipschitz condition using posterior probabilities
        optimizer = NLopt.GN_ORIG_DIRECT_L()

        # without constraint
        optimizer = NLopt.GN_DIRECT()
        optimizer = NLopt.GN_DIRECT_L()
        optimizer = NLopt.GN_DIRECT_L_RAND()
        optimizer = NLopt.GN_DIRECT_NOSCAL()
        optimizer = NLopt.GN_DIRECT_L_NOSCAL()
        optimizer = NLopt.GN_DIRECT_L_RAND_NOSCAL()
        optimizer = NLopt.GD_STOGO()
        optimizer = NLopt.GD_STOGO_RAND()
        optimizer = NLopt.GN_CRS2_LM()
        optimizer = NLopt.GN_MLSL()
        optimizer = NLopt.GD_MLSL()
        optimizer = NLopt.GN_MLSL_LDS()
        optimizer = NLopt.GD_MLSL_LDS()
        optimizer = NLopt.G_MLSL()
        optimizer = NLopt.G_MLSL_LDS()
        optimizer = NLopt.GN_ESCH()
    elseif Optlibary == :Optimisers              #(J!,H!,cons,lb,ub), for training NN
        optimizer = Optimisers.Descent(η=0.1)                # Classic gradient descent
        optimizer = Optimisers.Momentum(η=0.01,ρ=0.9)        # Classic gradient descent
        optimizer = Optimisers.Nesterov(η=0.001,ρ=0.9)       # Gradient descent
        # RProp(η,ℓ,Γ)              # RProp, a full-batch learning algorithm that depends only on the sign of the gradient.
        optimizer = Optimisers.RMSProp(0.001,ρ=0.9)          # RMSProp, a good choice for recurrent networks.
        # CentredRMSProp(η,ρ)       # a variant which normalises gradients by an estimate their variance, instead of their second moment.
        optimizer = Optimisers.Adam(η=0.001, β=(0.9, 0.999))   # Adam optimizer
        optimizer = Optimisers.RAdam(η=0.001, β=(0.9, 0.999))  # Rectified Adam
        optimizer = Optimisers.OAdam(η=0.001, β=(0.5, 0.9))    # Optimistic Adam
        optimizer = Optimisers.AdaMax(η=0.001, β=(0.9, 0.999)) # AdaMax
        optimizer = Optimisers.ADAGrad(η=0.1)                      # ADAGrad
        optimizer = Optimisers.ADADelta(ρ=0.9)                     # ADADelta
        optimizer = Optimisers.AMSGrad(η=0.001, β=(0.9, 0.999))    # AMSGrad
        optimizer = Optimisers.NAdam(η=0.001, β=(0.9, 0.999))      # Nesterov variant of the Adam
        optimizer = Optimisers.AdamW(η=0.001, β=(0.9, 0.999), λ=0) # AdamW
        optimizer = Optimisers.ADABelief(η=0.001, β=(0.9, 0.999))  # ADABelief variant of Adam
        # Lion(η=0.001, β=(0.9, 0.999))                # Lion
    elseif Optlibary == :Polyalgorithms          #(J!,H!,cons,lb,ub),   # Adam + BFGS
        optimizer = PolyOpt
    elseif Optlibary == :MathOptInterface        #(J!,H!,cons)
        optimizer = Ipopt.Optimizer              # C++ 
        # optimizer = KNITRO.Optimizer()
        # optimizer = Juniper.Optimizer()
    else
        if Optlibary == :Manopt                  #(J!,lb,ub,Global, Monifold)
            optimizer = OptimizationManopt.NelderMeadOptimizer                 # Nelder-Mead (NM)
            optimizer = OptimizationManopt.QuasiNewtonOptimizer                # quasi-Newton 
            optimizer = OptimizationManopt.GradientDescentOptimizer            # gradient descent (GD)
            optimizer = OptimizationManopt.ConjugateGradientDescentOptimizer   # conjugate gradient (CG) descent 
            optimizer = OptimizationManopt.CMAESOptimizer                      # CMAES, convex bundle method (CBM)
            optimizer = OptimizationManopt.ParticleSwarmOptimizer              # particle swarm Optimization (PSO)
            optimizer = OptimizationManopt.ConvexBundleOptimizer               # 
            optimizer = OptimizationManopt.FrankWolfeOptimizer                 # Frank-Wolfe 
        elseif Optlibary == :Metaheuristics          #((cons), lb,ub,Global), Metaheuristic algorithm
            optimizer = OptimizationMetaheuristics.ECA()         # Evolutionary Centers Algorithm
            optimizer = OptimizationMetaheuristics.DE()          # Differential Evolution
            optimizer = OptimizationMetaheuristics.PSO()         # Particle Swarm Optimization
            optimizer = OptimizationMetaheuristics.ABC()         # Artificial Bee Colony 
            optimizer = OptimizationMetaheuristics.CGSA()        # Gravitational Search Algorithm
            optimizer = OptimizationMetaheuristics.SA()          # Simulated Annealing
            optimizer = OptimizationMetaheuristics.WOA()         # Whale Optimization Algorithm
        elseif Optlibary == :BlackBoxOptim           #(lb,ub,Global), Metaheuristic/stochastic algorithms
            optimizer = ()       # Natural Evolution Strategy (NES)
            optimizer = ()       # Differential Evolution (DE) Strategy
            optimizer = ()       # Direct search (DS)
            optimizer = ()       # Resampling Memetic Searchers
            optimizer = ()       # Stochastic Approximation
            optimizer = ()       # Random Search (RS)
        elseif Optlibary == :Evolutionary            #(lb,ub,Global)
            optimizer = Evolutionary.GA()   # Genetic Algorithm optimizer
            optimizer = Evolutionary.DE()   # Differential Evolution optimizer
            optimizer = Evolutionary.ES()   # Evolution Strategy algorithm
            optimizer = Evolutionary.CMAES()# Covariance Matrix Adaptation Evolution Strategy (CMAES) algorithm
        elseif Optlibary == :CMAEvolutionStrategy    #(lb,ub,Global) Covariance Matrix Adaptation Evolution Strategy (CMAES) algorithm.
            optimizer = CMAEvolutionStrategyOpt() #
        elseif Optlibary == :GCMAES                  #(J!,lb,ub,Global), Gradient-based Covariance Matrix Adaptation Evolutionary Strategy (CMAES)
            optimizer = GCMAESOpt()
        elseif Optlibary == :MultistartOptimization  #(lb,ub,Global)
            n_initial = 10
            optimizer = MultistartOptimization.TikTak(n_initial)
        elseif Optlibary == :NOMAD                   #(lb,ub,Global), C++, Mesh Adaptive Direct Search algorithm (MADS), designed for difficult blackbox optimization problems
            optimizer = NOMADOpt()
        elseif Optlibary == :SpeedMapping           #(J!,lb,ub) without constraint
            optimizer = SpeedMappingOpt()
        else
            dsfgfgf
        end
    end
    # optimizer = Optim.KrylovTrustRegion()
    # optimizer = Optim.Newton()            # Without bounds "lbs" and "ubs

elseif NL_solve == :NonlinearSolve 
    # # # # # # Valid choices are types from ADTypes.jl
    ADtype = nothing       # (default = nothing), meaning that a default is selected according to the problem specification!
    ADtype = :central    # A little more complicated but maybe obtaining no benefits.
    # [AutoForwardDiff(), AutoEnzyme(; mode=Enzyme.Forward),AutoZygote(),AutoSymbolics(), AutoSparse(AutoSymbolics()),AutoForwardDiff(),AutoFastDifferentiation(), AutoSparse(AutoFastDifferentiation()),AutoPolyesterForwardDiff(; tag=:hello),AutoMooncake()]
    # backends = [AutoReverseDiff(; compile=false), AutoReverseDiff(; compile=true),AutoFiniteDiff(),AutoFiniteDifferences(; fdm=FiniteDifferences.central_fdm(3, 1))]
    # second_order_backends = [SecondOrder(AutoForwardDiff(), AutoReverseDiff())]

    linsolve = nothing # (default = nothing), meaning that a default is selected according to the algorithm
    if Optlibary == :LeastSquaresOptim
        optimizer = :dogleg                      # Newton Trust Region dogleg method
        optimizer = :lm                          # Advanced Levenberg-Marquardt 
        # :newton
        # :anderson

        linsolve = :qr                      # default, LeastSquaresOptim.QR()
        linsolve = :cholesky
        linsolve = :lsmr

        ADtype = :central
        ADtype = :forward 
        
        algnl = LeastSquaresOptimJL(optimizer;linsolve=linsolve,autodiff=ADtype)
    elseif Optlibary == :NonlinearSolve   
        ##### linear solvers : LinearSolve.jl
        ## 1. `A` is Number, `u = b / A`
        ## 2. `A` is SMatrix, `u = b \ A`
        ## 3. `linsolve` is `\`, `\div!(u,A,b)`
        ## 4. `linsolve` is `alg` in LinearSolve.jl, we utilize `alg` to solve the linear system

        #### RecursiveFactorization.jl
        linsolve = RFUFactorization()                     # (Float64, fastest for dense, smaller matrices (<500x500)) Fast pure Julia LU-factorization.

        #### Base.LinearAlgebra, working for many array types, such as CuArrays for GPU-accelerated solving
        # linsolve = (x, A, b) -> copyto!(x, A\b)           # `lu!(A)`
        linsolve = LUFactorization()                      # `lu!(A)`
        linsolve = GenericLUFactorization()               # `generic_lufact!(A)` Supports arbitrary number types and good for small matrices
        linsolve = QRFactorization()                      # `qr!(A)` more precision for dense matrix, sparse matrix, CuMatrix, BandedMatrix and BlockBandedMatrix
        linsolve = SVDFactorization()                     # `svd!(A)` slowest but most precise for dense matrix, 
        linsolve = CholeskyFactorization()
        linsolve = BunchKaufmanFactorization()
        linsolve = CHOLMODFactorization()
        linsolve = NormalCholeskyFactorization()
        linsolve = NormalBunchKaufmanFactorization()

        #### LinearSolve.jl
        linsolve = SimpleLUFactorization()
        linsolve = DiagonalFactorization()
        linsolve = SimpleGMRES()

        #### FastLapackInterface.jl
        linsolve = FastLUFactorization()
        linsolve = FastQRFactorization()

        #### SuiteSparse.jl
        linsolve = KLUFactorization()
        linsolve = UMFPACKFactorization()

        #### Sparspak.jl
        linsolve = SparspakFactorization()

        #### Krylov.jl
        # KrylovJL(args...; KrylovAlg = Krylov.gmres!,Pl = nothing, Pr = nothing,gmres_restart = 0, window = 0,kwargs...)
        linsolve = KrylovJL_CG()
        linsolve = KrylovJL_MINRES()
        linsolve = KrylovJL_GMRES()
        linsolve = KrylovJL_BICGSTAB()
        linsolve = KrylovJL_LSMR()
        linsolve = KrylovJL_CRAIGMR()

        #### IterativeSolvers.jl
        linsolve = IterativeSolversJL()
        linsolve = IterativeSolversJL_CG()
        linsolve = IterativeSolversJL_MINRES()
        linsolve = IterativeSolversJL_GMRES()
        linsolve = IterativeSolversJL_BICGSTAB()
        linsolve = IterativeSolversJL_IDRS()

        #### KrylovKit.jl
        linsolve = KrylovKitJL()
        linsolve = KrylovKitJL_CG()
        linsolve = KrylovKitJL_GMRES()

        #### HYPRE.jl
        linsoltype1 = HYPRE.BiCGSTAB
        linsoltype1 = HYPRE.BoomerAMG
        linsoltype1 = HYPRE.FlexGMRES
        linsoltype1 = HYPRE.GMRES
        linsoltype1 = HYPRE.Hybrids
        linsoltype1 = HYPRE.ILU
        linsoltype1 = HYPRE.ParaSails                # (as preconditioner only)
        linsoltype1 = HYPRE.PCG

        linsolve = HYPREAlgorithm(linsoltype1)

        #### MKL.jl
        linsolve = MKLLUFactorization()

        #### CUDA.jl
        linsolve = CudaOffloadFactorization()

        #### Pardiso.jl
        linsolve = PardisoJL()
        linsolve = MKLPardisoFactorize()
        linsolve = MKLPardisoIterate()

        #### Metal.jl
        linsolve = MetalLUFactorization()

        #### AppleAccelerate.jl
        linsolve = AppleAccelerateLUFactorization()
        



        # nls_algs = (NewtonRaphson(), TrustRegion(), LevenbergMarquardt(), PseudoTransient(), 
        #             GeneralBroyden(), GeneralKlement(), DFSane(), nothing)
        # # # (linsolve=nothing,linesearch=missing,autodiff=nothing, 
        # # #  vjp_autodiff=nothing,jvp_autodiff=nothing,concrete_jac=nothing)

        # Adavanced solvers: General source schemes for nonlinear solvers presented in followsing sections.
        NonlinearSolveQuasiNewton.QuasiNewtonAlgorithm            # Using a Trust Region Method. 
        NonlinearSolveFirstOrder.GeneralizedFirstOrderAlgorithm   # Uses Jacobian.
        NonlinearSolveSpectralMethods.GeneralizedDFSane           # Jacobian-Free Spectral Method

        ##### Nonlinear solvers (NS)
        ############################## for large-scale and numerically-difficult nonlinear systems
        NonlinearSolveFirstOrder.NewtonRaphson             # Advanced NewtonRaphson, supporting for efficient handling of sparse matrices via colored AD and preconditioned linear solverss for high performance on large and sparse systems.
        ############################## for solving large-scale nonlinear systems of equations
        NonlinearSolveSpectralMethods.DFSane               # Low-overhead and allocation-free implementation of the df-sane method, 
        ############################## 
        NonlinearSolveQuasiNewton.Broyden                  # Broyden's Method [3] with resetting and line search. Fast method but unstable when the condition number of the Jacobian matrix is sufficiently large.
        NonlinearSolveQuasiNewton.LimitedMemoryBroyden
        NonlinearSolveQuasiNewton.Klement                  # Implementation of Klement [4] with line search, preconditioning and customizable linear solves. 
        
        ##### Nonlinear least squares solvers (NLSS)
        NonlinearSolveFirstOrder.GaussNewton               # Advanced GaussNewton for large-scale and numerically-difficult nonlinear systems.

        
        ##### Nonlinear (least squares) solvers (NS/NLSS) default method：
        ############################## for large-scale and numerically-difficult nonlinear systems
        NonlinearSolveFirstOrder.TrustRegion               # Advanced TrustRegion method dogleg method, handling of sparse matrices via colored AD and preconditioned linear solversfor high performance on large and sparse systems.
        NonlinearSolveFirstOrder.LevenbergMarquardt        # Advanced Levenberg-Marquardt implementation with the improvements suggested by Transtrum and Sethna for large-scale and numerically-difficult nonlinear systems
        ############################## To solve steady state problems in an accelerated manner using "switched evolution relaxation" [8] SER method.
        NonlinearSolveFirstOrder.PseudoTransient           # Using an adaptive time-stepping to switch over to Newton's method and gain a rapid convergence. Good for highly unstable systems.

        #### Polyalgorithms: Container for a tuple of NS or NLSS
        NonlinearSolve.RobustMultiNewton                   # Robustness, using a mixture of Newton methods with different globalizing techniques (trust region updates, line searches, etc.) 
        NonlinearSolve.FastShortcutNLLSPolyalg             # Less robust methods for more performance and then tries more robust techniques if the faster ones fail.
        NonlinearSolve.FastShortcutNonlinearPolyalg        # (default) Mixing fast methods with fallbacks to robust methods to allow for solving easy problems quickly without sacrificing robustness on the hard problems.
        NonlinearSolveBase.NonlinearSolvePolyAlgorithm(algs) # (alg1,alg2), Container for a tuple of algorithms for NonlinearProblem and NonlinearLeastSquaresProblem


        custom_polyalg = NonlinearSolvePolyAlgorithm((GeneralBroyden(), LimitedMemoryBroyden()))
    elseif Optlibary == :SimpleNonlinearSolve
        # # # # ADtype is provided in DifferentiationInterface.jl which is determined by ADTypes.jl
        optimizer = SimpleNewtonRaphson(;autodiff=ADtype)
        optimizer = SimpleBroyden(;linesearch = Val(false), alpha = nothing)       # `Val(true)` denotes using `LiFukushimaLineSearch` line search
        optimizer = SimpleHalley(;autodiff=ADtype)
        optimizer = SimpleKlement()
        optimizer = SimpleDFSane()          # df-sane method for solving large-scale nonlinear systems of equations
        # SimpleDFSane(;σ_min::Real = 1e-10, σ_max::Real = 1e10, σ_1::Real = 1.0,
        #               M::Union{Int, Val} = Val(10), γ::Real = 1e-4, τ_min::Real = 0.1, τ_max::Real = 0.5,
        #               nexp::Int = 2, η_strategy::Function = (f_1, k, x, F) -> f_1 ./ k^2
        #             )
        optimizer = SimpleLimitedMemoryBroyden(linesearch = Val(false), alpha = nothing)  # L-BFGS <-> SimpleBroyden
        optimizer = SimpleTrustRegion(;autodiff=ADtype)
        # SimpleTrustRegion(;
        #                   autodiff = AutoForwardDiff(), max_trust_radius = 0.0,
        #                   initial_trust_radius = 0.0, step_threshold = nothing,
        #                   shrink_threshold = nothing, expand_threshold = nothing,
        #                   shrink_factor = 0.25, expand_factor = 2.0, max_shrink_times::Int = 32,
        #                   nlsolve_update_rule = Val(false)
        #                   )

    elseif Optlibary == :DescentSubroutines
        optimizer = NonlinearSolveBase.NewtonDescent(;linsolve=linsolve)         # For non-square Jacobian problems, this is commonly referred to as the Gauss-Newton Descent.
        optimizer = NonlinearSolveBase.SteepestDescent(;linsolve=linsolve)       # The linear solver and preconditioner are only used if Jacobian is provided in the inverted form.
        optimizer = NonlinearSolveBase.DampedNewtonDescent(;linsolve=linsolve, initial_damping, damping_fn)   # A Newton descent algorithm with damping.
        optimizer = NonlinearSolveBase.Dogleg(;linsolve=linsolve)                # Switch between Newton's method and the steepest descent method depending on the size of the trust region.
        optimizer = NonlinearSolveBase.GeodesicAcceleration(;linsolve=linsolve)  # Special Levenberg Marquardt Descent Subroutines: Uses the descent algorithm to compute the velocity and acceleration terms for the geodesic acceleration method. The velocity and acceleration terms are then combined to compute the descent direction.

    elseif Optlibary == :GlobalizationSubroutines  #  Radius Update Schemes
        # # solve(prob,alg=TrustRegion(radius_update_scheme=optimizer))  # 
        # # NonlinearSolveFirstOrder.RadiusUpdateSchemes
        optimizer = RadiusUpdateSchemes.Simple         # (default) Follows the conventional approach to update the trust region radius
        optimizer = RadiusUpdateSchemes.Hei            # Radius depends on the size (norm) of the current step size
        optimizer = RadiusUpdateSchemes.Yuan           # Similar to "Hei" but can converges to zero
        optimizer = RadiusUpdateSchemes.Fan            # Similar to "Yuan", radius depend on the current size (norm) of the objective (merit) function itself
        optimizer = RadiusUpdateSchemes.Bastin         # A retrospective update scheme as it uses the model function at the current iteration to compute the ratio of the actual reduction and the predicted reduction in the previous trial step
        optimizer = RadiusUpdateSchemes.NLsolve        # Trust region dogleg in NLsolve.jl
        optimizer = RadiusUpdateSchemes.NocedalWright  # Given by Nocedal and Wright

    elseif Optlibary == :SteadyStateDiffEq
        optimizer = DynamicSS(Tsits())                # Coupled with OrdinaryDiffEq.jl
        optimizer = SSRootfind()                      # SSRootfind(alg=nothing)
        
    elseif Optlibary == :NLsolve  
        ADtype = :central
        ADtype = :forward 
        ADtype = :forward 

        optimizer = :trust_region       # Trust region Newton method (the default choice)
        optimizer = :newton             # Classical Newton method with an optional line search
        optimizer = :broyden            # Broyden's quasi-Newton method
        optimizer = :anderson           # Anderson-accelerated fixed-point iteration

        NLsolveJL(;method=optimizer,autodiff=ADtype,linesearch=linesearch,
                linsolve = (x, A, b) -> copyto!(x, A\b), factor = one(Float64), autoscale = true,m = 10, beta = one(Float64))
    elseif Optlibary == :SIAMFANLEquations    # C. T. Kelley
        optimizer = :newton             # Classical Newton method.
        optimizer = :pseudotransient    # Pseudo transient method.
        optimizer = :secant             # Secant method for scalar equations.
        optimizer = :anderson           # Anderson acceleration for fixed point iterations.
    else
        if Optlibary == :BracketingNonlinearSolve # Interval Methods  are suited for interval (scalar) root-finding problems
            # # ITP()
            # # Alefeld()
            # # Falsi()
            # # Ridder()
            Brent()
        elseif Optlibary == :FixedPointAcceleration     # `f(x) = x`
            optimizer = :Newton
            optimizer = :Aitken
            optimizer = :VEA
            optimizer = :SEA
            optimizer = :Simple
            optimizer = :RRE                  # Reduced Rank Extrapolation
            optimizer = :MPE                  # Minimal Polynomial Extrapolation
            optimizer = :Anderson             # Anderson (1965) acceleration
        elseif Optlibary == :MINPACK           # C/C++
            optimizer = :hybr                      # Modified version of Powell's algorithm
            optimizer = :lm                        # Levenberg-Marquardt 
            optimizer = :lmdif                     # Advanced Levenberg-Marquardt 
            optimizer = :hybrd                     # Advanced modified version of Powell's algorithm
        elseif Optlibary == :NLSolves   # Without AD
            # Sampling based
            SimulatedAnnealing
            PureRandomSearch
            ParticleSwarm
            
            # Direct search
            NelderMead
    
    
            # Quasi-Newton Line search
            DBFGS
            BFGS
            SR1
            DFP
            GradientDescent
            LBFGS
    
            # Conjugate Gradient Descent
            ConjugateGradient
    
            # Newton Line Search
            Newton
    
            # Gradient based (no line search)
            Adam
            AdaMax
    
            # Acceleration methods
            Anderson
    
            # Krylov
            InexactNewton
    
            # BB style
            BB
            DFSANE
    
            # Simple bounds (box)
            ParticleSwarm
            ActiveBox
            rtfghnb
        else
            ikggfffdsdds
        end
    end
elseif NL_solver == :JuMP 
else
    if 1 == 2
    elseif NL_solver == :NLPModels 
        # # # # non-linear programming (NLP)
    
        # # # # Nonlinear Least Squares (NLS) 
        tgfgfgfgg
    end
    ghhjjjfjffjk
    ADtype = :central    # A little more complicated but maybe obtaining no benefits.
end 

show_trace = false      # (=false, default) 
# ########## The initial solution noises
maxIterKing = 500       # (=200 default) The maximum inner iteration number for King's optimization
x_tol = epsT / 1000     # (=1e-13, default)
f_tol = epsT / 1000 
g_tol = epsT / 1000

 optimizer, ADtype
#  @show NL_solve,Optlibary,optimizer,ADtype, (is_C,is_constraint,is_bs,is_MTK)
 

