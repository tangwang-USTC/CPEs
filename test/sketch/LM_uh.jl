
using Plots
using LaTeXStrings 
    # abs_dfLM = rel_dfLM = 1e-10
    uhvec1 = [0.44, 0.22, 0.1, 0.022, 1e-2,  2.2e-3, 1e-3, 2.2e-4, 1e-4, 2.2e-5]
    lvec1 =  [ 11,    9,   7,   5,     4,     3,     3,     2,      2,    2]

    # abs_dfLM = rel_dfLM = 1e-14
    uhvec3 = [3.0,2.0,1.45,1.414, 1.0, 0.965, 0.707,0.682, 0.482,  0.341,0.3, 0.145, 0.1, 0.048,    3e-2,    1.45e-2,   4.82e-3,   1.45e-3,   1e-3, 4.82e-4, 1e-4, 4.82e-5, 1e-5, 4.82e-6, 1e-7, 4.82e-8, 1e-8]
    lvec3 =  [45, 34,  28,  27,    25 , 23,    20,   20,    17,    14,    14,   10,   10,    8,       7,        6,          5,        4,        4,      4,      3,     3,      2,     2,      2,     2,     1]
    # K:     9.0  4.0   9.0  2.0       1.0     4.0   0.5   2.0      1.0     0.5  9e-2   9e-2    1e-2
    # û:      3.0  2.0  1.45  1.414     1.0   0.965 0.707  0.682   0.482    0.341  0.3   0.145   0.1
    # L_limit: 45  34    28    27        23      25   20     20      17      14    14     10     10
    # index_L: 37  28    23  (21,-13)  (17,-13)  17   15   (13,-13) (11,-13) 10  (9,-12) (3,-12) (2~3,-12)
    # ℓ>(,)[1]: -2  -2   -2  (18,-3)   (18,-3)   -3   -3   (14,-3)  (12,-3)  -3    -3    -3~-6   (≥4,-3~-7)
    # index_uL:-13 -13   -13   -13      -13     -13   -13    -13     -13     -13   -13  (11,-13) (8+,-13)    ✓✓✓
    #                                                                                            (L≥10,1e-4)
    # K:        1e-2                 9e-4       1e-4     9e-4      9e-6    1e-4      9e-6     1e-6
    # û:         0.1       0.048      3e-2       1e-2    1.45e-2    3e-3   4.82e-3   1.45e-3   1e-3
    # L_limit:   10          8         7           6       6         5        5        4        4
    # index_L: (2~3,-12)  (1+,-12)   (0,-3)     (0,-3)   (0,-3)    (0,-5)   (0,-5)   (0,-7)  (ℓ=0,-7)(≥1,-13)  ✓✓
    #         (≥4,-3~-7) (≥3,-3~-7) (≥1,-5~-8) (-7~-10)  (-6~-10) (-9~-10) (-9~-11) (-10~-12) (≥1,-13)
    # index_uL: (8+,-13)   (2,-12)    (0,-3)    (0,-4)   (0,-3)   (0,-5)   (0,-5)   (0,-7)    -8            ✓✓
    #          (L≥10,1e-4) (3,1e-4)  (1,1e-4)   (-4~-5)  (-4~-5)  (-5~-7)  (-5~-6)  (-7~-8)
    # fMvp2,(ℓ=0):-1(abs-6) -2(-8)   -3(-8)    -5(-10)   -5(-10)  -7(-12)  -6(-12)  -9(-14)  -9(-14)
    # fMvp3,(ℓ=0):-1(abs-6) -2(-8)   -3(-9)    -5(-10)   -5(-10)  -7(-13)  -6(-12)  -9(-14)  -9(-14)
    #
    # K:        1e-6            1e-8         1e-10        1e-14        1e-16
    # û:         1e-3    4.82e-4 1e-4 4.82e-5 1e-5 4.82e-6 1e-7 4.82e-8 1e-8  4.82e-9
    # L_limit:    4        4      3     3      2     2      2     2     1     1
    # index_L: (≥1,-13) (≥1,-13) -13   -13    -13   -13    -13   -13   -13   -13     ✓✓✓
    #          (ℓ=0,-7) (ℓ=0,-9)
    # index_uL:  -8        -9    -11   -13    -13   -13    -13   -13   -13   -13     ✓✓
    # fMvp2,(ℓ=0): -9     -10    -13   -13    -14   -14    -14   -14   -14   -14
    # fMvp3,(ℓ=0): -9     -10    -13   -13    -14   -14    -14   -14   -14   -14
    
    # LM
    
    label = string(L"Atol_{df}=",1e-10)
    xlabel = L"\hat{u}_a"
    ylabel = L"L_M"
    xscale = :log10
    xscale = :norm
    pLMuh = plot(uhvec1,lvec1, ylabel=ylabel, label=label, line=(3,:solid),
                xlabel=xlabel,xscale=xscale)
    if xscale == :log10
        savefig(string("LM_uh1_log.pdf"))
    else
        savefig(string("LM_uh1.pdf"))
    end
    display(pLMuh)
    
    label = string(L"Atol_{df}=",1e-14)
    pLMuh = plot(uhvec3,lvec3, ylabel=ylabel, label=label, line=(3,:solid),
                xlabel=xlabel,xscale=xscale)
    if xscale == :log10
        savefig(string("LM_uh3_log.pdf"))
    else
        savefig(string("LM_uh3.pdf"))
    end
    display(pLMuh)

    # The relative error of fitting will be: `10^(index_)`, (ℓ→0,ℓₛ~LM1)

    # Interval `û ∈ [0.3,1e-3]` need special method to deciding the fitting process.
    # When `û ≤ 1e-3`, procedure `weightfunctions_v` is the best one.
    # When `û ≥ 0.3`, procedure `weightfunctions_uL` is the best one.
    # But when there are `û[1] ≤ 1e-3` and `û[2] ≥ 0.3` at the same time, how to decide the fitting process?

    