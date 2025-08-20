## Mathematic constants

if isdefined(Main,:sqrt2)
    if is_change_datatype
        Pi = π |> datatype
        const sqrt2 = datatype(2)^0.5
        const sqrtpi = √(Pi)
        const pi3 = Pi^3             # = π^3)
        const sqrtpi3 = Pi^1.5       # = √(π^3) = π^(3/2)
        const sqrt2pi = √(2Pi)             # √(2π)
        const sqrt2invpi = (2 / Pi) ^ 0.5
        
        const lnpi = log(Pi)
        const lnsqrtpi3 = log(sqrtpi3)
    end
else
    Pi = π |> datatype
    const sqrt2 = datatype(2)^0.5
    const sqrtpi = √(Pi)
    const pi3 = Pi^3             # = π^3)
    const sqrtpi3 = Pi^1.5       # = √(π^3) = π^(3/2)
    const sqrt2pi = √(2Pi)             # √(2π)
    const sqrt2invpi = (2 / Pi) ^ 0.5
    
    const lnpi = log(Pi)
    const lnsqrtpi3 = log(sqrtpi3)
end

