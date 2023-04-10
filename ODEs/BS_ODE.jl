##### Single Bick system #####
struct BS_ODE_Sngl <: ODEType
    K::Float64
    r::Float64
    a2::Float64
    a4::Float64
    P::SVector{3, Float64} #Pre-evluated parameters [K*cos(a4), 4r*cos(2*a2), -2*cos(a2)]

    function BS_ODE_Sngl(p::Dict{String, T}) where T<:Real
        BS_ODE_Sngl(p["K"], p["r"], p["a2"], p["a4"])
    end
    function BS_ODE_Sngl(K::Real, r::Real, a2::Real, a4::Real)
        P = SA[K * cos(a4), 4.0 * r * cos(2.0 * a2), -2.0 * cos(a2)]
        new(K, r, a2, a4, P)
    end

    function Base.show(io::IO, bs::BS_ODE_Sngl)
        for property in propertynames(bs)
            value = getproperty(bs, property)
            println(io, "[$property::$(typeof(value))]: $(value)")
        end
    end
end

@inbounds function (bs::BS_ODE_Sngl)(X, p = [], t = 0.0)
    phi1, phi2, phi3 = X
    A, B, C = bs.P

    dX1 = sin(phi1) * (A * (cos(phi3) - cos(phi2)) + B * cos(phi1) + C)
    dX2 = sin(phi2) * (A * (cos(phi1) - cos(phi3)) + B * cos(phi2) + C)
    dX3 = sin(phi3) * (A * (cos(phi2) - cos(phi1)) + B * cos(phi3) + C)
    return SA[dX1, dX2, dX3]
end





##### Couple functions for double Bick system #####
@inline Cos_Couple(x) = 1.0 - cos(x)
Base.print(io::IO, ::typeof(Cos_Couple)) = print(io, "[Cos_Couple]: 1 - cos(x)")


@inline Exp_Couple(x) = 1.0 / (1.0 + exp(10.0 * cos(x)))
Base.print(io::IO, ::typeof(Exp_Couple)) = print(io, "[Exp_Couple]: 1.0 / (1.0 + exp(10.0 * cos(x)))")





##### Double coupled Bick system #####
struct BS_ODE_Duo <: ODEType
    Forward_ODE::BS_ODE_Sngl
    Backward_ODE::BS_ODE_Sngl
    
    Eps::Float64
    Couple::Union{Function, typeof(Cos_Couple), typeof(Exp_Couple)}

    function BS_ODE_Duo(p::Dict{String, T}, Couple::String) where T<:Real
        BS_ODE_Duo(p, AnonFunc(Couple))
    end
    function BS_ODE_Duo(p::Dict{String, T}, Couple::Function) where T<:Real
        BS_ODE_Duo(p["K"], p["r"], p["a2"], p["a4"], p["Eps"], Couple)
    end
    function BS_ODE_Duo(K::Real, r::Real, a2::Real, a4::Real, Eps::Real, Couple::Function)
        Fwd_ODE = BS_ODE_Sngl(K, r, a2, a4)
        Bckwd_ODE = BS_ODE_Sngl(-K, r, a2, a4)

        new(Fwd_ODE, Bckwd_ODE, Eps, Couple)
    end

    function Base.show(io::IO, bs::BS_ODE_Duo)
        println(io, "[Forward_ODE::$(typeof(bs.Forward_ODE))]:\n$(bs.Forward_ODE)")
        println(io, "[Backward_ODE::$(typeof(bs.Backward_ODE))]:\n$(bs.Backward_ODE)")
        println(io, "[Eps::$(typeof(bs.Eps))]: $(bs.Eps)")
        println(io, "[Couple::$(typeof(bs.Couple))]: $(bs.Couple)")
    end
end

@inbounds function (bs::BS_ODE_Duo)(X, p = [], t = 0.0)
    Phi = SVector{3, Float64}(X[1:3])
    Psi = SVector{3, Float64}(X[4:6])

    return SVector{6}([bs.Forward_ODE(Phi); bs.Backward_ODE(Psi)] + bs.Eps * bs.Couple.([Psi - Phi; Phi - Psi]))
end





##### Synchronization #####
function GetSyncs(sol::py_sol)
    return BS_Syncs(sol)
end

struct BS_Syncs
    syncs::Vector{Pair}
    delay::Real

    function BS_Syncs(sol::py_sol)
        synphase_std = [[std(map(cos, Psi) - map(cos, Phi)) for Psi in sol.y[4:6]] for Phi in sol.y[1:3]]
        antiphase_std = [[std(map(cos, Psi) - map(x -> cos(x + pi), Phi)) for Psi in sol.y[4:6]] for Phi in sol.y[1:3]]
        if minimum(synphase_std[1]) < minimum(antiphase_std[1])
            delay = 0.0
            syncs = [Pair(i, argmin(s)) for (i, s) in enumerate(synphase_std)]
        else 
            delay = pi
            syncs = [Pair(i, argmin(s)) for (i, s) in enumerate(antiphase_std)]
        end
        new(syncs, delay)
    end
    function BS_Syncs(syncs::Vector{Pair}, delay::Real)
        new(syncs, delay)
    end

    function Base.show(io::IO, s::BS_Syncs)
        println(io, "Synchronization: $(["phi_$(sync[1]) => psi_$(sync[2])" for sync in s.syncs]) with delay: $(s.delay)")
    end
    
    function (s::BS_Syncs)(Phi)
        return [Phi[sync[2]] + s.delay for sync in s.syncs]
    end
end






##### Reduced Bick system #####
struct BS_ODE_Red <: ODEType
    BS_ODE::BS_ODE_Sngl
    Couple::Function
    Syncs::BS_Syncs
    Eps::Real

    function BS_ODE_Red(p::Dict{String, T}, Couple::String, Syncs::BS_Syncs)  where T<:Real
        BS_ODE_Red(p["K"], p["r"], p["a2"], p["a4"], p["Eps"], AnonFunc(Couple), Syncs)
    end
    function BS_ODE_Red(p::Dict{String, T}, Couple::Function, Syncs::BS_Syncs)  where T<:Real
        BS_ODE_Red(p["K"], p["r"], p["a2"], p["a4"], p["Eps"], Couple, Syncs)
    end
    function BS_ODE_Red(K::T, r::T, a2::T, a4::T, Eps::T, Couple::Function, Syncs::BS_Syncs)  where T<:Real
        new(BS_ODE_Sngl(K, r, a2, a4), Couple, Syncs, Eps)
    end

    function Base.show(io::IO, bs::BS_ODE_Red)
        println(io, "[BS_ODE::$(typeof(bs.BS_ODE))]:\n$(bs.BS_ODE)")
        println(io, "[Syncs::$(typeof(bs.Syncs))]:\n$(bs.Syncs)")
        println(io, "[Eps::$(typeof(bs.Eps))]: $(bs.Eps)")
        println(io, "[Couple::$(typeof(bs.Couple))]: $(bs.Couple)")
    end
end

function (bs::BS_ODE_Red)(X, p = [], t = 0.0)
    Y = bs.Syncs(X)
    dX = bs.BS_ODE(X, p, t) + bs.Eps * [bs.Couple(p) for p in X - Y]
    return SVector{3}(dX)
end