##### Single Bick system #####
struct BS_ODE_Sngl <: ODEType
    K::Real
    r::Real
    a2::Real
    a4::Real
    P::Vector{Real} #Pre-evluated parameters [K*cos(a4), 4r*cos(2*a2), -2*cos(a2)]

    function BS_ODE_Sngl(p::Dict{String, T}) where T<:Real
        BS_ODE_Sngl(p["K"], p["r"], p["a2"], p["a4"])
    end
    function BS_ODE_Sngl(K::Real, r::Real, a2::Real, a4::Real)
        P = [K * cos(a4), 4 * r * cos(2 * a2), -2 * cos(a2)]
        new(K, r, a2, a4, P)
    end
end

function Base.show(io::IO, bs::BS_ODE_Sngl)
    for property in propertynames(bs)
        value = getproperty(bs, property)
        println(io, "[$property::$(typeof(value))]: $(value)")
    end
end

function (bs::BS_ODE_Sngl)(dX, X, p = [], t = 0.0)
    A, B, C = bs.P
    dX[1] = sin(X[1]) * (A * (cos(X[3]) - cos(X[2])) + B * cos(X[1]) + C)
    dX[2] = sin(X[2]) * (A * (cos(X[1]) - cos(X[3])) + B * cos(X[2]) + C)
    dX[3] = sin(X[3]) * (A * (cos(X[2]) - cos(X[1])) + B * cos(X[3]) + C)
    return SVector{3}(dX)
end





##### Double coupled Bick system #####
struct BS_ODE_Duo <: ODEType
    Forward_ODE::BS_ODE_Sngl
    Backward_ODE::BS_ODE_Sngl
    
    Eps::Real
    Couple::Function

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
end

function Base.show(io::IO, bs::BS_ODE_Duo)
    println(io, "[Forward_ODE::$(typeof(bs.Forward_ODE))]:\n$(bs.Forward_ODE)")
    println(io, "[Backward_ODE::$(typeof(bs.Backward_ODE))]:\n$(bs.Backward_ODE)")
    println(io, "[Eps::$(typeof(bs.Eps))]: $(bs.Eps)")
    println(io, "[Couple::$(typeof(bs.Couple))]: $(bs.Couple)")
end

function (bs::BS_ODE_Duo)(dX, X, p = [], t = 0.0)
    dX[1:3] = bs.Forward_ODE(dX[1:3], X[1:3], [], t) + bs.Eps * [bs.Couple(p) for p in  X[4:6] - X[1:3]]
    dX[4:6] = bs.Backward_ODE(dX[4:6], X[4:6], [], t) + bs.Eps * [bs.Couple(p) for p in  X[1:3] - X[4:6]]
    return SVector{6}(dX)
end




##### Synchronization #####
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
end

function Base.show(io::IO, s::BS_Syncs)
    println(io, "Synchronization: $(["phi_$(sync[1]) => psi_$(sync[2])" for sync in s.syncs]) with delay: $(s.delay)")
end

function (s::BS_Syncs)(Phi::Vector)
    return [Phi[sync[2]] + s.delay for sync in s.syncs]
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
end

function Base.show(io::IO, bs::BS_ODE_Red)
    println(io, "[BS_ODE::$(typeof(bs.BS_ODE))]:\n$(bs.BS_ODE)")
    println(io, "[Syncs::$(typeof(bs.Syncs))]:\n$(bs.Syncs)")
    println(io, "[Eps::$(typeof(bs.Eps))]: $(bs.Eps)")
    println(io, "[Couple::$(typeof(bs.Couple))]: $(bs.Couple)")
end

function (bs::BS_ODE_Red)(dX, X, p = [], t = 0.0)
    Y = bs.Syncs(X)
    dX[1:3] = bs.BS_ODE(dX, X, [], t) + bs.Eps * [bs.Couple(p) for p in X - Y]
    return SVector{3}(dX)
end




##### Tangent Bick system #####
struct BS_ODE_Tan <: ODEType
    BS_ODE::BS_ODE_Duo
    Syncs::BS_Syncs

    function BS_ODE_Tan(p::Dict{String, T}, Couple::Union{String, Function}, Syncs::BS_Syncs) where T<:Real
        new(BS_ODE_Duo(p, Couple), Syncs)
    end
end

function Base.show(io::IO, bs::BS_ODE_Red)
    println(io, "[BS_ODE::$(typeof(bs.BS_ODE))]:\n$(bs.BS_ODE)")
end

function (bs::BS_ODE_Tan)(dW, W, X)
    Y = bs.Syncs(X)
    dW = jacobian((y, x) -> bs.BS_ODE(y, x, [], 0), dW, [X; Y]) * W
    return SVector{6}(dW)
end

















function Get_Tan_Red_Sys(dict::Dict, Couple_Func, syncWith)
    function BS_Tan_Red_Sys!(dW, W, X)
        phi1, phi2, phi3 = X
        psi1 = X[syncWith["phi1"]["psi"]] + syncWith["phi1"]["delay"]
        psi2 = X[syncWith["phi2"]["psi"]] + syncWith["phi2"]["delay"]
        psi3 = X[syncWith["phi3"]["psi"]] + syncWith["phi3"]["delay"]
        jacobian((y, x) -> Get_Uni_Sys(dict, Couple_Func)(y, x, [], 0), dW, [phi1, phi2, phi3, psi1, psi2, psi3]) * W
        #dW .= jacobian(ContinuousDynamicalSystem(Get_Uni_Sys(dict, Couple_Func), zeros(6), []), [phi1, phi2, phi3, psi1, psi2, psi3]) * W
        return SVector{6}(dW)
    end
    return BS_Tan_Red_Sys!
end

function Get_Tan_Red_Sys_2(dict::Dict, Couple_Func, syncWith)
    function BS_Tan_Red_Sys!(dW, W, X)
        jacobian((y, x) -> Get_Reduced_Sys(dict, Couple_Func, syncWith)(y, x, [], 0), dW, X) * W
        #dW .= jacobian(ContinuousDynamicalSystem(Get_Uni_Sys(dict, Couple_Func), zeros(6), []), X) * W
        return SVector{3}(dW)
    end
    return BS_Tan_Red_Sys!
end