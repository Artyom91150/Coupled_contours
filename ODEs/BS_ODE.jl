##### Single Bick system #####
struct BS_ODE_Sngl
    K::Float64
    r::Float64
    a2::Float64
    a4::Float64
    P::SVector{3, Float64} #Pre-evluated parameters [K*cos(a4), 4r*cos(2*a2), -2*cos(a2)]

    function BS_ODE_Sngl(p::Dict{String, T}) where T <: Real
        BS_ODE_Sngl(p["K"], p["r"], p["a2"], p["a4"])
    end
    function BS_ODE_Sngl(K::Float64, r::Float64, a2::Float64, a4::Float64)
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

@inbounds function (bs::BS_ODE_Sngl)(X::T, p, t = 0.0) where T
    phi1, phi2, phi3 = X
    A, B, C = bs.P

    dX1 = sin(phi1) * (A * (cos(phi3) - cos(phi2)) + B * cos(phi1) + C)
    dX2 = sin(phi2) * (A * (cos(phi1) - cos(phi3)) + B * cos(phi2) + C)
    dX3 = sin(phi3) * (A * (cos(phi2) - cos(phi1)) + B * cos(phi3) + C)
    return T([dX1, dX2, dX3])
end





##### Couple functions for double Bick system #####
Dummy_Couple(x::String) = "undef couple"

Cos_Couple(x::Float64) = (1.0 - cos(x))/2.0
Cos_Couple(x::String) = "\\frac{1 - cos($x)}{2}"

Sin_Couple(x::Float64) = sin(x)/2.0
Sin_Couple(x::String) = "sin($x)"

Exp_Couple(x, k = 10.0, zero = 0.0) = 1.0 / (1.0 + exp(k * cos(x))) - zero
Exp_Couple(x::String) = "\\frac{1}{1 + e^{k cos($x)}}"

function Get_Exp_Couple(k = 10.0, zero = nothing)
    zero = isnothing(zero) ? Exp_Couple(0.0, k) : zero
    function f(x)
        return Exp_Couple(x, k, zero)
    end
    return f
end




##### Double coupled Bick system #####
struct BS_ODE_Duo
    Forward_ODE::BS_ODE_Sngl
    Backward_ODE::BS_ODE_Sngl
    
    Eps::Float64
    Couple::Union{Function, typeof(Cos_Couple), typeof(Exp_Couple)}

    function BS_ODE_Duo(p::Dict{String, T}, Couple::String) where T <: Real
        BS_ODE_Duo(p, AnonFunc(Couple))
    end
    function BS_ODE_Duo(p::Dict{String, T}, Couple::Function) where T <: Real
        BS_ODE_Duo(p["K"], p["r"], p["a2"], p["a4"], p["Eps"], Couple)
    end
    function BS_ODE_Duo(K::Float64, r::Float64, a2::Float64, a4::Float64, Eps::Float64, Couple::Function)
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

@inbounds function (bs::BS_ODE_Duo)(X::T, p = [], t = 0.0) where T
    Phi = X[1:3]
    Psi = X[4:6]

    return SVector{6}([bs.Forward_ODE(Phi); bs.Backward_ODE(Psi)] + bs.Eps * bs.Couple.([Psi - Phi; Phi - Psi]))
end

function Get_Fast_BS(Couple)
    @inbounds function BS_ODE_Duo_Fast(X::T, p, t = 0.0) where T
        phi1, phi2, phi3, psi1, psi2, psi3 = X
        A, B, C, Eps = p
    
        dPhi1 = sin(phi1) * (A * (cos(phi3) - cos(phi2)) + B * cos(phi1) + C) + Eps * Couple(psi1 - phi1)
        dPhi2 = sin(phi2) * (A * (cos(phi1) - cos(phi3)) + B * cos(phi2) + C) + Eps * Couple(psi2 - phi2)
        dPhi3 = sin(phi3) * (A * (cos(phi2) - cos(phi1)) + B * cos(phi3) + C) + Eps * Couple(psi3 - phi3)
    
        dPsi1 = sin(psi1) * (-A * (cos(psi3) - cos(psi2)) + B * cos(psi1) + C) + Eps * Couple(phi1 - psi1)
        dPsi2 = sin(psi2) * (-A * (cos(psi1) - cos(psi3)) + B * cos(psi2) + C) + Eps * Couple(phi2 - psi2)
        dPsi3 = sin(psi3) * (-A * (cos(psi2) - cos(psi1)) + B * cos(psi3) + C) + Eps * Couple(phi3 - psi3)
        return T([dPhi1, dPhi2, dPhi3, dPsi1, dPsi2, dPsi3])
    end
    return BS_ODE_Duo_Fast
end