struct SL_ODE <: ODEType
    r::Real
    omega::Real
    s::Real

    function SL_ODE(p::Dict)
        SL_ODE(p["r"], p["omega"], p["s"])
    end
    function SL_ODE(r::Real, omega::Real, s::Real)
        new(r, omega, s)
    end
end

function Base.show(io::IO, sl::SL_ODE)
    for property in propertynames(sl)
        value = getproperty(sl, property)
        println(io, "[$property::$(typeof(value))]: $(value)")
    end
end

function (sl::SL_ODE)(dX, X, p = [], t = 0.0)
    x, y = X
    dX[1] = sl.omega * y + sl.s * x * (sl.r * sl.r - x * x - y * y)
    dX[2] = -sl.omega * x + sl.s * y * (sl.r * sl.r - x * x - y * y)
    return SVector{2}(dX)
end