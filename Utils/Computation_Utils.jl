struct py_sol{T<:Any}
    t::Vector{Float64}
    y::Vector{Vector{Float64}}
    t_events::Vector{Vector{Float64}}
    y_events::Vector{Matrix{Float64}}
    retcode::Union{String, Nothing}

    function py_sol{T}(sol::ODESolution; event_arr::Union{MySavedValues, Nothing} = nothing) where T
        py_sol_t = copy(sol.t)
        py_sol_y = Vector{Vector{Float64}}(undef, size(sol)[1])
        for i = 1 : size(sol)[1]
            py_sol_y[i] = Vector{Float64}(undef, size(sol)[2])
            for j = 1 : size(sol)[2]
                py_sol_y[i][j] = sol.u[j][i]
            end
        end

        if (event_arr !== nothing) && (!isempty(event_arr.t))
            py_sol_t_events = copy(event_arr.t)
            py_sol_y_events = Matrix{Float64}(undef, length(event_arr.saveval), length(event_arr.saveval[1]))
            for i = 1 : size(py_sol_y_events)[1]
                for j = 1 : size(py_sol_y_events)[2]
                    py_sol_y_events[i, j] = event_arr.saveval[i][j]
                end
            end
        else
            py_sol_t_events = Vector{Float64}(undef, 0)
            py_sol_y_events = Matrix{Float64}(undef, 0, 0)
        end

        py_sol_retcode = string(sol.retcode)

        new{T}(py_sol_t, py_sol_y, [py_sol_t_events], [py_sol_y_events], py_sol_retcode)
    end
end

function Base.show(io::IO, sol::py_sol)
    println(io, "[t]: $(length(sol.t))-element $(typeof(sol.t))")
    println(io, "[y]: $(length(sol.y))-element $(typeof(sol.y)) with size $(length(sol.y[1]))")
    println(io, "[t_events]: $(length(sol.t_events))-element $(typeof(sol.t_events)) with size $(length(sol.t_events[1]))")
    println(io, "[y_events]: $(length(sol.y_events))-element $(typeof(sol.y_events)) with size $(size(sol.y_events[1])[1])x$(size(sol.y_events[1])[2])")
    println(io, "[retcode]: $(sol.retcode)")
end

function SolveODE(ODESystem, x_0, time_span; p = [], alg = nothing, trans_time = 0.0, callback = nothing, kwargs = Dict())

    # Transit time process integration
    prob_julia = ODEProblem{true, SciMLBase.FullSpecialize}(ODESystem, x_0, (0.0, trans_time), p)
    sol_julia = solve(prob_julia, alg; save_everystep = false, kwargs...)

    # Define callback function
    cb, observer_array = callback === nothing ? (nothing, nothing) : Get_callback(callback, length(x_0))

    # Main time peocess integration
    prob_julia = ODEProblem{true, SciMLBase.FullSpecialize}(ODESystem, sol_julia.u[end], (trans_time, time_span + trans_time), p)
    sol_julia = solve(prob_julia, alg; callback = cb, kwargs...)

    # Solution convertation into python style
    return py_sol{typeof(ODESystem)}(sol_julia; event_arr = observer_array)
end








mutable struct TangentODE <: ODEType
    ODE::Union{Function, ODEType}
    Jacobian::Union{Function, ODEType}
    ODEDim::Int
    JacDim::Int

    function TangentODE(ODE::Union{Function, ODEType}, ODEDim::Int)
        Jacobian = z -> jacobian((y, x) -> ODE(y, x, [], 0), zeros(ODEDim), z)
        new(ODE, Jacobian, ODEDim, ODEDim)
    end
    function TangentODE(ODE::Union{Function, ODEType}, custom_jacobian::Union{Function, ODEType}, ODEDim::Int, JacDim::Int)
        Jacobian = z -> jacobian((y, x) -> custom_jacobian(y, x, [], 0), zeros(ODEDim), z)
        new(ODE, Jacobian, ODEDim, JacDim)
    end
end

function Base.show(io::IO, ode::TangentODE)
    for property in propertynames(ode)
        value = getproperty(ode, property)
        println(io, "[$property::$(typeof(value))]: $(value)")
    end
end

function (ode::TangentODE)(dX, X, p = [], t = 0.0)
    dX[1 : ode.ODEDim] .= ode.ODE(dX[1 : ode.ODEDim], X[1 : ode.ODEDim], p, t)

    Jac = ode.Jacobian(X[1 : ode.ODEDim])
    for i in 0 : ode.JacDim - 1
        dX[ode.ODEDim + 1 + ode.JacDim * i : ode.ODEDim + ode.JacDim * (i + 1)] .=
                Jac * X[ode.ODEDim + 1 + ode.JacDim * i : ode.ODEDim + ode.JacDim * (i + 1)]
    end

    return SVector{ode.ODEDim + ode.JacDim^2}(dX)
end