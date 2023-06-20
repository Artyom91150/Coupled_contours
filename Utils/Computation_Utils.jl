    ## Structure of solution like solve_ivp from SyPy ##
struct py_sol{T<:Any}
    t::Vector{Float64}
    y::Vector{Vector{Float64}}
    t_events::Vector{Vector{Float64}}
    y_events::Vector{Matrix{Float64}}
    retcode::Union{String, Nothing}

    function py_sol{T}(sol::ODESolution; event_arr::Union{MySavedValues, Nothing} = nothing) where T
        py_sol_t = copy(sol.t)
        py_sol_y = [[u[i] for u in sol.u] for i in range(1, length(sol.u[1]))]

        if (event_arr !== nothing) && (!isempty(event_arr.t))
            py_sol_t_events = copy(event_arr.t)
            py_sol_y_events = reduce(vcat, transpose.([[v for v in tup] for tup in event_arr.saveval]))
        else
            py_sol_t_events = Vector{Float64}(undef, 0)
            py_sol_y_events = Matrix{Float64}(undef, 0, 0)
        end

        py_sol{T}(py_sol_t, py_sol_y, [py_sol_t_events], [py_sol_y_events], string(sol.retcode))
    end
    function py_sol{T}(py_sol_t::Vector{S}, py_sol_y::Vector{Vector{S}},
                        py_sol_t_events::Vector{Vector{S}} = [Vector{S}(undef, 0)],
                        py_sol_y_events::Vector{Matrix{S}} = [Matrix{S}(undef, 0, 0)],
                        py_sol_retcode::Union{String, Nothing} = nothing) where {T <: Union{ODEType, Any}, S <: Float64}
        new{T}(py_sol_t, py_sol_y, py_sol_t_events, py_sol_y_events, py_sol_retcode)
    end

    function Base.show(io::IO, sol::py_sol)
        println(io, "[t]: $(length(sol.t))-element $(typeof(sol.t))")
        println(io, "[y]: $(length(sol.y))-element $(typeof(sol.y)) with size $(length(sol.y[1]))")
        println(io, "[t_events]: $(length(sol.t_events))-element $(typeof(sol.t_events)) with size $(length(sol.t_events[1]))")
        println(io, "[y_events]: $(length(sol.y_events))-element $(typeof(sol.y_events)) with size $(size(sol.y_events[1])[1])x$(size(sol.y_events[1])[2])")
        println(io, "[retcode]: $(sol.retcode)")
    end
end





    ## Result of ODE solution task ##
mutable struct Trajectory{T, P}
    ## Trajectory info ##
    ODE::T # Function with signature f(dX, u, p, t) or f(u, p, t)
    solution::py_sol
    p::P

    ## Numerical info ##
    alg::Union{Any, Nothing}
    kwargs::Union{Dict{Symbol, Any}, Nothing}
    callback::Union{Function, String, Expr, Nothing, Any}

    function Trajectory(ODE::T, sol::py_sol, p::P;
                        alg = nothing,
                        callback = nothing,
                        kwargs = nothing) where {T <: Union{ODEType, Any}, P <: AbstractArray}
        new{T, P}(ODE, sol, p, alg, kwargs, callback)
    end

    function Base.show(io::IO, t::Trajectory{T}) where T
        for property in propertynames(t)
            value = getproperty(t, property)
            println(io, "[$property::$(typeof(value))]: $(value)")
        end
    end

    function (trj::Trajectory)(t)
        kwargs = Dict(:save_everystep => false)
        if trj.kwargs !== nothing merge(trj.kwargs, Dict(:save_everystep => false)) end
        
        nearest_point_idx = BinarySearchLeftBorder(trj.solution.t, t)
        t0 = trj.solution.t[nearest_point_idx]
        u0 = [y[nearest_point_idx] for y in trj.solution.y]
        
        new_trj = SolveODE(trj.ODE, u0, abs(t - t0); p = trj.p, alg = trj.alg, init_time = t0, kwargs = kwargs)
        return [y[end] for y in new_trj.solution.y]
    end
end

function BinarySearchLeftBorder(arr, key)
    left = 2
    right = length(arr)
    while (left <= right)
        mid = round(Int, (right + left) / 2)
        if (arr[mid] == key)
            return mid
        elseif (arr[mid] < key)
            left = mid + 1
        else
            right = mid - 1
        end
    end
    return left - 1
end





    ## ODE solving method ##
function SolveODE(ODE, x_0, time_span; p = [], alg = nothing, trans_time = 0.0, init_time = 0.0, callback = nothing, kwargs = Dict())

    # Transit time process integration
    transit_time = (init_time, init_time + trans_time)
    #x_0 = typeof(x_0) <: SVector ? x_0 : SVector{length(x_0)}(x_0)

    prob_julia = ODEProblem{false, SciMLBase.FullSpecialize}(ODE, x_0, transit_time, p)
    sol_julia = solve(prob_julia, alg; save_everystep = false, kwargs...)

    # Define callback function
    cb, observer_array = callback === nothing ? (nothing, nothing) : Get_callback(callback, length(x_0))

    # Main time process integration
    integration_time = transit_time .+ (trans_time, time_span)
    prob_julia = ODEProblem{false, SciMLBase.FullSpecialize}(ODE, sol_julia.u[end], integration_time, p)
    sol_julia = solve(prob_julia, alg; callback = cb, kwargs...)

    # Solution convertation into python style
    pysol = py_sol{typeof(ODE)}(sol_julia; event_arr = observer_array)

    return Trajectory(ODE, pysol, p, alg = alg, callback = callback, kwargs = kwargs)
end






mutable struct TangentODE <: ODEType
    ODE::Union{Function, ODEType}
    Jacobian::Union{Function, ODEType}
    ODEDim::Int
    JacDim::Int

    function TangentODE(ODE::Union{Function, ODEType}, ODEDim::Int)
        Jacobian = z -> jacobian(x -> ODE(x, [], 0.0), z)
        new(ODE, Jacobian, ODEDim, ODEDim)
    end
    function TangentODE(ODE::Union{Function, ODEType}, custom_jacobian::Union{Function, ODEType}, ODEDim::Int, JacDim::Int)
        new(ODE, custom_jacobian, ODEDim, JacDim)
    end
end

function Base.show(io::IO, ode::TangentODE)
    for property in propertynames(ode)
        value = getproperty(ode, property)
        println(io, "[$property::$(typeof(value))]: $(value)")
    end
end

function (ode::TangentODE)(X, p = [], t = 0.0)
    dX = ode.ODE(X[1 : ode.ODEDim], p, t)

    Jac = ode.Jacobian(X[1 : ode.ODEDim])
    dW = zeros(ode.JacDim^2)

    for i in 0 : ode.JacDim - 1
        dW[1 + ode.JacDim * i : ode.JacDim * (i + 1)] .=
                Jac * X[ode.ODEDim + 1 + ode.JacDim * i : ode.ODEDim + ode.JacDim * (i + 1)]
    end

    return SVector{ode.ODEDim + ode.JacDim^2}([dX; dW])
end

function splitSol(sol::py_sol{T}) where T <: TangentODE
    ode_dim, tan_dim = getTanODEDims(length(sol.y))

    ode_sol_y = sol.y[1 : ode_dim]
    tan_sol_y = sol.y[ode_dim + 1 : end]

    ode_sol_events = !isempty(sol.y_events[1]) ? sol.y_events[1][1 : ode_dim] : sol.y_events
    tan_sol_events = !isempty(sol.y_events[1]) ? sol.y_events[1][ode_dim + 1 : end] : sol.y_events

    ode_sol = py_sol{T}(sol.t, ode_sol_y, sol.t_events, ode_sol_events, sol.retcode)
    tan_sol = py_sol{T}(sol.t, tan_sol_y, sol.t_events, tan_sol_events, sol.retcode)

    return ode_sol, tan_sol
end
function getTanODEDims(sol_len::Int)
    tan_dim = round(Int, sqrt(sol_len), RoundDown)^2
    if (tan_dim == sol_len) tan_dim -= 1 end
    ode_dim = sol_len - tan_dim
    return ode_dim, tan_dim
end
function getTanEigens(tan_sol::py_sol)
    tan_size = round(Int, sqrt(length(tan_sol.y)))
    if tan_size^2 != length(tan_sol.y) throw(DimensionMismatch) end

    eigens = [eigen(reshape([t[i] for t in tan_sol.y], (tan_size, tan_size))).values for i in range(1, length(tan_sol.t))]
    eigens = [[e[i] for e in eigens] for i in range(1, length(eigens[1]))]
    return eigens
end








function GetLEs(ds, u0, p, Eps_mesh, time_span, trans_time, diffeq)
    lvDs = CoupledODEs(ds, u0, p; diffeq)
    λs = Matrix{Float64}(undef, dimension(lvDs), length(Eps_mesh))

    @showprogress "Computing..." for (i, eps) in enumerate(Eps_mesh)
        set_parameter!(lvDs, 4, eps)
        set_state!(lvDs, u0)
        tanDs = TangentDynamicalSystem(lvDs)

        λ = lyapunovspectrum(tanDs, time_span; Ttr = trans_time)
        if all(isapprox.(ds(get_state(tanDs), current_parameters(lvDs)), 0.0; atol = 1e-8, rtol = 1e-8))
            #@warn "Trajectory get into equilibria point"
            λs[:, i] = [NaN for _ in 1:dimension(lvDs)]
        else
            λs[:, i] = λ
        end
    end
    return λs
end