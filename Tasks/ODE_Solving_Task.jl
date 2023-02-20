mutable struct ODESolve_Task{T<:ODEType} <: AbstractTask
    ## Required fileds ##
    ODE::T # Function with signature f(dX, u, p, t) or f(u, p, t)
    time_span::Float64
    init_cond::Vector{Float64}

    ## Optional fields ##
    trans_time::Union{Float64, Nothing}
    alg::Union{Any, Nothing}
    callback::Union{Function, String, Expr, Nothing}
    kwargs::Union{Dict{Symbol, Any}, Nothing}

    solution::Union{py_sol, Nothing}

    function ODESolve_Task(ODE::T, time_span, init_cond;
                            trans_time = 0.0,
                            alg = nothing,
                            callback = nothing,
                            kwargs = Dict()) where T <: ODEType
        new{T}(ODE, time_span, init_cond, trans_time, alg, callback, kwargs, nothing)
    end
end

function (Task::ODESolve_Task)()
    sol = SolveODE(Task.ODE, Task.init_cond, Task.time_span;
                    alg = Task.alg,
                    trans_time = Task.trans_time,
                    kwargs = Task.kwargs,
                    callback = Task.callback);
    Task.solution = sol
    return sol
end

function Base.show(io::IO, t::ODESolve_Task{T}) where T
    
    for property in propertynames(t)
        value = getproperty(t, property)
        println(io, "[$property::$(typeof(value))]: $(value)")
    end
end








mutable struct TanODESolve_Task <: AbstractTask
    ## Required fileds ##
    ODE::TangentODE
    time_span::Float64
    u0::Vector{Float64}
    Q0::Matrix{Float64}

    ## Optional fields ##
    alg::Union{Any, Nothing}
    kwargs::Union{Dict{Symbol, Any}, Nothing}

    solution::Union{py_sol, Nothing}

    function TanODESolve_Task(ODE::TangentODE, time_span::Float64, u0::Vector{Float64}, Q0::Matrix{Float64};
                            alg = nothing,
                            kwargs = Dict())
        new(ODE, time_span, u0, Q0, alg, kwargs, nothing)
    end
end

function (Task::TanODESolve_Task)()
    ic = cat(u0, vec(Q0), dims = 1)

    sol = SolveODE(Task.ODE, ic, Task.time_span;
                    alg = Task.alg,
                    kwargs = Task.kwargs);
    Task.solution = sol
    return sol
end

function Base.show(io::IO, t::TanODESolve_Task)
    for property in propertynames(t)
        value = getproperty(t, property)
        println(io, "[$property::$(typeof(value))]: $(value)")
    end
end