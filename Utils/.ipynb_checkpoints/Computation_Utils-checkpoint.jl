module Computation_Utils
using DifferentialEquations
using DynamicalSystems
include("Callback_Utils.jl")

struct py_sol
    t::Vector{Float64}
    y::Vector{Vector{Float64}}
    t_events::Vector{Vector{Float64}}
    y_events::Vector{Matrix{Float64}}
end   

function jl2pySol(sol, obs_arr = [])
    py_sol_y = Vector{Vector{Float64}}(undef, size(sol)[1])
    for i = 1 : size(sol)[1]
        py_sol_y[i] = Vector{Float64}(undef, size(sol)[2])
        for j = 1 : size(sol)[2]
            py_sol_y[i][j] = sol.u[j][i]
        end
    end
    
    if obs_arr != [] && obs_arr.t != []
        py_sol_y_events = Matrix{Float64}(undef, length(obs_arr.saveval), length(obs_arr.saveval[1]))
        for i = 1 : size(py_sol_y_events)[1]
            for j = 1 : size(py_sol_y_events)[2]
                py_sol_y_events[i, j] = obs_arr.saveval[i][j]
            end
        end
        return py_sol(sol.t, py_sol_y, [obs_arr.t], [py_sol_y_events])
    else
        return py_sol(sol.t, py_sol_y, [], [])
    end
end

function SolveODE(ODESystem, x_0, time_span; trans_time = 0.0, callback_func = nothing, kwargs = Dict())
    # Transit time process integration
    prob_julia = ODEProblem(ODESystem, x_0, (0.0, trans_time))
    sol_julia = solve(prob_julia; kwargs = merge(kwargs, Dict("save_everystep" => false)))

    # Define callback function
    cb, observer_array = callback_func === nothing ? (nothing, []) : Callback_Utils.GetCallBack(callback_func, length(x_0))

    # Main time peocess integration
    prob_julia = ODEProblem(ODESystem, sol_julia.u[end], (0.0, time_span))
    sol_julia = solve(prob_julia, callback = cb; kwargs = kwargs)

    # Solution convertation into python style
    return jl2pySol(sol_julia, observer_array)
end

end