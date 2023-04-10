using PyCall

pushfirst!(PyVector(pyimport("sys")."path"), "")
MU = pyimport("Utils.Graphics_Utils")

struct TimeSeries
    projFunc::Union{Any, Nothing}
    varNames::Union{Vector, Nothing}
    title::Union{String, Nothing}
    plotKwargs::Union{Dict, Nothing}

    function TimeSeries(properties::Dict = Dict())
        projFunc = get(properties, "projFunc", "projNone")
        varNames = get(properties, "varNames", nothing)
        title = get(properties, "title", nothing)
        plotKwargs = get(properties, "plotKwargs", Dict())

        TimeSeries(; projFunc, varNames, title, plotKwargs)
    end
    function TimeSeries(; projFunc = nothing, varNames = nothing, title = nothing, plotKwargs = Dict())
        projFunc = projFunc === nothing ? MU["projNone"] : MU[projFunc]
        new(projFunc, varNames, title, plotKwargs)
    end

    function (p::TimeSeries)(sol::py_sol; fig = nothing, savePath = nothing)
        if p.varNames !== nothing && length(p.varNames) < length(sol.y) throw(DimensionMismatch("varNames ($(length(p.varNames)) mismatch sol.y ($(length(sol.y))))")) end
        MU.plotTimeSeries(sol, fig = fig, savePath = savePath, projFunc = p.projFunc, varNames = p.varNames, title = p.title, plotKwargs = p.plotKwargs)
    end
end




struct Poincare
    projFunc::Union{Any, Nothing}
    varNames::Union{Vector, Nothing}
    varPairs::Union{Vector, Nothing}
    title::Union{String, Nothing}
    plotKwargs::Union{Dict, Nothing}
    showEvents::Union{Bool, Nothing}

    function Poincare(properties::Dict = Dict())
        projFunc = get(properties, "projFunc", "projNone")
        varNames = get(properties, "varNames", nothing)
        varPairs = get(properties, "varPairs", nothing)
        title = get(properties, "title", nothing)
        plotKwargs = get(properties, "plotKwargs", Dict())
        showEvents = get(properties, "showEvents", nothing)

        Poincare(; varPairs, projFunc, varNames, title, plotKwargs, showEvents)
    end
    function Poincare(; varPairs = nothing, projFunc = nothing, varNames = nothing, title = nothing, plotKwargs = Dict(), showEvents = nothing)
        projFunc = projFunc === nothing ? MU["projNone"] : MU[projFunc]
        if varPairs === nothing @warn "You should define varPairs of Poincare projections" end
        new(projFunc, varNames, varPairs, title, plotKwargs, showEvents)
    end

    function (p::Poincare)(sol::py_sol; fig = nothing, savePath = nothing)
        if p.varNames !== nothing && length(p.varNames) < length(sol.y) throw(DimensionMismatch("varNames ($(length(p.varNames)) mismatch sol.y ($(length(sol.y))))")) end
        MU.plotPoincare(sol, p.varPairs, fig = fig, savePath = savePath, projFunc = p.projFunc, varNames = p.varNames, showEvents = p.showEvents, title = p.title, plotKwargs = p.plotKwargs)
    end
end




struct ReturnTime
    normFunc::Union{Any, Nothing}
    title::Union{String, Nothing}
    plotKwargs::Union{Dict, Nothing}

    function ReturnTime(properties::Dict = Dict())
        normFunc = get(properties, "normFunc", "normNone")
        title = get(properties, "title", nothing)
        plotKwargs = get(properties, "plotKwargs", Dict())

        ReturnTime(; normFunc, title, plotKwargs)
    end
    function ReturnTime(; normFunc = nothing, title = nothing, plotKwargs = Dict())
        normFunc = normFunc === nothing ? MU["normNone"] : MU[normFunc]
        new(normFunc, title, plotKwargs)
    end

    function (p::ReturnTime)(sol::py_sol; fig = nothing, savePath = nothing)
        if length(sol.t_events[1]) < 2 @warn "Cannot plot return time (perhaps, not enough event points)" end
        MU.plotReturnTime(sol, fig = fig, savePath = savePath, normFunc = p.normFunc, title = p.title, plotKwargs = p.plotKwargs)
    end
end




struct ActivationDiagram
    varNames::Union{Vector, Nothing}
    title::Union{String, Nothing}
    continuous::Union{Bool, Nothing}

    function ActivationDiagram(properties::Dict = Dict())
        varNames = get(properties, "varNames", nothing)
        title = get(properties, "title", nothing)
        continuous = get(properties, "continuous", false)

        ActivationDiagram(; varNames, title, continuous)
    end
    function ActivationDiagram(; varNames = nothing, title = nothing, continuous = false)
        new(varNames, title, continuous)
    end

    function (p::ActivationDiagram)(sol::py_sol; fig = nothing, savePath = nothing)
        if p.continuous
            MU.plotActivationDiagram_continuos(sol, fig = fig, savePath = savePath, title = p.title, varNames = p.varNames)
        else
            MU.plotActivationDiagram(sol, fig = fig, savePath = savePath, title = p.title, varNames = p.varNames)
        end
        
    end
end




struct Eigens
    projFunc::Union{Function, Nothing}
    title::Union{String, Nothing}
    plotKwargs::Union{Dict, Nothing}

    function Eigens(properties::Dict = Dict())
        projFunc = get(properties, "projFunc", nothing)
        title = get(properties, "title", nothing)
        plotKwargs = get(properties, "plotKwargs", Dict())

        Eigens(; projFunc, title, plotKwargs)
    end
    function Eigens(; projFunc = identity, title = nothing, plotKwargs = Dict())
        new(projFunc, title, plotKwargs)
    end

    function (p::Eigens)(sol::py_sol; fig = nothing, savePath = nothing)
        plotEigens(sol; fig, savePath, p.projFunc, p.title, p.plotKwargs)
    end
end

function plotEigens(sol; fig = nothing,
                        savePath = nothing,
                        projFunc = identity,
                        title = L"$\bf{Собственные~числа~матрицы~W}$",
                        plotKwargs = Dict())
    if (fig === nothing) fig = figure(figsize = (8, 4)); ax = fig.gca() else ax = fig.subplots(1, 1) end

    eigs = getTanEigens(sol)

    for Λ in eigs
        Λ = [projFunc(λ) for λ in Λ]
        ax.scatter(sol.t, Λ, s = 1)
    end
    ax.plot([sol.t[begin], sol.t[end]], [0.0, 0.0], color = "black", ls = "--", alpha = 0.5)

    ax.set_xlabel(L"t", fontsize=16)
    ax.set_ylabel(L"λ_i", fontsize=16)

    ax.set_title(title)
    if savePath !== nothing  fig.savefig(savePath) end

    return fig
end

function plotText(str; fig = nothing,
    savePath = nothing,
    title = nothing,
    plotKwargs = Dict())
    if (fig === nothing) fig = figure(); ax = fig.gca() else ax = fig.subplots(1, 1) end

    ax.set_axis_off()
    ax.text(0, 0, str)
end