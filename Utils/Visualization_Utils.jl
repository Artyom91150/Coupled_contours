using PyCall

pushfirst!(PyVector(pyimport("sys")."path"), "")
MU = pyimport("Utils.Graphics_Utils")

struct TimeSeries
    projFunc::Union{Any, Nothing}
    varNames::Union{Vector, Nothing}
    title::Union{String, Nothing}
    plotKwargs::Union{Dict, Nothing}

    function TimeSeries(properties::Dict = Dict())
        projFunc = MU[get(properties, "projFunc", "projNone")]
        varNames = get(properties, "varNames", nothing)
        title = get(properties, "title", nothing)
        plotKwargs = get(properties, "plotKwargs", Dict())

        TimeSeries(; projFunc, varNames, title, plotKwargs)
    end
    function TimeSeries(; projFunc = nothing, varNames = nothing, title = nothing, plotKwargs = Dict())
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
        projFunc = MU[get(properties, "projFunc", "projNone")]
        varNames = get(properties, "varNames", nothing)
        varPairs = get(properties, "varPairs", nothing)
        title = get(properties, "title", nothing)
        plotKwargs = get(properties, "plotKwargs", Dict())
        showEvents = get(properties, "showEvents", nothing)

        Poincare(; varPairs, projFunc, varNames, title, plotKwargs, showEvents)
    end
    function Poincare(; varPairs = nothing, projFunc = nothing, varNames = nothing, title = nothing, plotKwargs = Dict(), showEvents = nothing)
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
        normFunc = MU[get(properties, "normFunc", "normNone")]
        title = get(properties, "title", nothing)
        plotKwargs = get(properties, "plotKwargs", Dict())

        ReturnTime(; normFunc, title, plotKwargs)
    end
    function ReturnTime(; normFunc = nothing, title = nothing, plotKwargs = Dict())
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






struct BS_FullPlot
    TS::TimeSeries
    RT::ReturnTime
    AD::ActivationDiagram
    Pcr::Poincare

    function BS_FullPlot(TS::Union{TimeSeries, Nothing} = nothing,
                         RT::Union{ReturnTime, Nothing} = nothing,
                         Pcr::Union{Poincare, Nothing} = nothing,
                         AD::Union{ActivationDiagram, Nothing} = nothing)
        # Time series
        if isnothing(TS)
            projFunc = "projCos"
            varNames = ["\\phi_1", "\\phi_2", "\\phi_3", "\\psi_1", "\\psi_2", "\\psi_3"]
            title = "Time Series"
            plotKwargs = Dict(:linewidth => 3)
    
            TS = TimeSeries(Dict("projFunc" => projFunc, "varNames" => varNames, "title" => title, "kwargs" => plotKwargs));
        end

        # Poincare & projections
        if isnothing(Pcr)
            title = "Poincare"
            plotKwargs = Dict(:linewidth => 1)
            showEvents = true
            varPairs = [[[0, 1], [0, 2], [1, 2]],
                        [[3, 4], [3, 5], [4, 5]],
                        [[0, 3], [0, 4], [0, 5]],
                        [[1, 3], [1, 4], [1, 5]],
                        [[2, 3], [2, 4], [2, 5]]]

            Pcr = Poincare(Dict("projFunc" => projFunc,
                                    "varNames" => varNames,
                                    "title" => title,
                                    "showEvents" => showEvents,
                                    "varPairs" => varPairs,
                                    "kwargs" => plotKwargs))
        end

        # Return time
        if isnothing(RT)
            normFunc = "normNone"
            title = "Return Time"
            plotKwargs = Dict(:linewidth => 3)

            RT = ReturnTime(Dict("normFunc" => normFunc, "title" => title, "plotKwargs" => plotKwargs))
        end

        # Activation diagram
        if isnothing(AD)
            title = "Activation Diagram"
            continuous = true

            AD = ActivationDiagram(Dict("varNames" => varNames, "title" => title, "continuous" => continuous))
        end

        new(TS, RT, AD, Pcr)
    end

    function (p::BS_FullPlot)(trj::Trajectory; fig = nothing, savePath = nothing)
        if (fig === nothing) fig = plt.figure(figsize=(22, 16)); ax = fig.gca() else ax = fig.subplots(1, 1) end
        subfigs1 = fig.subfigures(1, 2, width_ratios=[2, 1], wspace=0)

        left_subfigs = subfigs1[1].subfigures(4, 1, wspace=0.07, height_ratios = [1, 1, 1, 1.5])
        right_subfigs = subfigs1[2].subfigures(2, 1, height_ratios = [5, 1])

        sol = trj.solution

        p.TS(sol; fig = left_subfigs[1])
        p.AD(sol; fig = left_subfigs[2])
        p.RT(sol; fig = left_subfigs[3])
        plotText(repr(trj); fig = left_subfigs[4])

        p.Pcr(sol; fig = right_subfigs[1])
        plotNorms(sol; fig = right_subfigs[2], skip_points = 10, norm_arg = 2)
    end
end







struct BS_FullPlot_Tan
    TS::TimeSeries
    EG::Eigens

    function BS_FullPlot_Tan(;TS::Union{TimeSeries, Nothing} = nothing,
                         EG::Union{Eigens, Nothing} = nothing)
        # Time series
        if isnothing(TS)
            projFunc = "projCos"
            varNames = ["\\phi_1", "\\phi_2", "\\phi_3", "\\psi_1", "\\psi_2", "\\psi_3"]
            title = "Time Series"
            plotKwargs = Dict(:linewidth => 3)
    
            TS = TimeSeries(Dict("projFunc" => projFunc, "varNames" => varNames, "title" => title, "kwargs" => plotKwargs));
        end

        # Eigen values
        if isnothing(EG)
            EG = Eigens()
        end

        new(TS, EG)
    end

    function (p::BS_FullPlot_Tan)(trj::Trajectory; fig = nothing, savePath = nothing)
        if (fig === nothing) fig = plt.figure(figsize=(16, 18)); ax = fig.gca() else ax = fig.subplots(1, 1) end
        subfigs1 = fig.subfigures(3, 1)

        ode_sol, tan_sol = splitSol(tan_trj.solution)

        p.TS(ode_sol; fig = subfigs1[1])
        p.EG(tan_sol; fig = subfigs1[2])
        subfigs1[2].gca().axvline(ode_sol.t[end], color = "red", ls = "--", alpha = 0.5)
        #subfigs1[2].gca().set_ylabel(L"log_{10}(|\lambda_i|)")

        subfigs2 = subfigs1[3].subfigures(1, 2, width_ratios = [2, 1])
        plotText(repr(trj); fig = subfigs2[1])
        plotText(repr("text/plain", [y[end] for y in getTanEigens(tan_sol)]); fig = subfigs2[2])
    end
end