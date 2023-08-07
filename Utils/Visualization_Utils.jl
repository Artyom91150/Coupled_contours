#using PyCall

#pushfirst!(PyVector(pyimport("sys")."path"), "")
#MU = pyimport("Utils.Graphics_Utils")
#mpc = pyimport("matplotlib.colors")
#gridspec = pyimport("matplotlib.gridspec")

#font_family = "times"
#font_family = "georgia"
font_family = "arial"
#font_family = "tahoma"
#font_family = "verdana"
#font_family = "sans"
#font_family = "impact"
#font_family = "courier"
#font_family = "lucida"

default(titlefont = (12, font_family), guidefont = (16, font_family), tickfont = (10, font_family), color_palette = :tab10)

function MakeLabels(::Nothing; proj_func = identity, dims = 1)
    return MakeLabels("y"; proj_func, dims)
end

function MakeLabels(label::String; proj_func = identity, dims = 1)
    return [latexstring(proj_func(label * "_$i")) for i in 1 : dims]
end

function MakeLabels(labels::Vector{String}; proj_func = identity, dims = 1)
    if length(labels) != dims @warn "labels dimension mismatch." end
    return [latexstring(proj_func(l)) for l in labels]
end

function set_margin!(plotattributes; scale = 1.0)
    plotattributes[:left_margin] = get(plotattributes, :left_margin, 10Plots.Measures.mm * scale)
    plotattributes[:right_margin] = get(plotattributes, :left_margin, 5Plots.Measures.mm * scale)
    plotattributes[:top_margin] = get(plotattributes, :left_margin, 2Plots.Measures.mm * scale)
    plotattributes[:bottom_margin] = get(plotattributes, :left_margin, 2Plots.Measures.mm * scale)

    
end

function getTrigTiks()
    return ([0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 2*pi], [L"0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi", L"\frac{5\pi}{4}", L"\frac{3\pi}{2}", L"\frac{7\pi}{4}", L"2\pi"])
end

makeLabel(x) = permutedims(hcat(x))

projCos(x::Float64) = cos(x)
projCos(x::String) = "\\cos{$x}"

projCos2(x::Float64) = cos(x/2.0)
projCos2(x::String) = "\\cos{$x/2}"

projLogCos(x::Float64) = x > 0 ? -log10(1e-15 + 1 - cos(x)) : log10(1e-15 + 1 + cos(x))
projLogCos(x::String) = L"\pm log_{10} (1 \pm cos(x))"

projLogSin2(x::Float64) = x > 0 ? -log10(1e-15 + 1 - sin(x/2)) : log10(1e-15 + 1 + sin(x/2))
projLogSin2(x::String) = L"\pm log_{10} (1 \pm sin(x/2))"

#=
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
=#

@userplot TimeSeries

@recipe function f(ts::TimeSeries; var_labels = nothing, proj_func = identity, enable_margin = true)
    if length(ts.args) != 1 || !(typeof(ts.args[1]) <: py_sol) 
        error("Time Series should be given py_sol object.  Got: $(typeof(ts.args))")
    else
        sol = ts.args[1]
    end

    var_labels = MakeLabels(var_labels; proj_func, dims = length(sol.y))

    # set up properties for all subplots
    linewidth := get(plotattributes, :linewidth, 3)
    size := get(plotattributes, :size, (800, 800))
    legend := get(plotattributes, :legend, :none)
    framestyle := get(plotattributes, :framestyle, :box)
    layout := (length(sol.y), 1)

    #if enable_margin set_margin!(plotattributes) end
    
    left_margin := get(plotattributes, :left_margin, 10Plots.Measures.mm)
    right_margin := get(plotattributes, :left_margin, 5Plots.Measures.mm)

    for (i, y) in enumerate(sol.y)
        # X axis properties
        xguide := ifelse(i == length(sol.y), "\$t\$", "")
        xticks := ifelse(i == length(sol.y), :auto, :none)
        xgrid := :none

        # Y axis properties
        yguide := var_labels[i]

        # Other axes properties
        left_margin := (50, :px)
        right_margin := (50, :px)

        # set up properties for current subplot
        @series begin
            seriestype := :path
            subplot := i
            
            # Data for visualization
            sol.t, proj_func.(sol.y[i])
        end
    end
end



#=
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
=#

@userplot Poincare

@recipe function f(p::Poincare; var_labels = nothing, proj_func = identity, enable_events = false, enable_margin = true)
    if length(p.args) != 2 || !(typeof(p.args[1]) <: py_sol)  || !(typeof(p.args[2]) <: Matrix{Tuple{Int64, Int64}})
        error("Poincare should be given py_sol object and matrix with variables pairs (Matrix{Tuple{Int64, Int64}}).  Got: $(typeof(p.args))")
    else
        sol = p.args[1]
        var_pairs = p.args[2]
    end

    # Data preprocessing
    var_labels = MakeLabels(var_labels; proj_func, dims = length(sol.y))
    y = [proj_func.(y) for y in sol.y]
    ax_lims = (1.2 * minimum(minimum.(y)), 1.2* maximum(maximum.(y)))

    # set up properties for all subplots
    if enable_margin set_margin!(plotattributes; scale = 0.25) end
    linewidth := get(plotattributes, :linewidth, 3)
    size := get(plotattributes, :size, (800, 800))
    legend := get(plotattributes, :legend, :none)
    framestyle := get(plotattributes, :framestyle, :box)
    layout := size(var_pairs)


    #link := :none
    aspect_ratio := :equal

    for (i, pair) in enumerate(var_pairs)

        # X axis properties
        xguide := var_labels[pair[1]]

        # Y axis properties
        yguide := var_labels[pair[2]]

        # Other axes properties
        xlims := get(plotattributes, :xlims, ax_lims)
        ylims := get(plotattributes, :ylims, ax_lims)

        # set up properties for current subplot
        @series begin
            seriestype := :path
            subplot := i

            #alpha := ifelse(enable_events, 0.5, get(plotattributes, :alpha, 1.0))
            
            # Data for visualization
            y[pair[1]], y[pair[2]]
        end

        if enable_events && length(sol.t_events[1]) > 0
            @series begin
                seriestype := :scatter
                subplot := i
                
                markeralpha := get(plotattributes, :markeralpha, :true)
                markercolor := palette(get(plotattributes, :markercolor, :auto))[2]
                markersize := get(plotattributes, :markersize, 6.0)
                markerstrokewidth := get(plotattributes, :markerstrokewidth, 0.0)

                proj_func.(sol.y_events[1][:, pair[1]]), proj_func.(sol.y_events[1][:, pair[2]])
            end
        end
    end
end



#=
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
=#



@userplot ReturnTime

@recipe function f(rt::ReturnTime; enable_margin = true)
    if length(rt.args) != 1 || !(typeof(rt.args[1]) <: py_sol) 
        error("ReturnTime should be given py_sol object.  Got: $(typeof(rt.args))")
    else
        sol = rt.args[1]
    end

    # Data preprocessing
    diffs = ifelse(length(sol.t_events[1]) > 2, sol.t_events[1][2 : end] - sol.t_events[1][1 : end - 1], [])

    # set up properties for plot
    size := get(plotattributes, :size, (800, 400))
    legend := get(plotattributes, :legend, :none)
    framestyle := get(plotattributes, :framestyle, :box)

    if enable_margin set_margin!(plotattributes) end

    @series begin
        seriestype := :scatter

        # X axis properties
        xguide := "\$t\$"

        # Y axis properties
        yguide := "\$R\$"

        markeralpha := get(plotattributes, :markeralpha, :true)
        markercolor := get(plotattributes, :markercolor, :auto)
        markersize := get(plotattributes, :markersize, 6.0)
        markerstrokewidth := get(plotattributes, :markerstrokewidth, 0.0)

        # Data for visualization
        sol.t_events[1][2 : end], diffs
    end
end



#=
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
=#

@userplot ActivationDiagram

@recipe function f(ad::ActivationDiagram; var_labels = nothing, enable_margin = true, color_eps = nothing)
    if length(ad.args) != 1 || !(typeof(ad.args[1]) <: py_sol) 
        error("ActivationDiagram should be given py_sol object.  Got: $(typeof(ad.args))")
    else
        sol = ad.args[1]
    end

    # Data preprocessing
    proj_func = projCos
    var_labels = reverse(MakeLabels(var_labels; proj_func, dims = length(sol.y)))
    z = (1.0 .+ cos.(reduce(hcat, reverse(sol.y)))') ./ 2.0

    # set up properties for plot
    size := get(plotattributes, :size, (1200, 600))
    legend := get(plotattributes, :legend, :none)
    framestyle := get(plotattributes, :framestyle, :box)
    ytickfontsize := get(plotattributes, :guidefontsize, default(:yguidefontsize))


    if enable_margin set_margin!(plotattributes) end

    @series begin
        seriestype := :heatmap

        # X axis properties
        xguide := L"t"

        # Y axis properties
        
        # Colorbar properties
        colorbar_tickfontsize := 1
        if !isnothing(color_eps)
            color := cgrad(get(plotattributes, :color, :greys), [color_eps, 1.0 - color_eps], rev = true, categorical = true, scale = :log)
            #colorbar_ticks := [0.0, color_eps, 1.0 - color_eps, 1.0]
        end
        

        # Data for visualization
        sol.t, var_labels, z
    end
end


#=
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
                        title = latexify("\bf{Собственные~числа~матрицы~W}"),
                        plotKwargs = Dict())
    if (fig === nothing) fig = figure(figsize = (8, 4)); ax = fig.gca() else ax = fig.subplots(1, 1) end

    eigs = getTanEigens(sol)

    for Λ in eigs
        Λ = [projFunc(λ) for λ in Λ]
        ax.scatter(sol.t, Λ, s = 1)
    end
    ax.plot([sol.t[begin], sol.t[end]], [0.0, 0.0], color = "black", ls = "--", alpha = 0.5)

    ax.set_xlabel(latexify("t"), fontsize=16)
    ax.set_ylabel(latexify("λ_i"), fontsize=16)

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
=#

@userplot EigensPlot

@recipe function f(ep::EigensPlot; enable_margin = true)
    if length(ep.args) != 1 || !(typeof(ep.args[1]) <: py_sol) 
        error("EigensPlot should be given py_sol object.  Got: $(typeof(ep.args))")
    else
        sol = ep.args[1]
    end

    # Data preprocessing
    eigs = [abs.(eig) .+ 1e-15 for eig in getTanEigens(sol)]

    # set up properties for plot
    size := get(plotattributes, :size, (800, 400))
    legend := get(plotattributes, :legend, :none)
    framestyle := get(plotattributes, :framestyle, :box)
    yscale := :log10

    if enable_margin set_margin!(plotattributes) end

    # X axis properties
    xguide := "\$t\$"

    # Y axis properties
    yguide := "\$λ\$"

    @series begin
        seriestype := :scatter

        markeralpha := get(plotattributes, :markeralpha, :true)
        markercolor := get(plotattributes, :markercolor, :auto)
        markersize := get(plotattributes, :markersize, 1.0)
        markerstrokewidth := get(plotattributes, :markerstrokewidth, 0.0)

        # Data for visualization
        sol.t, eigs
    end

    @series begin
        seriestype := :hline
        color := :black
        linestyle := :dash
        [1.0]
    end
end




#=
function LLE_plot(λs, Eps_mesh)
    fig = figure(figsize=(8, 4))
        
    for i in 1:size(λs)[1]
        scatter(Eps_mesh, λs[i, :], zorder = 101, s = 5)
        plot(Eps_mesh, λs[i, :], linewidth = 2, alpha = 0.5)
    end
    plot([Eps_mesh[1], Eps_mesh[end]], [0, 0], "k--", alpha = 0.2, zorder = 99)
    
    myCmap = mpc.ListedColormap([[1, 0.7, 0.7], [1, 1, 1], [0.7, 1, 0.7]]);
    ColorVal = zeros(length(Eps_mesh), 2)

    for i = 1 : length(Eps_mesh)
        if any(i -> abs(i) < 1e-3, λs[:, i])
            ColorVal[i, :] = ColorVal[i, :] + [1, 1]
            if any(i -> i > 1e-3, λs[:, i])
                ColorVal[i, :] = ColorVal[i, :] + [1, 1]
            end
        end
    end

    min_λ = minimum(x -> isnan(x) ? Inf : x, λs)
    max_λ = maximum(x -> isnan(x) ? -Inf : x, λs)
    
    plt.pcolor(Eps_mesh, [1.1 * min_λ, 1.1 * max_λ], ColorVal', cmap = myCmap, vmin=minimum(0), vmax=maximum(2))
    plot([Eps_mesh[1], Eps_mesh[end]], [max_λ, max_λ], "g--", alpha = 0.4, zorder = 99)
    
    xlabel("\$\\varepsilon\$", fontsize=20)
    ylabel("\$\\lambda\$", fontsize=20)

    fig.tight_layout(pad=0.3)
    fig.gca().set_xscale("log")

    tick_params(which="major", width=1.0, labelsize=14)
    tick_params(which="major", length=5, labelsize=14)

    plt.grid(true, alpha = 0.2)

    return fig
end
=#

@userplot LEPlot

@recipe function f(le::LEPlot; enable_margin = true)
    if length(le.args) != 2 || !(typeof(le.args[1]) <: AbstractArray) || !(typeof(le.args[2]) <: AbstractArray) 
        error("LEPlot should be given array of LEs and mesh vector.  Got: $(typeof(le.args))")
    else
        λs = le.args[1]
        mesh = le.args[2]
    end

    # Data preprocessing
    λs_vec = [λs[i, :] for i in 1 : size(λs)[1]]
    # set up properties for plot
    size := get(plotattributes, :size, (800, 400))
    legend := get(plotattributes, :legend, :none)
    framestyle := get(plotattributes, :framestyle, :box)
    xscale := get(plotattributes, :xsacle, :log10)

    if enable_margin set_margin!(plotattributes, scale = 0.5) end

    # X axis properties
    xguide := "\$\\varepsilon\$"

    # Y axis properties
    yguide := "\$λ\$"

    color_palette  := palette(get(plotattributes, :color_palette, default(:color_palette)))[1 : length(λs_vec)]


    @series begin
        seriestype := :line

        alpha := 0.8

        # Data for visualization
        mesh, λs_vec
    end

    @series begin
        seriestype := :scatter

        markeralpha := get(plotattributes, :markeralpha, :true)
        markercolor := get(plotattributes, :markercolor, :auto)
        markersize := get(plotattributes, :markersize, 3.0)
        markerstrokewidth := get(plotattributes, :markerstrokewidth, 0.0)

        # Data for visualization
        mesh, λs_vec
    end

    yscale = get(plotattributes, :yscale, :none)
    zero_line = [0.0]
    if yscale == :log10
        zero_line = [1.0]
    end

    @series begin
        seriestype := :hline
        color := :black
        linestyle := :dash
        zero_line
    end
end