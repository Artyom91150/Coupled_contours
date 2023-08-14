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