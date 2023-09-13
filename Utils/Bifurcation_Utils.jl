using MAT
using Plots
using YAML: load_file

struct MatcontCurve
    x::Matrix{T} where T <: Number
    s::Dict{String, Any}

    function MatcontCurve(x::Matrix, s::Dict)
        new(x, s)
    end

    function MatcontCurve(matcont_data::M) where M <: MAT.MAT_v5.Matlabv5File
        x = read(matcont_data, "x")
        s = read(matcont_data, "s")
        MatcontCurve(x, s)
    end

    function MatcontCurve(matcont_filename::String)
        file = matopen(matcont_filename)
        MatcontCurve(file)
    end
end

function (mc::MatcontCurve)(idx::Int)
    if idx == 0
        p_1 = mc.s["data"][begin]["T"]
        p_2 = mc.s["data"][end]["T"]
    else
        p_1 = mc.s["data"][begin]["parametervalues"][idx]
        p_2 = mc.s["data"][end]["parametervalues"][idx]
    end

    if p_1 == p_2
        return fill(p_1, size(mc.x)[2])
    else
        for i in size(mc.x)[1] : -1 : 1
            if(mc.x[i, begin] == p_1 && mc.x[i, end] == p_2)
                return mc.x[i, :]
            end
        end
    end
end

@userplot BifCurve

@recipe function f(bc::BifCurve)
    if length(bc.args) != 2 || !(typeof(bc.args[1]) <: MatcontCurve) || !(typeof(bc.args[2]) <: Pair)
        error("Bifurcation curve should be given matcont_curve object and parameter indexes pair.  Got: $(typeof(bc.args))")
    else
        mc = bc.args[1]
        p = bc.args[2]
    end

    # set up properties for all subplots
    linewidth := get(plotattributes, :linewidth, 3)
    size := get(plotattributes, :size, (800, 800))
    legend := get(plotattributes, :legend, :none)
    framestyle := get(plotattributes, :framestyle, :box)
    
    xguide := "x"
    yguide := "y"

    @series begin
        seriestype := :path
        
        # Data for visualization
        mc(p[1]), mc(p[2])
    end
end








struct MatcontData
    curves::Vector{MatcontCurve}
    parameters::Dict{Any, Any}

    function MatcontData(dir_name::String = pwd())
        cd(dir_name)

        curves = Vector{MatcontCurve}()
        parameters = Dict{Any, Any}()

        for file_name in readdir()
            if contains(file_name, ".mat")
                push!(curves, MatcontCurve(file_name))
            elseif contains(file_name, ".yaml")
                parameters = load_file(file_name)
            end
        end
    
        cd("../")
        cd("../")
        new(curves, parameters)
    end
end


@userplot BifCurves

@recipe function f(bc::BifCurves)
    if length(bc.args) != 2 || !(typeof(bc.args[1]) <: MatcontData) || !(typeof(bc.args[2]) <: Pair)
        error("Bifurcation curve should be given matcont_curve object and parameter indexes pair.  Got: $(typeof(bc.args))")
    else
        md = bc.args[1]
        p_names = bc.args[2]
        p_idx = Pair(md.parameters[p_names[1]], md.parameters[p_names[2]])
    end

    # set up properties for all subplots
    linewidth := get(plotattributes, :linewidth, 3)
    size := get(plotattributes, :size, (800, 600))
    legend := get(plotattributes, :legend, :none)
    framestyle := get(plotattributes, :framestyle, :box)
    
    xguide := p_names[1]
    yguide := p_names[2]

    for mc in md.curves
        @series begin
            seriestype := :path
            
            # Data for visualization
            mc(p_idx[1]), mc(p_idx[2])
        end
    end
end