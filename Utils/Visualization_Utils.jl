using PyCall

pushfirst!(PyVector(pyimport("sys")."path"), "")
MU = pyimport("Utils.Graphics_Utils")

function plotTimeSeries(sol; plot_properties = Dict(), Folder_path::String = "", plot_Name::String = "TimeSeries")
    projFunc = haskey(plot_properties, "projFunc") ? MU[plot_properties["projFunc"]] : nothing
    varNames = haskey(plot_properties, "varNames") ? plot_properties["varNames"] : nothing
    title = haskey(plot_properties, "title") ? plot_properties["title"] : nothing
    plotKwargs = haskey(plot_properties, "kwargs") ? plot_properties["kwargs"] : nothing
    savePath = Folder_path === "" ? nothing : Folder_path * "\\" * plot_Name

    MU.plotTimeSeries(sol, savePath = savePath, projFunc = projFunc, varNames = varNames, title = title, plotKwargs = plotKwargs)
end

function plotPoincare(sol; plot_properties = Dict(), Folder_path::String = "", plot_Name::String = "Poincare")
    projFunc = haskey(plot_properties, "projFunc") ? MU[plot_properties["projFunc"]] : nothing
    varNames = haskey(plot_properties, "varNames") ? plot_properties["varNames"] : nothing
    title = haskey(plot_properties, "title") ? plot_properties["title"] : nothing
    plotKwargs = haskey(plot_properties, "kwargs") ? plot_properties["kwargs"] : nothing
    savePath = Folder_path === "" ? nothing : Folder_path * "\\" * plot_Name

    showEvents = haskey(plot_properties, "showEvents") ? plot_properties["showEvents"] : nothing
    varPairs = plot_properties["varPairs"]

    MU.plotPoincare(sol, varPairs, savePath = savePath, projFunc = projFunc, varNames = varNames, showEvents = showEvents, title = title, plotKwargs = plotKwargs)
end

function plotReturnTime(sol; plot_properties = Dict(), Folder_path::String = "", plot_Name::String = "Return Time")
    normFunc = haskey(plot_properties, "normFunc") ? MU[plot_properties["normFunc"]] : nothing
    title = haskey(plot_properties, "title") ? plot_properties["title"] : nothing
    plotKwargs = haskey(plot_properties, "kwargs") ? plot_properties["kwargs"] : nothing
    savePath = Folder_path === "" ? nothing : Folder_path * "\\" * plot_Name

    try
        MU.plotReturnTime(sol, savePath = savePath, normFunc = normFunc, title = title, plotKwargs = plotKwargs)
    catch
        @warn "Cannot plot return time (perhaps, not enough event points)"
    end
end

function plotActivationDiagram(sol; plot_properties = Dict(), Folder_path::String = "", plot_Name::String = "Activation Diagram")
    varNames = haskey(plot_properties, "varNames") ? plot_properties["varNames"] : nothing
    title = haskey(plot_properties, "title") ? plot_properties["title"] : nothing
    savePath = Folder_path === "" ? nothing : Folder_path * "\\" * plot_Name

    MU.plotActivationDiagram(sol, varNames = varNames, savePath = savePath, title = title)
end