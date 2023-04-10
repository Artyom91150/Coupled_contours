function normArctan(x::Vector{T}) where T <: Real
    return [normArctan(y) for y in x]
end

function normArctan(x::T) where T <: Real
    return atan(sin(x), cos(x))
end

function circleDist(x::Vector{T}, y::Vector{T}) where T <: Real
    if length(x) != length(y) throw(DimensionMismatch) end
    return [circleDist(x[i], y[i]) for i = 1:length(x)]
end

function circleDist(x::T, y::T) where T <: Real
    x = normArctan(x)
    y = normArctan(y)
    return min(abs(x - y), 2*pi - abs(x - y))
end



function getNorms(sol; skip_points = 0, norm_arg = 2)
    norms = fill(NaN, length(sol.t))
    for i = 1 + skip_points : length(norms)
        norms[i] = log10(norm([circleDist(y[1], y[i]) for y in sol.y], norm_arg)) 
    end
    return norms
end

function getPeriodTime(sol; skip_points = 0, norm_arg = 2)
    norms = getNorms(sol;  skip_points = skip_points, norm_arg = norm_arg)
    sums = similar(norms)
    min_filter = [20.0, 10.0, -60.0, 10.0, 20.0]
    for i = 1 : length(norms) - length(min_filter) + 1
        sums[i] = sum(norms[i : i + length(min_filter) - 1] .* min_filter)
    end
    return mean([p for p in diff([sol.t[i] for (i, s) in enumerate(sums) if s >= 3]) if p >= 1])
end

function RefinePeriod(trj::Trajectory, t)
    func = x -> norm(circleDist(trj(trj.solution.t[begin]), trj(trj.solution.t[begin] + x[1])))
    s = optimize(func, [t])
    return s#.minimizer[1]
end

function plotNorms(sol; fig = nothing, skip_points = 0, norm_arg = 2, savePath = nothing)
    if (fig === nothing) fig = figure(figsize = (8, 4)); ax = fig.gca() else ax = fig.subplots(1, 1) end
    norms = getNorms(sol; skip_points = skip_points, norm_arg = norm_arg)
    plot(sol.t, norms)
    ax.set_xlabel("t")
    ax.set_ylabel("log10( || circleDist(y0, yt) || )")
    ax.set_title("~" * string(getPeriodTime(sol; skip_points = skip_points, norm_arg = norm_arg)))

    savePath !== nothing ? fig.savefig(savePath) : nothing
end