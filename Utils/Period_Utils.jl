function normArctan(x)
    return atan(sin(x), cos(x))
end

function circleDist(x, y)
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

function plotNorms(sol; skip_points = 0, norm_arg = 2, savePath = nothing)
    fig = plt.figure(figsize=(6, 4))
    norms = getNorms(sol; skip_points = skip_points, norm_arg = norm_arg)
    plot(sol.t, norms)
    gca().set_xlabel("t")
    gca().set_ylabel("log10( || circleDist(y0, yt) || )")
    plt.title(string(getPeriodTime(sol; skip_points = skip_points, norm_arg = norm_arg)))

    savePath !== nothing ? fig.savefig(savePath) : nothing
end