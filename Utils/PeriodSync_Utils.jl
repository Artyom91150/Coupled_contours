function normArctan(x)
    return atan(sin(x), cos(x))
end

function circleDist(x, y)
    x = normArctan(x)
    y = normArctan(y)
    return min(abs(x - y), 2*pi - abs(x - y))
end

function getNorms(sol)
    norms = similar(sol.t)
    for i = 2 : length(sol.t)
        norms[i] = log10(norm([circleDist(y[1], y[i]) for y in sol.y], Inf)) 
    end
    return norms
end

function getPeriodTime(sol)
    norms = getNorms(sol)
    sums = similar(norms)
    min_filter = [5, -10, 5]
    for i = 1 : length(norms) - 2
        sums[i] = sum(norms[i : i + 2] .* min_filter)
    end
    return mean([p for p in diff([sol.t[i] for (i, s) in enumerate(sums) if s >= 1]) if p >= 1])
end

function plotNorms(sol; savePath = nothing)
    fig = plt.figure(figsize=(6, 4))
    norms = getNorms(sol)
    plot(sol.t, norms)
    gca().set_xlabel("t")
    gca().set_ylabel("log10( || circleDist(y0, yt) || )")
    plt.title(string(getPeriodTime(sol)))

    savePath !== nothing ? fig.savefig(savePath) : nothing
end