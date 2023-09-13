function get_ls_ip(start, stop, center, size)
    left_scale = abs(center - start) / (stop - start)

    left_size = max(round(Int, size * left_scale), 2)
    right_size = size - left_size
    
    return ([range(start, center, length = left_size+1)[begin:end-1]; center; range(center, stop, length = right_size+1)[begin+1:end]], left_size + 1)
end

mutable struct ICPMEsh
    p_mesh::Matrix
    ics::Matrix
    idx0::Tuple{Int, Int}

    function ICPMEsh(
        left::Number,
        right::Number,
        bottom::Number,
        top::Number,
        center::Tuple{Number, Number, T},
        mesh_size::Int
        ) where T

        frst_ls, i = get_ls_ip(left, right, center[1], mesh_size)
        scnd_ls, j = get_ls_ip(bottom, top, center[2], mesh_size)
        mesh = [[x, y] for x in frst_ls, y in scnd_ls]
        ics = fill(T(zeros(length(center[3]))), size(mesh))
        ics[i, j] = center[3]
        new(mesh, ics, (i, j))
    end
end

function get_prev_idx(i, j)
    i_sign = sign(i); i = abs(i)
    j_sign = sign(j); j = abs(j)
    if i > j
        return [i_sign * (i - 1), j_sign * j]
    elseif i < j
        return [i_sign * i, j_sign * (j - 1)]
    else
        return [i_sign * (i - 1), j_sign * (j - 1)]
    end
end

function get_radius(icpm::ICPMEsh)
    return maximum([icpm.idx0[1] - 1, icpm.idx0[2] - 1, size(icpm.p_mesh)[1] - icpm.idx0[1], size(icpm.p_mesh)[1] - icpm.idx0[2]])
end

function get_rad_idx(r)
    pts = Vector{Vector{Int}}()

    for i in -r:1:r
        for j in -r:1:r
            if abs(i) == r || abs(j) == r
                push!(pts, [i, j])
            end
        end
    end
    return pts
end

import Base: getindex
function Base.getindex(cm::ICPMEsh, idx...)
    i = idx[1]; j = idx[2]
    shifted_i = i + cm.idx0[1]
    shifted_j = j + cm.idx0[2]
    if 0 < shifted_i <= size(cm.p_mesh)[1] && 0 < shifted_j <= size(cm.p_mesh)[2]
        prev_i, prev_j = get_prev_idx(i, j) .+ cm.idx0
        return (cm.p_mesh[shifted_i, shifted_j], cm.ics[prev_i, prev_j])
    else
        return nothing
    end
end

import Base: setindex!
function Base.setindex!(cm::ICPMEsh, val, idx...)
    shifted_idx = idx .+ cm.idx0
    if 0 < shifted_idx[1] <= size(cm.p_mesh)[1] && 0 < shifted_idx[2] <= size(cm.p_mesh)[2]
        return setindex!(cm.ics, val, shifted_idx...)
    else
        return nothing
    end
end