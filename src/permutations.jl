export getperm

# computes permutation to cartesian index for elements
# needs lin = reshape(collect(1:length(xC)), (ex, ey, ez))
# but then xC[p[lin[i,j,k]]] works where i,j,k are cartesian
function getperm(xC, yC, zC, ex, ey, ez)
    pz = sortperm(zC)
    tmpY = reshape(yC[pz], (ex*ey, ez))
    tmp_py = [sortperm(tmpY[:,i]) for i in 1:ez]
    py = zeros(Int64,length(pz))
    for i in eachindex(tmp_py)
        n = length(tmp_py[i]) 
        ii = (i-1) * n + 1
        py[ii : ii + n-1] .= tmp_py[i] .+ ii .- 1
    end
    tmpX = reshape(xC[pz][py], (ex, ey*ez))
    tmp_px = [sortperm(tmpX[:,i]) for i in 1:ey*ez]
    px = zeros(Int64,length(pz))
    for i in eachindex(tmp_px)
        n = length(tmp_px[i]) 
        ii = (i-1) * n + 1
        px[ii : ii + n-1] .= tmp_px[i] .+ ii .- 1
    end
    p = [pz[py[px[i]]] for i in eachindex(px)]
    return p
end
