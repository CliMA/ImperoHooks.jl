export findelement

# 3D version
function findelement(xC, yC, zC, location, p, lin)
    ex, ey, ez = size(lin)
    # i 
    currentmin = ones(1)
    minind = ones(Int64, 1)
    currentmin[1] = abs.(xC[p[lin[1,1,1]]] .- location[1])
    for i in 2:ex
        current = abs.(xC[p[lin[i,1,1]]] .- location[1])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    i = minind[1]
    # j 
    currentmin[1] = abs.(yC[p[lin[1,1,1]]] .- location[2])
    minind[1] = 1
    for i in 2:ey
        current = abs.(yC[p[lin[1,i,1]]] .- location[2])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    j = minind[1]
    # k 
    currentmin[1] = abs.(zC[p[lin[1,1,1]]] .- location[3])
    minind[1] = 1
    for i in 2:ez
        current = abs.(zC[p[lin[1,1, i]]] .- location[3])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    k = minind[1]
    return p[lin[i,j,k]]
end

# 2D version
function findelement(xC, yC, location, p, lin)
    ex, ey = size(lin)
    # i 
    currentmin = ones(1)
    minind = ones(Int64, 1)
    currentmin[1] = abs.(xC[p[lin[1,1]]] .- location[1])
    for i in 2:ex
        current = abs.(xC[p[lin[i,1]]] .- location[1])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    i = minind[1]
    # j 
    currentmin[1] = abs.(yC[p[lin[1,1]]] .- location[2])
    minind[1] = 1
    for i in 2:ey
        current = abs.(yC[p[lin[1,i]]] .- location[2])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    j = minind[1]
    return p[lin[i,j]]
end