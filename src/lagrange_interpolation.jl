export lagrange_eval, getvalue

function checkgl(x, rx)
    for i in eachindex(rx)
        if abs(x - rx[i]) ≤ eps(rx[i])  return i end
    end
    return 0
end

function lagrange_eval(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
    icheck = checkgl(newx, rx)
    jcheck = checkgl(newy, ry)
    kcheck = checkgl(newz, rz)
    numerator = zeros(1)
    denominator = zeros(1)
    for k in eachindex(rz)
        if kcheck ==0
            Δz = (newz .- rz[k])
            polez = ωz[k] ./ Δz
            kk = k
        else
            polez = 1.0
            k = eachindex(rz)[end]
            kk = kcheck
        end
        for j in eachindex(ry)
            if jcheck ==0
                Δy = (newy .- ry[j])
                poley = ωy[j] ./ Δy
                jj = j
            else
                poley = 1.0
                j = eachindex(ry)[end]
                jj = jcheck
            end
            for i in eachindex(rx)
                if icheck == 0
                    Δx = (newx .- rx[i])
                    polex = ωx[i] ./ Δx
                    ii = i
                else
                    polex = 1.0
                    i = eachindex(rx)[end]
                    ii = icheck
                end
                numerator[1] += f[ii,jj,kk] * polex * poley * polez
                denominator[1] += polex * poley * polez
            end
        end
    end
    return numerator[1] / denominator[1]
end

function lagrange_eval(f, newx, newy, rx, ry, ωx, ωy)
    icheck = checkgl(newx, rx)
    jcheck = checkgl(newy, ry)
    numerator = zeros(1)
    denominator = zeros(1)

    for j in eachindex(ry)
        if jcheck ==0
            Δy = (newy .- ry[j])
            poley = ωy[j] ./ Δy
            jj = j
        else
            poley = 1.0
            j = eachindex(ry)[end]
            jj = jcheck
        end
        for i in eachindex(rx)
            if icheck == 0
                Δx = (newx .- rx[i])
                polex = ωx[i] ./ Δx
                ii = i
            else
                polex = 1.0
                i = eachindex(rx)[end]
                ii = icheck
            end
            numerator[1] += f[ii,jj] * polex * poley 
            denominator[1] += polex * poley 
        end
    end
    return numerator[1] / denominator[1]
end


function lagrange_eval_nocheck(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
    numerator = zeros(1)
    denominator = zeros(1)
    for k in eachindex(rz)
            Δz = (newz .- rz[k])
            polez = ωz[k] ./ Δz
        for j in eachindex(ry)
                Δy = (newy .- ry[j])
                poley = ωy[j] ./ Δy
            for i in eachindex(rx)
                    Δx = (newx .- rx[i])
                    polex = ωx[i] ./ Δx
                    numerator[1] += f[i,j,k] * polex * poley * polez
                    denominator[1] += polex * poley * polez
            end
        end
    end
    return numerator[1] / denominator[1]
end

function lagrange_eval_nocheck(f, newx, newy, rx, ry, ωx, ωy)
    numerator = zeros(1)
    denominator = zeros(1)
    for j in eachindex(ry)
            Δy = (newy .- ry[j])
            poley = ωy[j] ./ Δy
        for i in eachindex(rx)
                Δx = (newx .- rx[i])
                polex = ωx[i] ./ Δx
                numerator[1] += f[i,j] * polex * poley
                denominator[1] += polex * poley
        end
    end
    return numerator[1] / denominator[1]
end

function lagrange_eval_nocheck(f, newx, rx, ωx)
    numerator = zeros(1)
    denominator = zeros(1)
    for i in eachindex(rx)
            Δx = (newx .- rx[i])
            polex = ωx[i] ./ Δx
            numerator[1] += f[i] * polex * poley
            denominator[1] += polex * poley
    end
    return numerator[1] / denominator[1]
end

# 3D
function getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
    e = findelement(xC, yC, zC, location, p, lin)
    # need bounds to rescale
    
    xmax = x[linlocal[length(rx),1,1], e]
    xmin = x[linlocal[1,1,1], e]
    ymax = y[linlocal[1,length(ry),1], e]
    ymin = y[linlocal[1,1,1], e]
    zmax = z[linlocal[1,1,length(rz)], e]
    zmin = z[linlocal[1,1,1], e]

    # rescale new point to [-1,1]³
    newx = 2 * (location[1] - xmin) / (xmax - xmin) - 1
    newy = 2 * (location[2] - ymin) / (ymax - ymin) - 1
    newz = 2 * (location[3] - zmin) / (zmax - zmin) - 1

    return lagrange_eval(view(fl,:,:,:,e), newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
end

# 2D
function getvalue(fl, xC, yC, location, p, lin, linlocal, x, y, rx, ry,  ωx, ωy)
    e = findelement(xC, yC, location, p, lin)
    # need bounds to rescale
    xmax = x[linlocal[length(rx),1,1], e]
    xmin = x[linlocal[1,1,1], e]
    ymax = y[linlocal[1,length(ry),1], e]
    ymin = y[linlocal[1,1,1], e]

    # rescale new point to [-1,1]²
    newx = 2 * (location[1] - xmin) / (xmax - xmin) - 1
    newy = 2 * (location[2] - ymin) / (ymax - ymin) - 1

    return lagrange_eval(view(fl,:,:,e), newx, newy, rx, ry, ωx, ωy)
end