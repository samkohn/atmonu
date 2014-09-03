# Unit conventions:
# SI + MeV
# use "natural" units for values, physical units for names.

# Some constants
const gevtomev = 1000
const hbar = 6.582119e-22
const c = 299792458.
const hbarc = hbar * c
const mec2 = 0.510998928
const mmuc2 = 106
const mnpc2 = 1.29333217
const f = 1.6857
const taun = 880.

# Parse 
function csv2tuple(filename)
    open(filename) do f
        lines = readlines(f)
        result = similar(lines, (Float64,Float64))
        for i in 1:length(lines)
            line = lines[i]
            comma = search(line, ',')
            x = line[1:comma-1]
            y = line[comma + 2:end-1]
            x = parsefloat(x)
            y = parsefloat(y)
            result[i] = (x,y)
        end
        result
    end
end

# Calculate flux from flux * E^2
function getflux(filename)
    data = csv2tuple(filename)
    for i in 1:length(data)
        x, y = data[i]
        y /= x * x
        y /= gevtomev
        x *= gevtomev
        data[i] = x, y
    end
    data
end

function getcrosssection(filename)
    data = csv2tuple(filename)
    for i in 1:length(data)
        x, y = data[i]
        y *= 1e-43
        x *= gevtomev
        data[i] = x, y
    end
    data
end

function interactionrate()
    fluxfile = "fluxe2.csv"
    crosssectionfile = "totalcrosssection_ccqe.csv"
    fluxes = getflux(fluxfile)
    css = getcrosssection(crosssectionfile)
    N = 3e31
    W = 0.0
    for i in 1:length(fluxes)-1
        E, flux = fluxes[i]
        if E > 1000 # 1GeV
            break
        end
        nextE, nextflux = fluxes[i+1]
        dE = nextE - E
        #println("$E, $nextE, $dE")
        cs = 0
        for j in 1:length(css)
            other_energy, cs = css[j]
            #println("$E, $other_energy")
            if other_energy > E
                break
            end
        end
        dW = flux * dE * 4pi * N * cs
        println("$dW, $E, $dE, $cs")
        W += dW
    end
    W
end


