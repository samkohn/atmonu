# Unit conventions:
# SI - kg + MeV

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

# Calculate the cross section for an inverse beta decay, in m^2
function crossSection(Enu)
    prefactor_num = 2pi^2 * hbarc^3
    prefactor_den = mmuc2^5 * f * taun
    prefactor = prefactor_num/prefactor_den
    result = prefactor * electronEnergy(Enu) * electronMomentum(Enu)
end

# Calculate the energy of an electron given the neutrino energy
function electronEnergy(Enu)
    Enu - mnpc2
end

# Calculate the momentum of an electron given the neutrino energy
function electronMomentum(Enu)
    sqrt(electronEnergy(Enu)^2 - mmuc2 * mmuc2)/c
end

# Calculate flux from flux * e^2
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


