module CrossSection

using DataFrames

export getflux, getcrosssection, interactionrate, divideBySquare

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


# Calculate flux from flux * E^2
# Returns a DataFrame with columns :energy and :flux
function getflux(filename)
    data = readtable(filename, names=[:energy, :fluxe2])
    nrows = size(data, 1)
    data[:flux] = zeros(Float64, nrows)
    for i in 1:nrows
        # retrieve the data
        energy = data[:energy][i] # GeV
        fluxe2 = data[:fluxe2][i] # m^-2 s^-1 sr^-1 GeV
        # calculate the flux and convert to appropriate units
        data[:flux][i] = divideBySquare(fluxe2, energy) / gevtomev # m^-2 s^-1 sr^-1 MeV^-1
        # convert the energy to appropriate units
        data[:energy][i] *= gevtomev # MeV
    end
    # only return the energy and the flux (not flux * E^2)
    data[:,[:energy, :flux]]
end

function divideBySquare(a, b)
    a/(b * b)
end

function getcrosssection(filename)
    sigma = :σ
    energy = :energy
    data = readtable(filename, names=[energy, sigma])
    nrows = size(data, 1)
    for i in 1:nrows
        data[energy][i] *= gevtomev
        data[sigma][i] *= 1e-43
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
    for i in 1:size(fluxes, 1) - 1
        E = fluxes[:energy][i]
        flux = fluxes[:flux][i]
        if E > 1000 # 1GeV
            break
        end
        nextE = fluxes[:energy][i+1]
        nextflux = fluxes[:flux][i+1]
        dE = nextE - E
        #println("$E, $nextE, $dE")
        cs = 0
        for j in 1:size(css, 1)
            other_energy = css[:energy][j]
            cs = css[:σ][j]
            if other_energy > E
                break
            end
        end
        dW = flux * dE * 4pi * N * cs
        #println("$dW, $E, $dE, $cs")
        W += dW
    end
    W
end

end
