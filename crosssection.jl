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

function interactionrate(fluxfile, crosssectionfile, num_targets)
    fluxes = getflux(fluxfile)
    css = getcrosssection(crosssectionfile)
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
        dW = flux * dE * 4pi * num_targets * cs
        #println("$dW, $E, $dE, $cs")
        W += dW
    end
    W
end

end

using CrossSection
using ArgParse

settings = ArgParseSettings()
@add_arg_table settings begin
    "--flux-file", "-f"
        help = "a CSV file with E [GeV], Phi*E^2 [m^-2 s^-1 sr^-1 GeV] pairs"
        arg_type = String
        default = "fluxe2.csv"
    "--cross-section-file", "-c"
        help = "a CSV file with Sigma [10^-39 cm^2], E [GeV] pairs"
        arg_type = String
        default = "totalcrosssection_ccqe.csv"
    "--num-targets", "-n"
        help = "the number of target particles present"
        arg_type = Float64
        default = 3e31
end
args = parse_args(settings)
rate = interactionrate(args["flux-file"], args["cross-section-file"], args["num-targets"])
println("$rate s^-1")
