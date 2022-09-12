# IonoParamTest.jl
# Todd Anderson
# 12 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 6. Interpreting h' and β' (and f, and ν)
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/interpretinghpbeta/
#
# Explore the influence of h', β, electron collision frequency (ν), and transmitter 
# frequency (f) on amplitude curves.

using Printf, Plots, Interpolations
using LongwaveModePropagator
const LMP = LongwaveModePropagator

# Constant values
const BFIELD = BField(50e-6, π/2, 0)
const GND = GROUND[5]

const TX = Transmitter(20e3)
const RX = GroundSampler(0:5e3:3000e3, Fields.Ez)

## Varying h'
# Set f = 20 kHz, β = 0.4 km⁻¹; vary h'

function varyhp(hprimes)
    amps = Vector{Vector{Float64}}(undef, length(hprimes))
    for i in eachindex(hprimes)
        electrons = Species(LMP.QE, LMP.ME, z->waitprofile(z, hprimes[i], 0.4), 
            electroncollisionfrequency)
        waveguide = HomogeneousWaveguide(BFIELD, electrons, GND)
        E, a, p = propagate(waveguide, TX, RX)
        amps[i] = a
    end
    return amps
end

hprimes = 72:1:78
amps = varyhp(hprimes)

p = plot();
function buildplots!(p, amps)
    cmap = palette(:amp, length(hprimes)+1) # +1 allows a darker lightest color

    for i in eachindex(hprimes)
        plot!(p, RX.distance/1000, amps[i];
                label=hprimes[i], color=cmap[i+1]);
    end
end
buildplots!(p, amps);
plot!(p; size=(600,400), ylims=(22, 95),
        xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="h' (km)")

