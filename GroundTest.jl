# MultiSpeciesTest.jl
# Todd Anderson
# 12 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 8. Ground
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/ground/
#
# Compare the effect of ground conductivity on amplitude along the ground in the waveguide 
# for day and nighttime ionospheres

using Plots, Printf
using LongwaveModePropagator
const LMP = LongwaveModePropagator

# LWPC standard ground indices:
sort(GROUND)

# define constants
const BFIELD = BField(50e-6, π/2, 0)

const TX = Transmitter(20e3)
const RX = GroundSampler(0:5e3:3000e3, Fields.Ez)

const DAY = Species(LMP.QE, LMP.ME, z->waitprofile(z, 75, 0.3), electroncollisionfrequency)
const NIGHT = Species(LMP.QE, LMP.ME, z->waitprofile(z, 82, 0.6), electroncollisionfrequency)

# define function that takes day or night species and returns electric field amplitude
function varyground(prf)
    amps = Vector{Vector{Float64}}(undef, length(GROUND))
    for i = 1:length(GROUND)
        wg = HomogeneousWaveguide(BFIELD, prf, GROUND[i])
        E, a, p = propagate(wg, TX, RX)
        amps[i] = a
    end
    return amps
end

# define function to plot amplitude curves
function buildplots!(p, amps)
    cmap = palette(:thermal, length(GROUND)+1) # +1 so we don't get to yellow!

    for i = 1:length(GROUND)
        plot!(p, RX.distance/1000, amps[i];
            label=@sprintf("%d, %.1g", GROUND[i].ϵᵣ, GROUND[i].σ), color=cmap[i]);
    end
end

# daytime results
amps = varyground(DAY)

p = plot();
buildplots!(p, amps);
plot!(p; size=(600,400), ylims=(0,95), title="Day", legend=(0.85, 1.02),
    xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="ϵᵣ, σ")

# nighttime results
amps = varyground(NIGHT)

p = plot();
buildplots!(p, amps);
plot!(p; size=(600,400), ylims=(0,95), title="Night", legend=(0.85, 1.02),
    xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="ϵᵣ, σ")