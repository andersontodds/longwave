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

## Varying β
# fix h' at 78 km, vary β

function varybeta(betas)
    amps = Vector{Vector{Float64}}(undef, length(betas))
    for i in eachindex(betas)
        electrons = Species(LMP.QE, LMP.ME, z->waitprofile(z, 78, betas[i]), 
            electroncollisionfrequency)
        waveguide = HomogeneousWaveguide(BFIELD, electrons, GND)
        E, a, p = propagate(waveguide, TX, RX)
        amps[i] = a
    end
    return amps   
end

betas = [0.3, 0.4, 0.5, 0.7, 0.9, 2.0]
amps = varybeta(betas)

p = plot();
function buildplots!(p, amps)
    cmap = palette(:algae, length(betas)+1)

    for i in eachindex(betas)
        plot!(p, RX.distance/1000, amps[i];
            label=betas[i], color=cmap[i+1]);
    end
end
buildplots!(p, amps);
plot!(p; size=(600,400), ylims=(22, 95), 
    xlabel="Range (km", ylabel="Amplitude (dB)", legendtitle="β (km⁻¹)")

## Varying frequency
# set h' and β at typical daytime Wait ionosphere values, vary frequency between
# 5 kHz and 50 kHz.

function varyfreq(freqs)
    electrons = Species(LMP.QE, LMP.ME, z->waitprofile(z, 75, 0.35), 
        electroncollisionfrequency)
    waveguide = HomogeneousWaveguide(BFIELD, electrons, GND)

    amps = Vector{Vector{Float64}}(undef, length(freqs))
    for i in eachindex(freqs)
        tx = Transmitter(freqs[i])
        E, a, p = propagate(waveguide, tx, RX)
        amps[i] = a
    end
    return amps
end

freqs = 5e3:5e3:50e3
amps = varyfreq(freqs)

p = plot();
function buildplots!(p, amps)
    cmap = palette(:oslo, rev=true, length(freqs)+1)
    
    for i in eachindex(freqs)
        fkHz = trunc(Int, freqs[i]/1000)
        λkm = trunc(LMP.C0/freqs[i]/1000; digits=1)
        plot!(p, RX.distance/1000, amps[i] .+ (10*i);
            label=string(fkHz, ", ", λkm), color=cmap[i+1]);
    end
end
buildplots!(p, amps);
plot!(p; size=(600,400), 
    xlabel="Range, (km)", ylabel="Amplitude (dB)", legendtitle="f kHz, λ km")
    
## Collision frequency
# So far we have used the electron collision frequency from Morfitt and Shellman (1976). 
# Instead, let's look at different electron collision frequency profiles.

# Plot electron collision frequency profiles
collisionfrequency(z, ν₀, a) = ν₀*exp(-a*z/1000) # 1/1000 converts z to km

function reflectionheight(params, hp, β)
    species = Species(LMP.QE, LMP.ME, z->waitprofile(z, hp, β), 
            z->collisionfrequency(z, params...))
    ωr = LMP.waitsparameter.(alt, (TX.frequency, ), (BFIELD,), (species,))
    itp = LinearInterpolation(ωr, alt)
    eqz = itp(TX.frequency.ω)
    
    return eqz
end

alt = 40e3:500:110e3
params = [(1.816e11, 0.18),
          (1.816e11/5, 0.15),
          (1.816e11, 0.15),
          (1.816e11*5, 0.15),
          (1.816e11, 0.12)]

p = plot();
function buildplots!(p, params, hp, β)
    cmap = palette(:rainbow, length(params))

    for i in eachindex(params)
        eqz = reflectionheight(params[i], hp, β)

        ν₀label = params[i][1]/1.816
        plot!(p, collisionfrequency.(alt, params[i]...), alt/1000;
            label=@sprintf("%.0e, %.2f", ν₀label, params[i][2]), color=cmap[i]);
        hline!(p, [eqz/1000]; linestyle=:dash, color=cmap[i], linewidth=0.6, label=nothing);
    end
end
buildplots!(p, params, 75, 0.35)
annotate!(p, [(1e3, 64, text("ωᵣ = ω", 10))]);
plot!(p; size=(600,400), xscale=:log10, 
    xlabel="ν (s⁻¹)", ylabel="Altitude (km)", legendtitle="1.816 ν₀, a")

# typical daytime electron density profile
function varycollisionfreq(params, hp, β)
    amps = Vector{Vector{Float64}}(undef, length(params))
    for i in eachindex(params)
        electrons = Species(LMP.QE, LMP.ME, z->waitprofile(z, hp, β), 
            z->collisionfrequency(z, params[i]...))
        waveguide = HomogeneousWaveguide(BFIELD, electrons, GND)
        E, a, p = propagate(waveguide, TX, RX)
        amps[i] = a
    end
    return amps
end
amps = varycollisionfreq(params, 75, 0.35)

p = plot();
function buildplots!(p, amps)
    cmap = palette(:rainbow, length(params))

    for i in eachindex(params)
        ν₀label = params[i][1]/1.816
        plot!(p, RX.distance/1000, amps[i];
            label=@sprintf("%.0e, %.2f", ν₀label, params[i][2]), color=cmap[i]);
    end
end
buildplots!(p, amps);
plot!(p; size=(600, 400), ylims=(22, 95),
    xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="1.816 ν₀, a")

# typical nighttime electron density profile
amps = varycollisionfreq(params, 82, 0.55)

p = plot();
function buildplots!(p, amps)
    cmap = palette(:rainbow, length(params))

    for i in eachindex(params)
        v₀label = params[i][1]/1.816
        plot!(p, RX.distance/1000, amps[i];
              label=@sprintf("%.0e, %.2f", v₀label, params[i][2]), color=cmap[i]);
    end
end
buildplots!(p, amps);
plot!(p; size=(600,400), ylims=(22, 95),
      xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="1.816 ν₀, a")