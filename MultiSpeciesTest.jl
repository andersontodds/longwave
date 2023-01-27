# MultiSpeciesTest.jl
# Todd Anderson
# 12 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 7. Multiple ionospheric species
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/multiplespecies/
#
# Look at the influence of massive ions on wave properties

using Plots
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME

QE
ME

# define electron and ion density profiles.  Use a Wait profile for electrons, somewhat
# unrealistic exponential profile for ions
Ne(z) = waitprofile(z, 75, 0.32)
Np(z) = 2e6*exp(1e-4*z)

Nn(z) = Np(z) - Ne(z)

# for plotting, replace number densities below 1 m⁻³ with NaN
mask(x) = x < 1 ? NaN : x

z = 0:1e3:110e3
p = plot(mask.(Ne.(z)), z/1000; label="Ne", linewidth=1.5);
plot!(p, mask.(Np.(z)), z/1000; label="Np", linewidth=1.5);
plot!(p, mask.(Nn.(z)), z/1000; label="Nn", linewidth=1.5);
plot!(p, xlabel="Density (m⁻³)", ylabel="Altitude (km)",
    xscale=:log10, legend=:topleft)

p = plot(electroncollisionfrequency.(z), z/1000, label="νₑ", linewidth=1.5);
plot!(p, ioncollisionfrequency.(z), z/1000, label="νᵢ", linewidth=1.5);
plot!(p, xlabel="Collision frequency (s⁻¹)", ylabel="Altitude (km)", xscale=:log10)

# Species: LWPC allows electrons, and one species each of positive and negative ions.
# Each ion species has a mass of 58,000*mₑ, approximately the mass of O₂. 
electrons = Species(QE, ME, Ne, electroncollisionfrequency)
posions = Species(abs(QE), 58000*ME, Np, ioncollisionfrequency)
negions = Species(QE, 58000*ME, Nn, ioncollisionfrequency)

# Compare an electrons-only and electrons+ions ionosphere
tx = Transmitter(24e3)
rx = GroundSampler(0:5e3:2000e3, Fields.Ez)

bfield = BField(50e-6, π/2, 0)
ground = GROUND[10]

ewg = HomogeneousWaveguide(bfield, electrons, ground)
eiwg = HomogeneousWaveguide(bfield, (electrons, posions, negions), ground)

Ee, ae, pe = propagate(ewg, tx, rx);
Eei, aei, pei = propagate(eiwg, tx, rx);

p1 = plot(rx.distance/1000, ae; 
    label="elecrons", ylabel="Amplitude (dB μV/m)", xaxis=false);
plot!(p1, rx.distance/1000, aei, label="electrons & ions");
p2 = plot(rx.distance/1000, aei-ae;
    ylims=(-0.5, 0.5), xlabel="Range (km)", ylabel="Δ", legend=false);
plot(p1, p2, layout=grid(2,1,heights=[0.7, 0.3]))
nothing #hide



