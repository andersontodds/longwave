# LMPTest.jl
# Todd Anderson
# 8 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 1. Basic propagation simulation

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using Suppressor
using Plots

# vertical dipole transmitter at 24 kHz
tx = Transmitter(2.4e4)

# sample vertical electric field every 5 km out to 2000 km from tx
rx = GroundSampler(0:5e3:2000e3, Fields.Ez)

# vertical magnetic field
bfield = BField(50e-6, π/2, 0)

# daytime ionosphere
h = 75      # km
β = 0.35    # km⁻¹

electrons = Species(QE, ME, z->waitprofile(z, h, β), electroncollisionfrequency)

# "typical" earth ground 
ground = Ground(10,1e-4)

waveguide = HomogeneousWaveguide(bfield, electrons, ground)

# return the complex electric field, amplitude, and phase
@time @suppress E, a, p = propagate(waveguide, tx, rx);

plot(rx.distance/1000, a; 
    xlabel="distance (km)", ylabel="amplitude (dB)", 
    linewidth=1.5, legend=false)
plot(rx.distance/1000, rad2deg.(p); 
    xlabel="distance (km)", ylabel="phase (deg)", 
    linewidth=1.5, legend=false)
 
## SegmentedWaveguide example with 2 segments
h1 = 75     # km
β1 = 0.35   # km⁻¹
h2 = 82     # km
β2 = 0.5    # km⁻¹

distances = [0.0, 1000e3]
species = [ Species(QE, ME, z->waitprofile(z, h1, β1), electroncollisionfrequency), 
            Species(QE, ME, z->waitprofile(z, h2, β2), electroncollisionfrequency)]

waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield, species[i], ground, 
            distances[i]) for i in 1:2]);

E, a, p = propagate(waveguide, tx, rx);

plot(rx.distance/1000, a;
        xlabel="range (km)", ylabel="amplitude (dB)",
        linewidth=1.5, legend=false)

plot(rx.distance/1000, rad2deg.(p);
        xlabel="range (km)", ylabel="phase (∘)",
        linewidth=1.5, legend=false)
