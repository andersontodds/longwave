# LMPTest.jl
# Todd Anderson
# 8 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using Plots

# vertical dipole transmitter at 24 kHz
tx = Transmitter(2.4e4)

# sample vertical electric field every 5 km out to 2000 km from tx
rx = GroundSampler(0:5e3:2000e3, Fields.Ez)

# vertical magnetic field
bfield = BField(50e-6, π/2, 0)

# daytime ionosphere
electrons = Species(QE, ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency)

# "typical" earth ground 
ground = Ground(10,1e-4)

waveguide = HomogeneousWaveguide(bfield, electrons, ground)

# return the complex electric field, amplitude, and phase
E, a, p = propagate(waveguide, tx, rx);

plot(rx.distance/1000, a, xlabel="distance (km)", ylabel="amplitude (dB)")
plot(rx.distance/1000, p, xlabel="distance (km)", ylabel="phase (radians)")
 
