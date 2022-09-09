# ModeTest.jl
# Todd Anderson
# 8 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 2. Mode finding

using Plots
using Plots.Measures
using RootsAndPoles

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, solvemodalequation, trianglemesh

# define frequencies: 10 kHz, 20 kHz, 100 kHz
lowfrequency = Frequency(10e3)
midfrequency = Frequency(20e3)
highfrequency = Frequency(100e3)

# define day and night ionospheres
day = Species(QE, ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency)
night = Species(QE, ME, z->waitprofile(z, 85, 0.9), electroncollisionfrequency)

# waveguides
daywaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), day, GROUND[5])
nightwaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), night, GROUND[5])

# mode equations
day_low_me = PhysicalModeEquation(lowfrequency, daywaveguide)
day_mid_me = PhysicalModeEquation(midfrequency, daywaveguide)
day_high_me = PhysicalModeEquation(highfrequency, daywaveguide)
night_mid_me = PhysicalModeEquation(midfrequency, nightwaveguide)

# plot titles
day_low_title = "10 kHz\nh': 75 km, β: 0.35 km⁻¹"
day_mid_title = "20 kHz\nh': 75 km, β: 0.35 km⁻¹"
day_high_title = "100 kHz\nh': 75 km, β: 0.35 km⁻¹"
night_mid_title = "20 kHz\nh': 85 km, β: 0.9 km⁻¹"

# define rectangular mesh
Δr = 0.2
x = 0:Δr:90
y = -40:Δr:0
mesh = x .+ im*y';

function modeequationphase(me, mesh)
    phase = Vector{Float64}(undef, length(mesh))
    Threads.@threads for i in eachindex(mesh)
        f = solvemodalequation(deg2rad(mesh[i]), me)
        phase[i] = rad2deg(angle(f))
    end
    return phase
end

phase = modeequationphase(day_mid_me, mesh);

heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(0, 90), ylims=(-40, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=day_mid_title,right_margin=2mm)
