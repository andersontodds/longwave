# WavefieldTest.jl
# Todd Anderson
# 9 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 4. Wavefield integration
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/wavefieldintegration/

using CSV, Interpolations
using Plots
using Plots.Measures
using OrdinaryDiffEq

using TimerOutputs
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
const LMP = LongwaveModePropagator

## The ionosphere
# see model ionosphere in Piggot et al. 1965 (https://doi.org/10.1098/rsta.1965.0005)
#   f = 16 kHz
#   angle of incidence: 40 degrees from normal
#   B field magnitude = 50,000 nT
#           dip angle = 68 degrees
#           azimuth   = 111 degrees
#   ground conductivity index: 8 (see dictionary GROUND)
ea = EigenAngle(deg2rad(40))
frequency = Frequency(16e3)
bfield = BField(50e-6, deg2rad(68), deg2rad(111))
ground = GROUND[8]

# read in digitized version of Piggott electron density and collision frequency curves
root_dir = dirname(dirname(pathof(LMP)))
examples_dir = joinpath(root_dir, "examples")
data = CSV.File(joinpath(examples_dir, "piggott1965_data.csv"))

# interpolation object
day_itp = interpolate((data.day_ne_z,), data.day_ne, Gridded(Linear()))
day_etp = extrapolate(day_itp, Line())
dayfcn(z) = (v = day_etp(z); v > 0 ? v : 0.001)

# night has some rows of 'missing' data
filtnight_z = collect(skipmissing(data.night_ne_z))
filtnight = collect(skipmissing(data.night_ne))
night_itp = interpolate((filtnight_z,), filtnight, Gridded(Linear()))
night_etp = extrapolate(night_itp, Line())
nightfcn(z) = (v = night_etp(z); v > 0 ? v : 0.001)

# so does the collision frequency
filtnu_z = collect(skipmissing(data.nu_z))
filtnu = collect(skipmissing(data.nu))
nu_itp = interpolate((filtnu_z,), filtnu, Gridded(Linear()))
nu_etp = extrapolate(nu_itp, Line())
nufcn(z) = (v = nu_etp(z); v > 0 ? v : 0.001)

day = Species(QE, ME, dayfcn, nufcn)
night = Species(QE, ME, nightfcn, nufcn)

## Scaled, integrated wavefields 
# for a vertically-stratified, horizontally-homogeneous ionosphere, we only need to study the 
# horizontal components of the E and H fields.  Where x is in the propagation direction, z is
# vertically upwards, and x̂ × ŷ = ẑ, E_z and H_z are eliminated, and the differential 
# equations for a wave in the ionosphere can be written in matrix form as:
#
#   de/dz = -ikT̲e
#
# where
#       (  E_x)
#   e = ( -E_y)
#       (  H_x)
#       (  H_y)
#
# and T̲ is the 4×4 matrix consisting of elements of the susceptibility matrix for the
# ionosphere (calculated with LongwaveModePropagator.tmatrix).  
# See Clemmow and Heading (1954): https://doi.org/10.1017/S030500410002939X
#
# To compute wavefields, calculate them at some height z, use de/dz to step to a new height,
# and repeat from a point high in the ionosphere down to the ground.  The initial solution
# is from the Booker Quartic, which is a solution for the wavefields in a homogeneous 
# ionosphere.  At a great heigh in the ionosphere, we are interested in the two quartic
# solutions corresponding to upward-going waves.  Therefore, we are integrating two sets of 4
# complex variables simultaneously.  Throughout the integration, the two solutions lose their
# independence because of numerical accuracy limitations over a wide range of field magnitudes.
# To maintain accuracy, the wavefields are orthonormalized repeatedly during the downward 
# integration, and the scaling values are stored so that the fields can be "recorrected" after
# the integration is complete.

# compute the real components of Ex,1 and Hx,2 wavefields with and without recorrection
zs = 110e3:-50:0
zskm = zs/1000

e = LMP.integratewavefields(zs, ea, frequency, bfield, day; unscale=true)
e_unscaled = LMP.integratewavefields(zs, ea, frequency, bfield, day; unscale=false)

ex1 = getindex.(e,1)
ex1_unscaled = getindex.(e_unscaled, 1)
hx2 = getindex.(e,7)
hx2_unscaled = getindex.(e_unscaled, 7)

p1 = plot(real(ex1), zskm; title="\$E_{x,1}\$",
        ylims=(0, 90), xlims=(-1.2, 1.2), linewidth=1.5, ylabel="altitude (km)",
        label="corrected", legend=:topleft);
plot!(p1, real(ex1_unscaled), zskm;
        linewidth=1.5, label="scaled only");

p2 = plot(real(hx2), zskm; title="\$H_{x,2}\$",
        ylims=(0, 90), xlims=(-1.2, 1.2), linewidth=1.5, ylabel="altitude (km)",
        label="corrected", legend=:false);
plot!(p2, real(hx2_unscaled), zskm;
        linewidth=1.5, label="scaled only");

plot(p1, p2, layout=(1,2))
nothing #hide
# NOTE: looks like my hx2_unscaled differs slightly from the example; the lowest-altitude
# discontinuity occurs at a lower altitude (<60 km) than in the example (>60km).

## Differential Equations solver
# test different solvers to determine which makes the fewest function calls while maintaing
# good accuracy

# pass IntegrationParams through the LMPParams struct
TO = TimerOutput()

zs = 110e3:-50:0
solvers = [Tsit5(), BS5(lazy=false), OwrenZen5(),
           Vern6(lazy=false), Vern7(lazy=false), Vern8(lazy=false), Vern9(lazy=false)]

#solverstrings = replace.(string.(solvers), "OrdinaryDiffEq."=>"")
solverstrings = string.(first.(split.(string.(solvers), "(")))

day_es = []
night_es = []
for s in eachindex(solvers)
    ip = IntegrationParams(solver=solvers[s], tolerance=1e-6)
    params = LMPParams(wavefieldintegrationparams=ip)

    # make sure method is compiled
    LMP.integratewavefields(zs, ea, frequency, bfield, day; params=params);
    LMP.integratewavefields(zs, ea, frequency, bfield, night; params=params);

    solverstring = solverstrings[s]
    let day_e, night_e
        # repeat 200 times to average calls
        for i = 1:200
            # day ionosphere
            @timeit TO solverstring begin
                day_e = LMP.integratewavefields(zs, ea, frequency, bfield, day; params=params)
            end
            # night ionosphere
            @timeit TO solverstring begin
                night_e = LMP.integratewavefields(zs, ea, frequency, bfield, night; params=params)
            end
        end
        push!(day_es, day_e)
        push!(night_es, night_e)
    end
end

day_e1s = [getindex.(e,1) for e in day_es]

plot(real(day_e1s), zs/1000;
    label=permutedims(solverstrings), legend=:topleft)