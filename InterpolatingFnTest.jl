# MultiSpeciesTest.jl
# Todd Anderson
# 12 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 8. Multiple ionospheric species
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/interpolatingfunctions/
#
# Compare different interpolating functions for density and collision frequency

# Note: FaradayInternationalReferenceIonosphere.jl requires at least Julia v1.7.1!

using Printf, Statistics
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using Plots, Distances

# example code is "using FIRITools", but FIRITools.jl link points to
# FaradayInternationalReferenceIonosphere.jl repo.  Import this as FIRI.
using FaradayInternationalReferenceIonosphere
import FaradayInternationalReferenceIonosphere as FIRI

using Interpolations, NormalHermiteSplines

# discrete profile data
zs = 0:1e3:110e3
Ne = FIRI.extrapolate(firi(50,30), zs)

# construct interpolators
linear_itp = LinearInterpolation(zs, Ne)
cubic_itp = CubicSplineInterpolation(zs, Ne)

# monotonic interpolators from Interpolations.jl
fb_itp = Interpolations.interpolate(zs, Ne, FritschButlandMonotonicInterpolation())
fc_itp = Interpolations.interpolate(zs, Ne, FritschCarlsonMonotonicInterpolation())
s_itp = Interpolations.interpolate(zs, Ne, SteffenMonotonicInterpolation())

# Hermite splines 
spline = prepare(collect(zs), RK_H1())
spline = construct(spline, Ne)
hermite_itp(z) = evaluate_one(spline, z);

# note on spline kernels:
#   RK_H0 : continuous spline
#   RK_H1 : continuously-differentiable spline
#   RK_H2 : continuously-twice-differentiable spline
zs_fine = 40e3:100:110e3
Ne_fine = FIRI.extrapolate(firi(50,30), zs_fine);

linear_fine = linear_itp.(zs_fine);
cubic_fine = cubic_itp.(zs_fine);
fb_fine = fb_itp.(zs_fine);
fc_fine = fc_itp.(zs_fine);
s_fine = s_itp.(zs_fine);
hermite_fine = hermite_itp.(zs_fine);

# compare profiles using percentage difference relative to true profile Ne_fine
cmp(a,b) = (a-b)/b*100

dNe = cmp.(Ne_fine, Ne_fine);
dlinear = cmp.(linear_fine, Ne_fine);
dcubic = cmp.(cubic_fine, Ne_fine);
dfb = cmp.(fb_fine, Ne_fine);
dfc = cmp.(fc_fine, Ne_fine);
ds = cmp.(s_fine, Ne_fine);
dhermite = cmp.(hermite_fine, Ne_fine);

# plot number densities with log scale: set values less than 0.1 to NaN
cl(x) = replace(v->v <= 0.1 ? NaN : v, x)
lc(x) = replace(x, NaN => 0)

ne_c = cl([Ne_fine linear_fine cubic_fine hermite_fine fb_fine fc_fine s_fine])

p1 = plot(ne_c,
        zs_fine/1000; xscale=:log10, xlabel="Ne (m⁻³)", ylabel="Altitude (km)", 
        legend=:topleft, labels=["Truth" "Linear" "Cubic" "Hermite" "FritschButland" "FritschCarlson" "Steffen"]);
p2 = plot(lc([dNe dlinear dcubic dhermite dfb dfc ds]),
        zs_fine/1000; xlabel="% difference", legend=false);
plot!(p1, p2, layout=(1, 2), size=(800, 400), margin=3Plots.mm)
# not sure if this result agrees with the following text in the example:
# "Unsurprisingly, the error is highest at the cutoff altitude of 40 km where the
# densities below are 0."  The error appears to vary periodically with altitude, and
# is about the same between 40 and 60 km.  Unfortunately the example output plot is not 
# included online or in the examples directory.

# Compute the total absolute difference between each interpolator and the truth, and
# compute the average percentage difference for each
for (n, v) in ("linear"=>(linear_fine, dlinear), "cubic"=>(cubic_fine, dcubic), 
    "hermite"=>(hermite_fine, dhermite), "FritschButland"=>(fb_fine, dfb), 
    "FritschCarlson"=>(fc_fine, dfc), "Steffen"=>(s_fine, ds))
    @printf("%s:  %.3g  %.3g\n", n, cityblock(Ne_fine, v[1]), mean(abs, v[2]))
end

## Propagation results

interpolators = (
    "truth" => Interpolations.interpolate(0:100:110e3,
        FIRI.extrapolate(firi(50, 30), 0:100:110e3),
        FritschButlandMonotonicInterpolation()),
    "linear" => linear_itp,
    "cubic" => cubic_itp,
    "FritschButland" => fb_itp,
    "FritschCarlson" => fc_itp,
    "Steffen" => s_itp
)

function propagateitp(interpolators)
    bfield = BField(50e-6, deg2rad(68), deg2rad(111))
    ground = GROUND[5]

    tx = Transmitter(24e3)
    rx = GroundSampler(0:5e3:3000e3, Fields.Ez)

    results = Dict{String,Tuple{Float64,Vector{Float64},Vector{Float64}}}()
    for (n, itp) in interpolators
        species = Species(QE, ME, itp, electroncollisionfrequency)
        waveguide = HomogeneousWaveguide(bfield, species, ground)
        t0 = time()
        _, amp, phase = propagate(waveguide, tx, rx)
        runtime = time() - t0
        results[n] = (runtime, amp, phase)
    end
    return results
end

propagateitp(interpolators); # warmup
results = propagateitp(interpolators);
nothing #hide

d = 0:5:3000
p1 = plot(ylabel="Amplitude (dB μV/m)");
p2 = plot(xlabel="Range (km)", ylims=(-0.02, 0.02), ylabel="Δ", legend=false);
for (n, v) in results
    plot!(p1, d, v[2], label=n);
    plot!(p2, d, v[2]-results["truth"][2]);
end
plot(p1, p2, layout=grid(2,1,heights=[0.7,0.3]))
nothing #hide

# mean absolute amplitude difference between each technique and the truth profile amplitude 
for n in ("linear", "cubic", "FritschButland", "FritschCarlson", "Steffen")
    @printf("%s: %.3e\n",n, meanad(results[n][2], results["truth"][2]))
end

