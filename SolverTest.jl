# SolverTest.jl
# Todd Anderson
# 9 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 3. Solvers for ionosphere reflection coefficient
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/integratedreflection/
#
# LongwaveModePropagator uses the technique presented by Budden (1955) to compute the
# reflection coefficient of a horizontally stratified ionosphere consisting of an 
# anisotropic, collisional plasma. https://doi.org/10.1098/rspa.1955.0027
#
# Budden derives differential equation for the reflection coefficient R̲:
#
#   (2i/k)*(dR̲/dz) = W₂₁ + W₂₂R̲ - R̲W₁₁ - R̲W₁₂R̲
#
# where:
#   W: 4×4 matrix divided into four 2×2 submatrices each containing components of the 
#       T̲ matrix
#   R̲: 2×2 matrix of reflection coefficients
#
# NOTE: there are two differences in this code from the example code:
#   1. in the ODEProblem definition, LMPParams() is not passed to LMP.dRdz:
#   example:    prob = ODEProblem{false}(LMP.dRdz, Rtop, (topheight, 0.0), (me, LMPParams()))
#   here:       prob = ODEProblem{false}(LMP.dRdz, Rtop, (topheight, 0.0), me)
#   I found that the following line (sol = solve(prob ...)) fails with the following error:
#       "type Tuple has no field ea"
#   It looks like me and LMPParams() are not being combined correctly.  
#   The code runs with this change, but I'm not able to compare R̲ with the intended result.
#
#   2. solverstrings is defined differently:
#   example:    solverstrings = replace.(string.(solvers), "OrdinaryDiffEq."=>"")
#   here:       solverstrings = first.(split.(string.(solvers), "("))
#   The example version results in auto-filling of solver arguments, in the strings, resulting
#   in a very long string for Tsit5().  My version removes the substring beginning with "("
#   from each solverstring, so no arguments are displayed.  This fixes the following plot 
#   error: 
#       "GKS: Rectangle definition is invalid in routine SET_VIEWPORT"
#   and results in the heatmaps plotting as desired.

using Statistics
using Plots
using OrdinaryDiffEq
using Interpolations

using LongwaveModePropagator
using LongwaveModePropagator: StaticArrays, QE, ME
const LMP = LongwaveModePropagator

## R(z): ionosphere reflection coefficients as a function of height
species = Species(QE, ME, z->waitprofile(z, 75, 0.32), electroncollisionfrequency)
frequency = Frequency(24e3)
bfield = BField(50e-6, π/2, 0)
ea = EigenAngle(deg2rad(75));# wave angle of incidence on ionosphere

# start integration at a height significantly above reflection region, e.g. 110 km
topheight = 110e3
Mtop = LMP.susceptibility(topheight, frequency, bfield, species)
Rtop = LMP.bookerreflection(ea, Mtop)

# use OrdinaryDiffEq to integrate LMP.dRdz
ground = GROUND[1]
waveguide = HomogeneousWaveguide(bfield, species, ground)
me = PhysicalModeEquation(ea, frequency, waveguide);

# code as written in example errors on sol with error "type Tuple has no field ea"
#prob = ODEProblem{false}(LMP.dRdz, Rtop, (topheight, 0.0), (me, LMPParams()))
# replace (me, LMPParams()) with me; check results against example
prob = ODEProblem{false}(LMP.dRdz, Rtop, (topheight, 0.0), me)

sol = solve(prob, RK4(); abstol=1e-9, reltol=1e-9);
nothing #hide

# plot reflection coefficients with electron density and collision frequency curves.
# Reflection coefficients are complex, so plot their magnitude
zs = topheight:-1000:0

ne = species.numberdensity.(zs)
nu = species.collisionfrequency.(zs)
Wr = LMP.waitsparameter.(zs, (frequency,), (bfield,), (species,))

altinterp = LinearInterpolation(reverse(Wr), reverse(zs))
eqz = altinterp(frequency.ω)    # altitude where ω = ωᵣ

ne[end] = NaN # otherwise Plots errors
Wr[end] = NaN

p1 = plot([ne nu Wr], zs/1000;
        xlims=(10, 10^10), xaxis=(scale=:log10),
        ylabel="Altitude (km)",
        labels=["Nₑ (m⁻³)" "ν (s⁻¹)" "ωᵣ = ωₚ²/ν"], legend=:topleft,
        linewidth=1.5);

vline!(p1, [frequency.ω]; linestyle=:dash, color="gray", label="");
hline!(p1, [eqz/1000]; linestyle=:dash, color="gray", label="");
annotate!(p1, frequency.ω, 10, text(" ω", :left, 9));
annotate!(p1, 70, eqz/1000-3, text("ωᵣ = ω", :left, 9));

R11 = abs.(sol(zs; idxs=1))
R21 = abs.(sol(zs; idxs=2))
R12 = abs.(sol(zs; idxs=3))
R22 = abs.(sol(zs; idxs=4))

p2 = plot([R11 R21 R12 R22], zs/1000;
        xlims=(0, 1),
        yaxis=false, yformatter=_->"",
        legend=:right, labels=["R₁₁" "R₂₁" "R₁₂" "R₂₂"],
        linewidth=1.5);

hline!(p2, [eqz/1000]; linestyle=:dash, color="gray", label="");

plot(p1, p2; layout=(1,2), size=(800,400))
nothing # hide

## Generate random scenarios

function generatescenarios(N)
    eas = EigenAngle.(complex.(rand(N)*(π/2-π/6) .+ π/6, rand(N)*deg2rad(-10)))
    frequencies = Frequency.(rand(N)*50e3 .+ 10e3)

    B = rand(30e-6:5e-7:60e-6, N)
    # avoiding within 1ᵒ from 0ᵒ dip angle
    bfields = BField.(B, rand(N)*(π/2-0.018) .+ 0.018, rand(N)*2π)

    hps = rand(N)*20 .+ 69 #nice
    betas = rand(N)*0.8 .+ 0.2

    scenarios = Vector{PhysicalModeEquation}(undef, N)
    for i = 1:N
        species = Species(QE, ME, z-> waitprofile(z, hps[i], betas[i]),
                            electroncollisionfrequency)
        ground = GROUND[5]
        waveguide = HomogeneousWaveguide(bfields[i], species, ground)

        me = PhysicalModeEquation(eas[i], frequencies[i], waveguide)
        scenarios[i] = me
    end
    
    return scenarios
end

scenarios = generatescenarios(30);

## Reference solutions

ip = IntegrationParams(tolerance=1e-14, solver=RK4(), maxiters=1_000_000)
params = LMPParams(integrationparams=ip) 

Rrefs = [LMP.integratedreflection(scenario; params=params) for scenario in scenarios];

## Evaluate solvers
# compute and time results of different solver methods with a range of tolerances

function compute(scenarios, tolerances, solvers)
    dims = length(scenarios), length(tolerances), length(solvers)
    Rs = Array{StaticArrays.SMatrix{2,2,ComplexF64,4}}(undef, dims...)
    times = Array{Float64}(undef, dims...)

    for k in eachindex(solvers)
        for j in eachindex(tolerances)
            ip = IntegrationParams(tolerance=tolerances[j], solver = solvers[k])
            params = LMPParams(integrationparams=ip)

            for i in eachindex(scenarios)
                # warmup
                R = LMP.integratedreflection(scenarios[i]; params=params)

                # loop for average time
                N = 25
                t0 = time_ns()
                for n = 1:N
                    R = LMP.integratedreflection(scenarios[i]; params=params)
                end
                ttotal = time_ns() - t0
                
                Rs[i,j,k] = R
                times[i,j,k] = ttotal/N
            end
        end
    end
    return Rs, times
end

tolerances = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
tolerancestrings = string.(tolerances)

solvers = [RK4(), Tsit5(), BS5(), OwrenZen5(), Vern6(), Vern7(), Vern8()]
solverstrings = first.(split.(string.(solvers), "("))

Rs, times = compute(scenarios, tolerances, solvers);

# error in reflection coefficient matrices = maximum absolute difference of the four elements
# of the matrix compared to the reference reflection coefficient matrix

function differr(a, ref)
    differences = a .- ref
    aerror = similar(a, Float64)
    for i in eachindex(differences)
        absdiff = abs.(differences[i])
        aerror[i] = maximum(absdiff)
    end
    return aerror
end

Rerrs = differr(Rs, Rrefs)
mean_Rerrs = dropdims(mean(Rerrs; dims=1); dims=1)

heatmap(tolerancestrings, solverstrings, permutedims(log10.(mean_Rerrs));
        clims=(-9,-2),
        xlabel="tolerance", ylabel="solver",
        colorbar_title="log₁₀ max abs difference", colorbar=true)

# average runtimes
mean_times = dropdims(mean(times; dims=1); dims=1)

heatmap(tolerancestrings, solverstrings, permutedims(mean_times)/1e6;
        clims=(0,5),
        xlabel="tolerance", ylabel="solver",
        colorbar_title="time (μs)", colorbar=true)

