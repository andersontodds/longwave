# test DifferentialEquations.jl deprecation warning
# example from: https://docs.sciml.ai/DiffEqDocs/stable/examples/classical_physics/

using DifferentialEquations
using BenchmarkTools
using Suppressor
using CairoMakie

#Half-life of Carbon-14 is 5,730 years.
C₁ = 5.730

#Setup
u₀ = 1.0
tspan = (0.0, 1.0)

#Define the problem
radioactivedecay(u,p,t) = -C₁*u

#Pass to solver
prob = ODEProblem(radioactivedecay,u₀,tspan)
@time sol = solve(prob,Tsit5());
# throws warning with following versions:
# OrdinaryDiffEq v6.27.2
# DifferentialEquations v7.6.0


## plot simulated sferics from Dowden 2002

r = 10000e3
c = 2.99792458e8    # ms⁻¹
v_g = 0.9905*c      # speed of light in the EIWG

r/v_g

tᵣ = -0.2e-3:1e-6:1e-3;
t = tᵣ .+ r/c;
ωf = 2*pi*(2e3:2e2:24e3);
# ωf = 2*pi*7e3
ωₐ = 2*pi*14e3;     # frequency of peak spectral density: ~12 kHz
ωᵣ = 2*pi*11e3;     # tune this
ω₀ = 2*pi*1.6e3;    # waveguide cutoff frequency

A = (cos.(pi*(ωf.-ωₐ)./(2*ωᵣ))).^2
waveform = zeros(length(t));
for ω in eachindex(ωf)
    component = A[ω].*cos.(ωf[ω].*(t.-(r/c)*(1 - (ω₀^2)/(ωf[ω]^2))^(1/2)))
    waveform = waveform + component;
end

begin fig = Figure()
    fa = Axis(fig[1,1];
        xlabel="t - r/c (s)",
        ylabel="E (arbitrary units)",
        title="dispersion")
    # xlims!(fa, [minimum(tᵣ) maximum(tᵣ)])
    xlims!(fa, [-0.2e-3 1e-3])
    ylims!(fa, [-100 100])
    lines!(fa, tᵣ, waveform)
    fig
end