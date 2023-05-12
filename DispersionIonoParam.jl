# May 11 2023
# test waveguide impacts on broadband dispersion
# see dispersion methods in DispersionTest.jl and 
# ionosphere parameter variation in IonoParamTest.jl
#
# Test state 1: typical nightside
#   HomogeneousWaveguide
#   r = 0:10e3:5000e3 m → 0-5000 km in 10km increments
#   h', β, ν :  typical nightside ionosphere
#   ϵ, σ:       typical ocean ground
#
# Test state 2: typical dayside
#   HomogeneousWaveguide
#   r = 0:10e3:5000e3 m → 0-5000 km in 10km increments
#   h', β, ν :  typical dayside ionosphere
#   ϵ, σ:       typical ocean ground
#
# Test state 3: sharp conducting boundaries
#   HomogeneousWaveguide
#   r = 0:10e3:5000e3 m → 0-5000 km in 10km increments
#   h', β, ν :  typical nightside height, large β, nominal ν (ν shouldn't be important for very large β)
#       or, write new sharpreflector function to return delta-function electron number density
#   ϵ, σ :      typical ocean ground
#
# Parameter variation:
#   Wait ionosphere:
#       h'      height
#       β       sharpness
#       ν       neutral collision frequency: non-parameterized
#   

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using LMPTools
using LsqFit
using CairoMakie
using Printf
# using Dates

# constants
const c = 2.99792458e8;    # m s⁻¹
const v_g = 0.9905*c;      # speed of light in the EIWG
const BFIELD = BField(50e-6, π/2, 0);
const GND = GROUND[10];
const Nₑ = 1.0e20;           # m⁻³, sharp reflector electron density

const proprange = 5000e3;
const gsrange = 0:10e3:proprange;
const rx = GroundSampler(gsrange, Fields.Ez);

const freqs = 6e3:1e3:18e3;
const ωfreqs = 2*pi*freqs;

# 0. define functions
function varyfreq(waveguide, rx, freqs)
    amps = Vector{Vector{Float64}}(undef, length(freqs))
    phases = Vector{Vector{Float64}}(undef, length(freqs))
    # precompute waveguide modes and provide to propagate call?
    # modes depend on frequency → not sure this can be sped up
    for i in eachindex(freqs)
        tx = Transmitter(freqs[i])
        E, a, p = propagate(waveguide, tx, rx);
        amps[i] = a # this breaks for single frequency, need another method with amps[i] = [a]
        phases[i] = p # this breaks for single frequency
    end
    return amps, phases
end

function shiftfit(xdata, ydata, thres)

    @. model(x, p) = p[1]*x + p[2] + p[3]*(1/x)
    # bounds
    lb = [-Inf, -Inf, 0];
    ub = [Inf, Inf, Inf];
    p0 = [1E-6, 0.1, thres];
    # set up loop conditions
    xshift = copy(xdata); # can remove these if preserving original data is not important
    yshift = copy(ydata);
    # xout = Vector{Float64}(undef, 0);
    # yout = Vector{Float64}(undef, 0);
    for i in 2:length(ydata)
        n = ceil((yshift[i] - yshift[i-1])*1e3/(xshift[i] - xshift[i-1])) # number of 2π shifts
        yshift[i] = yshift[i] - n*2π
    end
    fit = curve_fit(model, xshift, yshift, p0, lower=lb, upper=ub)
    # sigma = stderror(fit)[3]
    # while sigma > thres
    #     out = findmax(abs.(fit.resid));
    #     push!(xout, xin[out[2]]);
    #     push!(yout, yin[out[2]]);
    #     popat!(xin, out[2]);
    #     popat!(yin, out[2]);
    #     fit = curve_fit(model, xin, yin, p0, lower=lb, upper=ub)
    #     sigma = stderror(fit)[3] 
    # end
    fit, xshift, yshift
end

function sharpreflector(z, hᵣ, n)
    # define delta-function ionosphere reflector with: 
    #   nₑ = n above hᵣ
    #   nₑ = 0 below hᵣ
    # inputs:
    #   z : altitude in meters
    #   hᵣ: reflector altitude in kilometers
    #   n : electron density above reflector altitude

    if z < hᵣ
        nₑ = 1.0;
    else
        nₑ = n
    end

    return nₑ

end

function varyionoparam(h, β)
    # run varyfreq() for different sets of ionosphere parameters
    a3_r = zeros(length(h), length(β));
    for i in eachindex(h), j in eachindex(β)
        hi = h[i];
        βi = β[j];
        @printf "Running varyfreq for h' = %i, β = %.2g" hi βi
        if βi > 10
            species = Species(QE, ME, z->sharpreflector(z, hi, Nₑ), electroncollisionfrequency)
        else
            species = Species(QE, ME, z->waitprofile(z, hi, βi), electroncollisionfrequency)
        end
        waveguide = HomogeneousWaveguide(BFIELD, species, GND, 0.0);
        ~, phases = varyfreq(waveguide, rx, freqs)
        r = length(gsrange);
        finalphase = [phases[i][r] for i in eachindex(freqs)];
        finalphasefit, ~, ~ = shiftfit([ωfreqs;], finalphase, proprange/50);
        a3_r[i,j] = finalphasefit.param[3]/proprange    # range-normalized dispersion

    end
    return a3_r
end

# 1. run test states 
nightside   = [87, 0.5];
dayside     = [70, 0.3];
pplate      = [87, 100.0];

@time night_a3_r = varyionoparam(nightside...)
@time day_a3_r = varyionoparam(dayside...)
# @time pplate_a3_r = varyionoparam(pplate...)
# TODO: fix sharpreflector; currently breaks integration

# 2. vary parameters
h = 60:5:100;
β = 0.1:0.1:1.0;

@time a3_r = varyionoparam(h,β);

# 3. make plots 

begin fig = Figure(resolution = (1000,1000))
    # global plot properties
    fontsize_theme = Theme(fontsize=20)
    set_theme!(fontsize_theme)

    fa1 = Axis(fig[1,1], title="Wait ionosphere parameters",
        xlabel="β (km⁻¹)",
        ylabel="h' (km)")

    # add more axes if needed
    
    # heatmap
    hm = heatmap!(fa1, β, h, a3_r)
    Colorbar(fig[:, end+1], hm)
    fig
end