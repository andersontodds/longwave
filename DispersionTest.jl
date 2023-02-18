# February 3 2023
# test broadband dispersion

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using LsqFit
#using Suppressor
using CairoMakie
using Printf

c = 2.99792458e8    # ms⁻¹
v_g = 0.9905*c      # speed of light in the EIWG

# simple SegmentedWaveguide example with 2 segments
h1 = 75     # km
β1 = 0.35   # km⁻¹
h2 = 82     # km
β2 = 0.5    # km⁻¹

# "typical" earth ground 
ground = Ground(10,1e-4)
# ground = Ground(81, 4.0)

# vertical magnetic field
bfield = BField(50e-6, π/2, 0)

# define waveguide
distances = [0.0, 2500e3]
species = [ Species(QE, ME, z->waitprofile(z, h1, β1), electroncollisionfrequency), 
            Species(QE, ME, z->waitprofile(z, h2, β2), electroncollisionfrequency)]

waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield, species[i], ground, 
            distances[i]) for i in 1:2]);

proprange = 5000e3;
rx = GroundSampler(0:10e3:proprange, Fields.Ez);

# vary frequency, propagate and sample only at station
function varyfreq(waveguide, rx, freqs)
    amps = Vector{Vector{Float64}}(undef, length(freqs))
    phases = Vector{Vector{Float64}}(undef, length(freqs))
    for i in eachindex(freqs)
        tx = Transmitter(freqs[i])
        E, a, p = propagate(waveguide, tx, rx);
        amps[i] = a
        phases[i] = p
    end
    return amps, phases
end

# fit phase dispersion to final phases
function phasefit(freqs, phases; p0 = [1E-6, 0, 0.5])

    @. model(x, p) = p[1]*x + p[2] + p[3]*(1/x)
    #p0 = [1E-6, 0, 0.5]
    fit = curve_fit(model, freqs, phases, p0)
    fit

end

# run broadband propagation
freqs = 6e3:1e3:18e3;
ωfreqs = 2*pi*freqs;
@time amps, phases = varyfreq(waveguide, rx, freqs);

# fit curve to final phase
finalphase = [phases[i][end] for i in eachindex(freqs)];
finalphasefit = phasefit(ωfreqs, finalphase)
finalphasecurve = finalphasefit.param[1].*ωfreqs .+ finalphasefit.param[2] .+ finalphasefit.param[3]./(ωfreqs);
ωₒ = (finalphasefit.param[3]*2*c/proprange)^(1/2)

# generate simulated and synthetic sferics
# simulated: use LMP-propagated amplitudes and phases, with Dowden inital amplitudes
# synthetic: use Dowden amplitudes and phases
r = proprange;
tᵣ = -0.2e-3:1e-6:1e-3;
t = tᵣ .+ r/c;
# initial amplitude parameters; see synthetic sferics in scratch.jl and Dowden 2002
ωₐ = 2*pi*14e3;     # frequency of peak spectral density: ~12 kHz
ωᵣ = 2*pi*11e3;     # tune this
ω₀_syn = 2*pi*1.6e3;    # waveguide cutoff frequency

# get final amplitudes
finalamps = [amps[i][end] for i in eachindex(freqs)];
Aₒ = cos.(pi*(ωfreqs.-ωₐ)./(2*ωᵣ)).^2;
A = Aₒ.*10.0.^(finalamps./20);

# build synthetic and simulated waveforms
waveform = zeros(length(t));
waveform_syn = zeros(length(t));
for j in eachindex(freqs)
    # component = amps[j][end]*cos.(ωfreqs[j]*(t.-proprange/v_g));
    
    component = A[j]*cos.(ωfreqs[j].*(t.-(r/c)*(1 - ωₒ^2/ωfreqs[j]^2)^(1/2)))
    waveform = waveform + component;

    component_syn = Aₒ[j]*cos.(ωfreqs[j].*(t.-(r/c)*(1 - (ω₀_syn^2)/(ωfreqs[j]^2))^(1/2)))
    waveform_syn = waveform_syn + component_syn;

end

# normalize waveforms
waveform = waveform./length(ωfreqs);
waveform_syn = waveform_syn./length(ωfreqs);

# plot propagation path with waveguide segments
begin fig = Figure(resolution = (1200,1200))

    # amplitude and phase for each frequency
    freqcolors = cgrad(:thermal, length(freqs)+4, categorical=true, rev=false)

    fa1 = Axis(fig[1,1], title="amplitude relative to source",
        xlabel="distance (km)",
        ylabel="amplitude (dB)")

    fa2 = Axis(fig[2,1], title="phase relative to source",
        xlabel="distance (km)",
        ylabel="phase (∘)")

    fa3 = Axis(fig[1:2, 3], title=@sprintf("phase fit at r = %d km", proprange/1e3),
        xlabel="frequency (kHz)",
        ylabel="phase (∘)")    

    fa4 = Axis(fig[3,1:3], title=@sprintf("waveform at r = %d km", proprange/1e3),
        xlabel="t - r/c (ms)",
        ylabel="amplitude (normalized)")

    ylims_fa1 = [-30, 90];
    lines!(fa1, [distances[2]/1000, distances[2]/1000], ylims_fa1; color="gray")
    ylims!(fa1, ylims_fa1);
    xlims!(fa1, [0, 5000]);

    ylims_fa2 = [-400, 400];
    lines!(fa2, [distances[2]/1000, distances[2]/1000], ylims_fa2; color="gray")
    ylims!(fa2, ylims_fa2);
    xlims!(fa2, [0, 5000]);

    for i in eachindex(freqs)
        lines!(fa1, rx.distance/1000, amps[i];
            linewidth=2, color=freqcolors[i+1], 
            label=string(trunc(Int, freqs[i]/1000)))
        lines!(fa2, rx.distance/1000, rad2deg.(phases[i]);
            linewidth=2, color=freqcolors[i+1])
    end

    scatter!(fa3, freqs/1000, rad2deg.(phases[i][end] for i in eachindex(freqs));
        linewidth=2, color="black",
        label="measured phase")
    lines!(fa3, freqs/1000, rad2deg.(finalphasecurve);
        linewidth=2, linestyle="-", color="red",
        label = "phase fit")

    lines!(fa4, tᵣ.*1e3, waveform;
        linewidth=2, color="black",
        label="simulated (LMP)")
    lines!(fa4, tᵣ.*1e3, waveform_syn;
        linewidth=2, color="red",
        label="synthetic (Dowden+ 2002)")

    xlims!(fa4, [-0.2 1])
    # ylims!(fa4, [-100 100])

    # legends
    axislegend(fa3)
    axislegend(fa4)
    fig[1:2, 2] = Legend(fig, fa1, "frequency (kHz)", framevisible=false)

    supertitle = Label(fig[0, :], "Broadband sferic propagation\n segment 1: d = 2500 km, h = 75 km, β = 0.35 km⁻¹\n segment 2: d = 2500 km, h = 82 km, β = 0.50 km⁻¹"; fontsize=20)
    fig
    # save("sample_sferic_dispersion.png", fig, px_per_unit=1)
end

# "dispersion parameter"
c3 = finalphasefit.param[3]/proprange
f₀ = (finalphasefit.param[3]*2*c/proprange)^(1/2)/(2*pi)