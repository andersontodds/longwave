# February 3 2023
# test broadband dispersion

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using LsqFit
using CairoMakie

# simple SegmentedWaveguide example with 2 segments
h1 = 75     # km
β1 = 0.35   # km⁻¹
h2 = 82     # km
β2 = 0.5    # km⁻¹

# "typical" earth ground 
ground = Ground(10,1e-4)

# vertical magnetic field
bfield = BField(50e-6, π/2, 0)

# define waveguide
distances = [0.0, 2500e3]
species = [ Species(QE, ME, z->waitprofile(z, h1, β1), electroncollisionfrequency), 
            Species(QE, ME, z->waitprofile(z, h2, β2), electroncollisionfrequency)]

waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield, species[i], ground, 
            distances[i]) for i in 1:2]);

rx = GroundSampler(0:10e3:5000e3, Fields.Ez);

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
function phasefit(freqs, phases)

    @. model(x, p) = p[1]*x + p[2] + p[3]*(1/x)
    p0 = [1E-6, 0, 0.5]
    fit = curve_fit(model, freqs, phases, p0)
    fit

end

# run broadband propagation
freqs = 6e3:1e3:18e3;
@time amps, phases = varyfreq(waveguide, rx, freqs);

finalphase = [phases[i][end] for i in eachindex(freqs)];

finalphasefit = phasefit(freqs, finalphase)
finalphasecurve = finalphasefit.param[1].*freqs .+ finalphasefit.param[2] .+ finalphasefit.param[3]./freqs;


# plot propagation path with waveguide segments
begin fig = Figure(resolution = (1200,1200))

    # amplitude and phase for each frequency
    freqcolors = cgrad(:thermal, length(freqs)+4, categorical=true, rev=false)

    fa1 = Axis(fig[1,1], title="amplitude",
        xlabel="distance (km)",
        ylabel="amplitude (dB)")

    fa2 = Axis(fig[2,1], title="phase",
        xlabel="distance (km)",
        ylabel="phase (∘)")

    fa3 = Axis(fig[1:2, 3], title="dispersion",
        xlabel="frequency (kHz)",
        ylabel="phase (∘)")    

    for i in eachindex(freqs)
        lines!(fa1, rx.distance/1000, amps[i];
            linewidth=2, color=freqcolors[i+1], 
            label=string(trunc(Int, freqs[i]/1000)))
        lines!(fa2, rx.distance/1000, rad2deg.(phases[i]);
            linewidth=2, color=freqcolors[i+1])
    end
    
    lines!(fa3, freqs/1000, rad2deg.(phases[i][end] for i in eachindex(freqs));
        linewidth=2, color="black",
        label="measured phase")
    lines!(fa3, freqs/1000, rad2deg.(finalphasecurve);
        linewidth=2, linestyle="-", color="red",
        label = "phase fit")

    # legends
    axislegend(fa3)
    fig[1:2, 2] = Legend(fig, fa1, "frequency (kHz)", framevisible=false)

    fig
    # save("LSIpath_segments_amp_phase_freq_6-24.png", fig, px_per_unit=1)
end