# February 3 2023
# test broadband dispersion

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using LMPTools
using GeographicLib
using LsqFit
#using Suppressor
using CairoMakie
using Printf

c = 2.99792458e8    # ms⁻¹
v_g = 0.9905*c      # speed of light in the EIWG

# magnetic field: magnitude, dip (from vertical), azimuth (from north)
bfield = BField.(50e-6, [π/3, π/2], 0)

# define waveguide
# 1. two-segment model
h1 = 74     # km
β1 = 0.3   # km⁻¹
h2 = 87     # km
β2 = 0.5    # km⁻¹

# "typical" earth ground 
ground = Ground(10,1e-4)
# ground = Ground(81, 4.0)
distances = [0.0, 2500e3]
species = [ Species(QE, ME, z->waitprofile(z, h2, β2), electroncollisionfrequency), 
            Species(QE, ME, z->waitprofile(z, h2, β2), electroncollisionfrequency)]

waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield[i], species[i], ground, 
            distances[i]) for i in 1:2]);

proprange = 5000e3;
rx = GroundSampler(0:10e3:proprange, Fields.Ez);

# 2. Ferguson 1980/LWPC-type waveguide
dt = DateTime(2022, 11, 07, 15,00,00)
txlat, txlon = [10, 10];
rxlat, rxlon = [47.6543, -122.3083]; # seattle

# vary frequency, propagate
# TODO March 3-5:
# 1. check possible mode precomputation strategies
# 2. resolve GroundSampler distance resolution influence on final propagate results
#    → see source code for propagate, modefinder.jl, ...
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

function buildwaveguide(dt, txlat, txlon, rxlat, rxlon)
    tx = Transmitter("", txlat, txlon, 20e3) # precompute waveguide segments
    rx = Receiver("", rxlat, rxlon, 0.0, VerticalDipole())
    geoaz = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi
    line = GeodesicLine(tx, rx)
    hcoeff = (1.345, 0.668, -0.177, -0.248, 0.096)
    bcoeff = (0.01, 0.005, -0.002, 0.01, 0.015)
    grounds, dists = groundsegments(tx, rx; resolution=20e3)
    wvgs = Vector{HomogeneousWaveguide{Species}}(undef, length(dists))
    for i in eachindex(dists)
        dist = dists[i]
        wpt = forward(line, dist)

        bfield = igrf(geoaz, wpt.lat, wpt.lon, year(dt))

        sza = zenithangle(wpt.lat, wpt.lon, dt)
        h, b = ferguson(wpt.lat, sza, dt)
        h += fourierperturbation(sza, hcoeff)
        b += fourierperturbation(sza, bcoeff)

        species = Species(QE, ME, z->waitprofile(z, h, b), electroncollisionfrequency)

        wvgs[i] = HomogeneousWaveguide(bfield, species, grounds[i], dist)
    end
    wvg = SegmentedWaveguide(wvgs)

    return wvg
end

# fit phase dispersion to final phases
function phasefit(freqs, phases; p0 = [1E-6, 0.1, 0.5])

    @. model(x, p) = p[1]*x + p[2] + p[3]*(1/x)
    # bounds
    lb = [-Inf, -Inf, 0];
    ub = [Inf, Inf, Inf];
    fit = curve_fit(model, freqs, phases, p0, lower=lb, upper=ub)
    fit

end

function iterfit(xdata, ydata, thres)

    @. model(x, p) = p[1]*x + p[2] + p[3]*(1/x)
    # bounds
    lb = [-Inf, -Inf, 0];
    ub = [Inf, Inf, Inf];
    p0 = [1E-6, 0.1, thres];
    # set up loop conditions
    xin = copy(xdata); # can remove these if preserving original data is not important
    yin = copy(ydata);
    xout = Vector{Float64}(undef, 0);
    yout = Vector{Float64}(undef, 0);
    fit = curve_fit(model, xin, yin, p0, lower=lb, upper=ub)
    sigma = stderror(fit)[3]
    while sigma > thres
        out = findmax(abs.(fit.resid));
        push!(xout, xin[out[2]]);
        push!(yout, yin[out[2]]);
        popat!(xin, out[2]);
        popat!(yin, out[2]);
        fit = curve_fit(model, xin, yin, p0, lower=lb, upper=ub)
        sigma = stderror(fit)[3] 
    end
    fit, xin, yin, xout, yout
end

# build waveguide
waveguide = buildwaveguide(dt, txlat, txlon, rxlat, rxlon);
nsegments = length(waveguide.v)
gsdist = inverse(txlon, txlat, rxlon, rxlat).dist;
proprange = gsdist;
rx = GroundSampler(0:10e3:gsdist, Fields.Ez)

# run broadband propagation
freqs = 6e3:1e3:18e3;
ωfreqs = 2*pi*freqs;
@time amps, phases = varyfreq(waveguide, rx, freqs);
# @time E, a, p = propagate(waveguide, Transmitter(freqs[1]), GroundSampler(gsdist, Fields.Ez));
# amps_hires = amps;
# phases_hires = phases;
# _hires: rx sampling interval = 10e3 m = 10 km
# @time results: 2389.612713 seconds (8.46 G allocations: 184.949 GiB, 2.32% gc time)
# _lowres: rx sampling interval = 1000e3m = 1000 km
# @time results: 2055.894965 seconds (8.47 G allocations: 185.048 GiB, 2.37% gc time)
# → sampling interval does not significantly impact runtime

# fit curve to final phase
finalphase = [phases[i][end] for i in eachindex(freqs)];
# finalphasefit = phasefit(ωfreqs, finalphase)
finalphasefit, ωin, ϕin, ωout, ϕout = iterfit([ωfreqs;], finalphase, proprange/50);
finalphasecurve = finalphasefit.param[1].*ωfreqs .+ finalphasefit.param[2] .+ finalphasefit.param[3]./(ωfreqs);
ωₒ = (finalphasefit.param[3]*2*c/proprange)^(1/2)
c3 = finalphasefit.param[3]/proprange
f₀ = (c3*4*c)^(1/2)/(2*pi) # note! factor of 4 is a fudge; parallel-plate dispersion relation indicates it should be a 2
hᵢ = c/(2f₀)

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
    global waveform = waveform + component;

    component_syn = Aₒ[j]*cos.(ωfreqs[j].*(t.-(r/c)*(1 - (ω₀_syn^2)/(ωfreqs[j]^2))^(1/2)))
    global waveform_syn = waveform_syn + component_syn;

end

# normalize waveforms
# waveform = waveform./length(ωfreqs);
# waveform_syn = waveform_syn./length(ωfreqs);
waveform = waveform./maximum(abs.(waveform));
waveform_syn = waveform_syn./maximum(abs.(waveform_syn));

# "dispersion parameter"
c3 = finalphasefit.param[3]/proprange
f₀ = (c3*2*c)^(1/2)/(2*pi)
f₀ = (finalphasefit.param[3]*2*c/proprange)^(1/2)/(2*pi)
# waveguide effective height
hᵢ = c/(2f₀)
# f₀_dn

# plot propagation path with waveguide segments
latrange = -89.5:89.5
lonrange = -179.5:179.5
latmesh = ones(360)' .* latrange
lonmesh = lonrange' .* ones(180)
sigmamap = [get_sigma(lat, lon) for lat in latrange, lon in lonrange]
wpts = waypoints(GeodesicLine(txlon, txlat, lon2=rxlon, lat2=rxlat), dist=100e3);
wlats = zeros(length(wpts));
wlons = zeros(length(wpts));
for i in eachindex(wpts)
    wlats[i] = wpts[i].lat;
    wlons[i] = wpts[i].lon;
end


begin fig = Figure(resolution = (1000,1000))
    
    # global plot properties
    fontsize_theme = Theme(fontsize=20)
    set_theme!(fontsize_theme)

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

    # ga1 = GeoAxis(fig[4,1]; coastlines = true, title = "σ",
    #     dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))

    ylims_fa1 = [-30, 90];
    # lines!(fa1, [distances[2]/1000, distances[2]/1000], ylims_fa1; color="gray")
    # ylims!(fa1, ylims_fa1);
    xlims!(fa1, [0, proprange/1e3]);

    ylims_fa2 = [-400, 400];
    # lines!(fa2, [distances[2]/1000, distances[2]/1000], ylims_fa2; color="gray")
    # ylims!(fa2, ylims_fa2);
    xlims!(fa2, [0, proprange/1e3]);

    for i in eachindex(freqs)
        lines!(fa1, rx.distance/1000, amps[i];
            linewidth=2, color=freqcolors[i+1], 
            label=string(trunc(Int, freqs[i]/1000)))
        lines!(fa2, rx.distance/1000, rad2deg.(phases[i]);
            linewidth=2, color=freqcolors[i+1])
    end

    scatter!(fa3, ωin/(2*pi*1000), rad2deg.(ϕin);
        linewidth=2, color="black",
        label="measured phase")
    # scatter!(fa3, ωout/(2*pi*1000), rad2deg.(ϕout);
    #     linewidth=2, color="red",
    #     label="outliers")
    lines!(fa3, freqs/1000, rad2deg.(finalphasecurve);
        linewidth=2, linestyle="-", color="black",
        label = "phase fit")

    lines!(fa4, tᵣ.*1e3, waveform;
        linewidth=2, color="black",
        label="simulated (LMP)")
    # lines!(fa4, tᵣ.*1e3, waveform_syn;
    #     linewidth=2, color="red",
    #     label="synthetic (Dowden+ 2002)")

    xlims!(fa4, [-0.2 1])
    # ylims!(fa4, [-100 100])

    # sf1 = surface!(ga1, lonmesh, latmesh, log10.(sigmamap); # change to log colorscale; colormap to categorical
    #     colormap=cgrad(:darkterrain, rev = true), shading=false);
    # cb1 = Colorbar(fig[4,2], sf1; label = "log10(σ / S m⁻¹)", height = Relative(0.65)) # check units
    # lines!(ga1, wlons, wlats; color=:yellow, lineweight=10)
    # scatter!(ga1, txlon, txlat; color=:green, markersize=10)
    # scatter!(ga1, rxlon, rxlat; color=:red, markersize=10)

    # legends
    axislegend(fa3)
    axislegend(fa4)
    fig[1:2, 2] = Legend(fig, fa1, "frequency (kHz)", framevisible=false)

    # supertitle = Label(fig[0, :], "Broadband sferic propagation\n segment 1: d = 2500 km, h' = 75 km, β = 0.35 km⁻¹\n segment 2: d = 2500 km, h' = 82 km, β = 0.50 km⁻¹"; fontsize=20)
    superstr = @sprintf("Broadband sferic propagation\nnumber of segments = %g\nf₀ = %.2f kHz; hᵢ = %.2f km", nsegments, f₀/1e3, hᵢ/1e3)
    supertitle = Label(fig[0, :], superstr; fontsize=20)
    fig
    # save("figures/sample_sferic_dispersion_1000_LMPonly.png", fig, px_per_unit=1)
end
