# iteratively exclude data from fit until fit error is below threshold
# 2/27 iterfit issues:
# - initial fit curve can be nearer to phase-jump points than some non-jump points;
#   subsequent iterations will reject non-jump points
#   ↪ fix this with narrower inital conditions?  
using LsqFit
using SolarAngles

# data from longwave/DispersionTest.jl
# 2π corrections added
ωf = [2*pi*(6e3:1e3:18e3);];
ϕ = [   5.840479829632033   #+ 2*π
        3.7609869647037457  #+ 2*π
        8.513717721555789   
        7.0863922989541015  
       -0.35831763785328985 #+ 2*π
        5.149585430769938
        4.344839597842935
        9.638081754389809   #- 2*π
        8.942412906068212   #- 2*π
        1.9855723961852674
        1.398661708163204
        1.1841425901930442
        0.7130530291814566];
r = 5000e3;
thres = r/50; # maybe put a bit more thought into threshold choice

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

# test n*2π shifting
# rules for ϕ(ω) fit:
#   1. minimize slope of ϕ(ω)
#   2. ϕ(ω) must be concave up
#   3. no ω components may be further than <threshold> from fit
# if Δϕ/Δω > -2π s⁻¹ (i.e. -2π radians for ω steps of 2π rad*1kHz), can shift data  
# until this criterion is met 
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

# xyfit, xin, yin, xout, yout = iterfit(ωf, ϕ, thres);
ϕfit, ωshift, ϕshift = shiftfit(ωf, ϕ, thres);
fitcurve = ϕfit.param[1].*ωf .+ ϕfit.param[2] .+ ϕfit.param[3]./(ωf);
fit_f₀ = (ϕfit.param[3]*2*c/r)^(1/2)/(2*pi);

begin fig = Figure()
    fa = Axis(fig[1, 1], title=@sprintf("phase fit at r = %d km\nf₀ = %.2f kHz", r/1e3, fit_f₀/1e3),
        xlabel="frequency (kHz)",
        ylabel="phase (∘)")

    scatter!(fa, ωf./(2*pi*1000), rad2deg.(ϕ);
        linewidth=2, color="black",
        label="unshifted")

    scatter!(fa, ωshift./(2*pi*1000), rad2deg.(ϕshift);
        linewidth=2, color="red",
        label="shifted")

    # scatter!(fa, xout./(2*pi*1000), rad2deg.(yout);
    #     linewidth=2, color="red",
    #     label="out")

    lines!(fa, ωf/(2*pi*1000), rad2deg.(fitcurve);
        linewidth=2, linestyle="-", color="black",
        label = "phase fit")

    axislegend()
    
    fig
end

# waveguide building: identify ionosphere h', β and ground ϵᵣ, σ based on
# LWPM-style ionosphere model and land-sea-ice mask 
# Example from @fgasdia's LMPTools.jl documentation
using LongwaveModePropagator
using LongwaveModePropagator: ME, QE
using LMPTools
using GeographicLib
using Dates, Printf
using GeoMakie

dt = DateTime(2022, 11, 07, 15,00,00)

# Propagation path from NAA in Maine to Boulder, Colorado.
tx = TRANSMITTER[:NAA]
rx = Receiver("Boulder", 40.01, -105.244, 0.0, VerticalDipole())

# For efficiency, precompute the geographic azimuth between the transmitter and receiver.
geoaz = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi

# and precompute a line between the transmitter and receiver.
line = GeodesicLine(tx, rx)

# Fourier perturbation coefficients
hcoeff = (1.345, 0.668, -0.177, -0.248, 0.096)
bcoeff = (0.01, 0.005, -0.002, 0.01, 0.015)

# Use `groundsegments` from LMPTools to divide the propagation path into segments
grounds, dists = groundsegments(tx, rx; resolution=20e3)

# Preallocate a vector of `HomogeneousWaveguide` to construct a `SegmentedWaveguide`
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

gs = GroundSampler(range(tx, rx), Fields.Ez)
E, amplitude, phase = propagate(wvg, tx, gs)  # field at the receiver

# plot global waveguide parameters
latrange = -89.5:89.5
lonrange = -179.5:179.5
latmesh = ones(360)' .* latrange
lonmesh = lonrange' .* ones(180)

sigmamap = [get_sigma(lat, lon) for lat in latrange, lon in lonrange]
epsilonmap = [get_epsilon(lat, lon) for lat in latrange, lon in lonrange]

szamap = [zenithangle(lat, lon, dt) for lat in latrange, lon in lonrange]

ionomap = [ferguson(lat, zenithangle(lat, lon, dt), dt) for lat in latrange, lon in lonrange]
hmap = zeros(length(latrange), length(lonrange));
βmap = zeros(length(latrange), length(lonrange));
for i = 1:length(latrange), j = 1:length(lonrange)
    hmap[i,j] = ionomap[i,j][1]
    βmap[i,j] = ionomap[i,j][2]
end


begin fig = Figure(resolution = (1200,800))
    
    ga1 = GeoAxis(fig[1,1]; coastlines = true, title = "σ",
        dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))

    ga2 = GeoAxis(fig[2,1]; coastlines = true, title = "ϵᵣ",
        dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))

    ga3 = GeoAxis(fig[1,2]; coastlines = true, title = "h'",
        dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))
    
    ga4 = GeoAxis(fig[2,2]; coastlines = true, title = "β",
        dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))

    surface!(ga1, lonmesh, latmesh, log10.(sigmamap); 
        colormap=Reverse(:broc), shading=false);

    surface!(ga2, lonmesh, latmesh, epsilonmap; 
        colormap=Reverse(:vik), shading=false);

    surface!(ga3, lonmesh, latmesh, hmap; 
        colormap=Reverse(:tokyo), shading=false);

    surface!(ga4, lonmesh, latmesh, βmap; 
        colormap=Reverse(:oslo), shading=false);

    superstr = @sprintf("Global waveguide parameters\n%s",dt)
    supertitle = Label(fig[0, :], superstr; fontsize=20)

    fig
end

# Wait profile
# from Wait and Spies
h = 75:5:85 # km
β = 0.3:0.1:0.5 # km⁻¹
z = 1:1:100 # km
#Nₑ = 1.43e13*exp(-0.15*h).*exp.((β-0.15).*(z.-h))
Nₑ_LMP = waitprofile.(z.*1000, h[1], β[1])

begin fig = Figure()

    colors = cgrad(:matter, length(h)*length(β), categorical=true)

    fa = Axis(fig[1, 1], title=@sprintf("Wait ionosphere profiles"),
        xlabel="electron density (m⁻³)",
        ylabel="altitude (km)",
        xscale=log10)

    for j in eachindex(β), i in eachindex(h) 
        Nₑ = 1.43e13*exp(-0.15*h[i]).*exp.((β[j]-0.15).*(z.-h[i]))
        Nₑ_LMP = waitprofile.(z.*1000, h[i], β[j])
        # fa = Axis(fig[i, j], title=@sprintf("Wait profile with h'=%d km, β=%0.2g km⁻¹", h[i], β[j]),
        #     xlabel="electron density (m⁻³)",
        #     ylabel="phase (∘)",
        #     xscale=log10)

        # lines!(fa, Nₑ, z;
        #     linewidth=2, color = colors[i + (j-1)*length(h)],
        #     label=@sprintf("h'=%d km, β=%0.2g km⁻¹", h[i], β[j]))

        lines!(fa, Nₑ_LMP, z;
            linewidth=2, color=colors[i + (j-1)*length(h)],
            label=@sprintf("LMP; h'=%d km, β=%0.2g km⁻¹", h[i], β[j]))

        hlines!(fa, h;
            color="gray", linestyle="..")

    end

    axislegend(position = :rb)

    fig
end