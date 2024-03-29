# BuildWaveguide.jl
# Todd Anderson
# 17 September 2022
#
# Construct SegmentedWaveguide for any great circle propagation path based
# on land/sea/ice mask and day/night determination

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using MAT
using GeoMakie
using CairoMakie
import GMT.geodesic as geodesic
import Distances.haversine as haversine


# in case of InitError: try Pkg.build(jll_module), where jll_module is the 
# module called out at the end of the InitError.
# InitError tests:
#   import GMT -> subsequent error on using GeoMakie, using CairoMakie
#   using GMT  -> same error
#   import GMT.geodesic as geodesic -> same error
# Looks like GMT needs to imported after GeoMakie and CairoMakie.


const R_KM = 6371.8

begin # import .mat file data
    lsifilename = "LSI_mask.mat"
    lsifile = matopen(lsifilename)
    maskstruct = read(lsifile, "mask")
    close(lsifile)
end

latmesh = maskstruct["lat_mesh"];
lonmesh = maskstruct["lon_mesh"];
LSI = maskstruct["LSI"];

# # GeoMakie examples
# let fig = Figure()
#     ga = GeoAxis(
#         fig[1, 1]; # any cell of the figure's layout
#         dest = "+proj=wintri", # the CRS in which you want to plot
#         coastlines = true # plot coastlines from Natural Earth, as a reference.
#     );
#     scatter!(ga, -120:15:120, -60:7.5:60; color = -60:7.5:60, strokecolor = (:black, 0.2));
#     fig
# end

# begin # plot variables
#     fieldlons = -180:180; 
#     fieldlats = -90:90;
#     field = [exp(cosd(lon)) + 3(lat/90) for lon in fieldlons, lat in fieldlats];

#     img = rotr90(GeoMakie.earth());
#     land = GeoMakie.land();
# end

# begin fig = Figure(resolution = (1000, 1000));

#     ga1 = GeoAxis(fig[1, 1]; dest = "+proj=ortho", coastlines = true, lonlims = (-90, 90), title = "Orthographic\n ")
#     ga2 = GeoAxis(fig[1, 2]; dest = "+proj=moll", title = "Image of Earth\n ")
#     ga3 = GeoAxis(fig[2, 1]; coastlines = false, title = "Plotting polygons")
#     ga4 = GeoAxis(fig[2, 2]; dest = "+proj=natearth", title = "Auto limits") # you can plot geodata on regular axes too

#     surface!(ga1, fieldlons, fieldlats, field; colormap = :rainbow_bgyrm_35_85_c69_n256, shading = false);
#     image!(ga2, -180..180, -90..90, img; interpolate = false); # this must be included
#     poly!(ga3, land[50:100]; color = 1:51, colormap = (:plasma, 0.5));
#     poly!(ga4, land[22]); datalims!(ga4);

#     fig
# end

Tx = [50.0 10.0]
Rx = [47.6062 -122.3321]
sspath = geodesic([reverse(Tx); reverse(Rx)], step=10, unit=:k)

gnd = function getground(loc, latmesh, lonmesh, mask)
    #loc = sspath[1:2,:]
    lon = loc[:,1]
    lat = loc[:,2]
    gnd = Vector{Float64}(undef, length(lon))
    for i in eachindex(lon)
        meshlat = findmin(abs.(lat[i].-latmesh[:,1]))[2]
        meshlon = findmin(abs.(lon[i].-lonmesh[1,:]))[2]
        gnd[i] = mask[meshlat, meshlon]
    end
    return gnd
end

gnd1 = getground(sspath[1,:]', latmesh, lonmesh, LSI) # not sure why sspath[1,:] returns a column vector
pathgnd = getground(sspath, latmesh, lonmesh, LSI)

# plot LSI mask with Tx-Rx path, colored by LSI value
# let fig = Figure(resolution = (1200,800))
#     ga = GeoAxis(fig[1,1]; coastlines = true, title = "Land-sea-ice mask")

#     surface!(ga, lonmesh, latmesh, LSI; 
#         colormap="broc", colorrange=(-5,5), shading=false);
#     scatter!(ga, reverse(Tx); color=:green)
#     scatter!(ga, reverse(Rx); color=:red)
#     lines!(ga, sspath; color=pathgnd, colormap="berlin", colorrange=(-2,2))

#     fig
# end

# find contiguous segments on geodesic that share ground values
# define function that takes sspath and pathgnd as input and outputs Mx2 matrix
# representing the length and ground value of each of M segments.  check that 
# the sum of the segment lengths is equal to the geodesic length
function segments(path, ground; segmentlimit="none")
    grounddiff = findall(abs.(diff(ground)).>0)
    nsegments = length(grounddiff)+1;
    segment_length = Vector{Float64}(undef, nsegments)
    segment_ground_flag = Vector{Float64}(undef, nsegments)
    segment_ground = Vector{Ground}(undef, nsegments)
    
    # first segment
    segment_ground_flag[1] = ground[1]
    segment_length[1] = haversine(path[1,:], path[grounddiff[1],:], R_KM)
    # subsequent segments
    for i in 2:length(segment_length)-1
        segment_ground_flag[i] = ground[grounddiff[i]]
        segment_start = path[grounddiff[i-1],:]
        segment_end = path[grounddiff[i],:]
        segment_length[i] = haversine(segment_start, segment_end, R_KM)
    end
    # last segment
    segment_ground_flag[end] = ground[end]
    segment_length[end] = haversine(path[grounddiff[end],:], path[end,:], R_KM)

    # match segment mask values with GROUND states
    land = findall(x->x==1, segment_ground_flag)
    sea  = findall(x->x==-1, segment_ground_flag)
    ice  = findall(x->x==0, segment_ground_flag)

    # combine segments, starting with the shortest segment, until 
    # number of segments <= segmentlimit and all segments are longer 
    # than segmentminlength
    if typeof(segmentlimit) <: Number
        while nsegments > segmentlimit
            # reduce number of segments:
            #   find smallest segment
            #   set segment mask value equal to value of adjacent segment
            #   combine all adjacent segments that share mask values
        end
    end

    # TODO: add option to reduce segments to one of each different component type, 
    # with lengths equal to the total length of that waveguide component.
    # E.g. reduce to one land, one sea and one ice segment, where the length of the
    # land segment is equal to the sum of the lengths of the individual land segments
    if segmentlimit == "reduce"
        # TODO: order these in order of appearrance on propagation path
        segment_total_length[1] = sum(segment_length[land])
        segment_total_length[2] = sum(segment_length[sea])
        segment_total_length[3] = sum(segment_length[ice])

        segment_total_ground[1] = Ref(GROUND[5])
        segment_total_ground[2] = Ref(GROUND[10])
        segment_total_ground[3] = Ref(GROUND[1])

        segment_total_ground_flag[1] = 1;
        segment_total_ground_flag[1] = -1;
        segment_total_ground_flag[1] = 0;
    end


    

    segment_ground[land].= Ref(GROUND[5])
    segment_ground[sea] .= Ref(GROUND[10])
    segment_ground[ice] .= Ref(GROUND[1])
    

    return segment_length, segment_ground_flag, segment_ground
end

segment_length, segment_ground_flag, segment_ground = segments(sspath, pathgnd)

# check sum of segment lengths equals distance between Tx and Rx
sum(segment_length)
total_distance = haversine(sspath[1,:], sspath[end,:], R_KM)

## build SegmentedWaveguide
tx = Transmitter(24e3)
rx = GroundSampler(0:10e3:10000e3, Fields.Ez)
rx_station = GroundSampler(total_distance, Fields.Ez)

h1 = 75     # km
β1 = 0.35   # km⁻¹

bfield = BField(50e-6, π/2, 0)

distances = vcat(0.0, segment_length)
pop!(distances); 
ground = segment_ground
species = Species(QE, ME, z->waitprofile(z, h1, β1), electroncollisionfrequency)

waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield, species, ground[i], 
            distances[i]) for i in eachindex(distances)]);

#@time E, a, p = propagate(waveguide, tx, rx);
# timing results 
#   with 13 segments, tx = Transmitter(24e3); rx = GroundSampler(0:10e3:10000e3, Fields.Ez):
# 469.944245 seconds (938.84 M allocations: 20.459 GiB, 1.69% gc time, 2.04% compilation time)
# -> need to be much faster!
# TODO: combine segments, starting with the shortest segment, until there are no more than 
#   (5?) segments.
#   BUT FIRST: test effect of increasing number of segments vs. increasing total distance
# TODO: what we need is a prediction of amplitude and phase at the receiver location, not
#   at a range of distances including the receiver location.  Try timing the same 
#   SegmentedWaveguide run, but with rx = GroundSampler(distance_to_station, Fields.Ez)

#@time E_s, a_s, p_s = propagate(waveguide, tx, rx_station)
# timing results
#   with 13 segments, tx = Transmitter(24e3); rx = GroundSampler(total_distance, Fields.Ez):
# 270.606685 seconds (931.78 M allocations: 20.124 GiB, 2.25% gc time, 0.31% compilation time)
# -> not 1001 times faster! try fewer segments

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

freqs = 6e3:1e3:24e3

#@time amps, phases = varyfreq(waveguide, rx_station, freqs)
# timing results (13 segments, one sample location, 11 frequencies):
#   924.117816 seconds (2.98 G allocations: 74.082 GiB, 3.34% gc time, 0.01% compilation time)

@time amps, phases = varyfreq(waveguide, rx, freqs)
# timing results (13 segments, 1001 sample locations, 11 frequencies):
#   1209.574735 seconds (2.98 G allocations: 74.082 GiB, 1.52% gc time, 0.02% compilation time)
#   interesting! number of sample locations does not impact allocations at all!
# timing results (13 segments, 1001 sample locations, 19 frequencies (6-24 kHz)):
#   1948.517982 seconds (8.94 G allocations: 206.378 GiB, 2.05% gc time)
#   odd that allocations increased so much; look into how higher frequencies are handled
#   see also phase behavior differences between f<=16 kHz and f>=17 kHz


# plot propagation path with waveguide segments
begin fig = Figure(resolution = (1200,1200))
    
    ga = GeoAxis(fig[1:2,1:2]; coastlines = true, title = "Land-sea-ice mask",
        dest = "+proj=natearth", latlims = (40,90), lonlims = (-140, 30))

    surface!(ga, lonmesh, latmesh, LSI; 
        colormap="broc", colorrange=(-5,5), shading=false);
    lines!(ga, sspath; color=pathgnd, colormap="lisbon", colorrange=(-2,2),
        linewidth=10)
    scatter!(ga, reverse(Tx); color=:green, markersize=10)
    scatter!(ga, reverse(Rx); color=:red, markersize=10)
    text!(ga, "Tx", position = reverse(Tx), align=(:left, :top))
    text!(ga, "Rx", position = reverse(Rx), align=(:left, :top))

    fa1 = Axis(fig[3,1], title="waveguide segments",
        xlabel="segment length (km)",
        ylabel="segment number",
        yreversed=true)
    barplot!(fa1, 1:length(segment_length), segment_length; 
        direction=:x,
        color=segment_ground_flag, colormap="lisbon", colorrange=(-2,2))
    
    # ground legend!
    colors = cgrad(:lisbon, 5, categorical=true, rev=true)
    labels = ["land: 15, 1e-3", "ice:     5, 1e-5", "sea:  81, 4"]
    elements = [PolyElement(polycolor=colors[i+1]) for i in eachindex(labels)]
    title = "Ground: ϵᵣ, σ"
    Legend(fig[3,2],elements, labels, title, framevisible=false)

    # amplitude and phase for each frequency
    freqcolors = cgrad(:thermal, length(freqs)+4, categorical=true, rev=false)

    fa2 = Axis(fig[4,1], title="amplitude",
        xlabel="distance (km)",
        ylabel="amplitude (dB)")

    fa3 = Axis(fig[5,1], title="phase",
        xlabel="distance (km)",
        ylabel="phase (∘)")

    for i in eachindex(freqs)
        lines!(fa2, rx.distance/1000, amps[i];
            linewidth=2, color=freqcolors[i+1], 
            label=string(trunc(Int, freqs[i]/1000)))
        lines!(fa3, rx.distance/1000, rad2deg.(phases[i]);
            linewidth=2, color=freqcolors[i+1])
    end
    
    # frequency legend!
    fig[4:5, 2] = Legend(fig, fa2, "frequency (kHz)", framevisible=false)

    fig
    # save("LSIpath_segments_amp_phase_freq_6-24.png", fig, px_per_unit=1)
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
using CairoMakie

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


begin fig = Figure(resolution = (1000,600))

    # global plot properties
    fontsize_theme = Theme(fontsize=20)
    set_theme!(fontsize_theme)
    
    ga1 = GeoAxis(fig[1,1]; coastlines = true, title = "σ",
        dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))

    ga2 = GeoAxis(fig[2,1]; coastlines = true, title = "ϵᵣ",
        dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))

    ga3 = GeoAxis(fig[1,3]; coastlines = true, title = "h'",
        dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))
    
    ga4 = GeoAxis(fig[2,3]; coastlines = true, title = "β",
        dest = "+proj=natearth", latlims = (-90,90), lonlims = (-180, 180))

    sf1 = surface!(ga1, lonmesh, latmesh, log10.(sigmamap); # change to log colorscale; colormap to categorical
        colormap=cgrad(:darkterrain, rev = true), shading=false);
    cb1 = Colorbar(fig[1,2], sf1; label = "log10(σ / S m⁻¹)", height = Relative(0.65)) # check units

    # TODO: change colorbar scaling
    epscolors = cgrad(:vik, rev=true, 100, categorical=true)
    epscolormap = epscolors[sort(unique(epsilonmap))]
    sf2 = contourf!(ga2, lonrange, latrange, epsilonmap';
        colormap=epscolormap, levels = vcat(0,sort(unique(epsilonmap.+1))), shading=false);
    cb2 = Colorbar(fig[2,2], sf2; label = "F m⁻¹", height = Relative(0.65), ticks = sort(unique(epsilonmap))) # check units

    sf3 = surface!(ga3, lonmesh, latmesh, hmap; 
        colormap=Reverse(:tokyo), shading=false);
    cb3 = Colorbar(fig[1,4], sf3; label = "km", height = Relative(0.65))

    sf4 = surface!(ga4, lonmesh, latmesh, βmap; 
        colormap=Reverse(:oslo), shading=false);
    cb4 = Colorbar(fig[2,4], sf4; label = "km⁻¹", height = Relative(0.65))

    superstr = @sprintf("Global waveguide parameters\n%s",dt)
    supertitle = Label(fig[0, :], superstr; fontsize=20)

    fig
    # save("figures/sample_waveguide_1000.png", fig, px_per_unit=1)
end

