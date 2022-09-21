# BuildWaveguide.jl
# Todd Anderson
# 17 September 2022
#
# Construct SegmentedWaveguide for any great circle propagation path based
# on land/sea/ice mask and day/night determination

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using MAT
using GMT, Distances
using GeoMakie, GLMakie

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

# GeoMakie examples
let fig = Figure()
    ga = GeoAxis(
        fig[1, 1]; # any cell of the figure's layout
        dest = "+proj=wintri", # the CRS in which you want to plot
        coastlines = true # plot coastlines from Natural Earth, as a reference.
    );
    scatter!(ga, -120:15:120, -60:7.5:60; color = -60:7.5:60, strokecolor = (:black, 0.2));
    fig
end

begin # plot variables
    fieldlons = -180:180; 
    fieldlats = -90:90;
    field = [exp(cosd(lon)) + 3(lat/90) for lon in fieldlons, lat in fieldlats];

    img = rotr90(GeoMakie.earth());
    land = GeoMakie.land();
end

let fig = Figure(resolution = (1000, 1000));

    ga1 = GeoAxis(fig[1, 1]; dest = "+proj=ortho", coastlines = true, lonlims = (-90, 90), title = "Orthographic\n ")
    ga2 = GeoAxis(fig[1, 2]; dest = "+proj=moll", title = "Image of Earth\n ")
    ga3 = GeoAxis(fig[2, 1]; coastlines = false, title = "Plotting polygons")
    ga4 = GeoAxis(fig[2, 2]; dest = "+proj=natearth", title = "Auto limits") # you can plot geodata on regular axes too

    surface!(ga1, fieldlons, fieldlats, field; colormap = :rainbow_bgyrm_35_85_c69_n256, shading = false);
    image!(ga2, -180..180, -90..90, img; interpolate = false); # this must be included
    poly!(ga3, land[50:100]; color = 1:51, colormap = (:plasma, 0.5));
    poly!(ga4, land[22]); datalims!(ga4);

    fig
end


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
function segments(path, ground)
    grounddiff = findall(abs.(diff(ground)).>0)
    segment_length = Vector{Float64}(undef, length(grounddiff)+1)
    segment_ground_flag = Vector{Float64}(undef, length(grounddiff)+1)
    segment_ground = Vector{Ground}(undef, length(grounddiff)+1)
    segment_ground_flag[1] = ground[1]
    segment_length[1] = haversine(path[1,:], path[grounddiff[1],:], R_KM)
    for i in 2:length(segment_length)-1
        segment_ground_flag[i] = ground[grounddiff[i]]
        segment_start = path[grounddiff[i-1],:]
        segment_end = path[grounddiff[i],:]
        segment_length[i] = haversine(segment_start, segment_end, R_KM)
    end
    segment_ground_flag[end] = ground[end]
    segment_length[end] = haversine(path[grounddiff[end],:], path[end,:], R_KM)

    land = findall(x->x==1, segment_ground_flag)
    sea  = findall(x->x==-1, segment_ground_flag)
    ice  = findall(x->x==0, segment_ground_flag)

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

@time E, a, p = propagate(waveguide, tx, rx);
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

@time E_s, a_s, p_s = propagate(waveguide, tx, rx_station)
# timing results
#   with 13 segments, tx = Transmitter(24e3); rx = GroundSampler(total_distance, Fields.Ez):
# 270.606685 seconds (931.78 M allocations: 20.124 GiB, 2.25% gc time, 0.31% compilation time)
# -> not 1001 times faster! try fewer segments

# plot propagation path with waveguide segments
let fig = Figure(resolution = (1200,1200))
    
    ga = GeoAxis(fig[1:2,1:2]; coastlines = true, title = "Land-sea-ice mask",
        dest = "+proj=natearth", latlims = (40,90), lonlims = (-140, 30))

    surface!(ga, lonmesh, latmesh, LSI; 
        colormap="broc", colorrange=(-5,5), shading=false);
    lines!(ga, sspath; color=pathgnd, colormap="berlin", colorrange=(-2,2),
        linewidth=10)
    scatter!(ga, reverse(Tx); color=:green, markersize=10)
    scatter!(ga, reverse(Rx); color=:red, markersize=10)
    GLMakie.text!(ga, "Tx", position = reverse(Tx), align=(:left, :top))
    GLMakie.text!(ga, "Rx", position = reverse(Rx), align=(:left, :top))

    fa1 = Axis(fig[3,1], title="waveguide segments",
        xlabel="segment length (km)",
        ylabel="segment number",
        yreversed=true)
    barplot!(fa1, 1:length(segment_length), segment_length; 
        direction=:x, stack = [1,1,1,1,1,2,2,2,2,3,3,3,3],
        color=segment_ground_flag, colormap="berlin", colorrange=(-2,2))
    # make a legend!!!
    colors = cgrad(:berlin, 5, categorical=true, rev=true)
    labels = ["land: 15, 1e-3", "ice:     5, 1e-5", "sea:  81, 4"]
    elements = [PolyElement(polycolor=colors[i+1]) for i in eachindex(labels)]
    title = "Ground: ϵᵣ, σ"
    Legend(fig[3,2],elements, labels, title)

    fa2 = Axis(fig[4,1:2], title="amplitude: f = 24 kHz",
        xlabel="distance (km)",
        ylabel="amplitude (dB)")
    lines!(rx.distance/1000, a;
        linewidth=1.5)

    fa3 = Axis(fig[5,1:2], title="phase: f = 24 kHz",
        xlabel="distance (km)",
        ylabel="phase (∘)")
    lines!(rx.distance/1000, rad2deg.(p);
        linewidth=1.5)
    

    fig
    #save("LSIpath_segments_amp_phase.png", fig, px_per_unit=2)
end
