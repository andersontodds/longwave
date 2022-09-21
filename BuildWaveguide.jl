# BuildWaveguide.jl
# Todd Anderson
# 17 September 2022
#
# Construct SegmentedWaveguide for any great circle propagation path based
# on land/sea/ice mask and day/night determination

using LongwaveModePropagator
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
    for i in 1:length(lon)
        meshlat = findmin(abs.(lat[i].-latmesh[:,1]))[2]
        meshlon = findmin(abs.(lon[i].-lonmesh[1,:]))[2]
        gnd[i] = mask[meshlat, meshlon]
    end
    return gnd
end

gnd1 = getground(sspath[1,:]', latmesh, lonmesh, LSI) # not sure why sspath[1,:] returns a column vector
pathgnd = getground(sspath, latmesh, lonmesh, LSI)

# plot LSI mask with Tx-Rx path, colored by LSI value
let fig = Figure(resolution = (1200,800))
    ga = GeoAxis(fig[1,1]; coastlines = true, title = "Land-sea-ice mask")

    surface!(ga, lonmesh, latmesh, LSI; 
        colormap="broc", colorrange=(-5,5), shading=false);
    scatter!(ga, reverse(Tx); color=:green)
    scatter!(ga, reverse(Rx); color=:red)
    lines!(ga, sspath; color=pathgnd, colormap="berlin", colorrange=(-2,2))

    fig
end

# find contiguous segments on geodesic that share ground values
# define function that takes sspath and pathgnd as input and outputs Mx2 matrix
# representing the length and ground value of each of M segments.  check that 
# the total length
function segments(path, ground)
    grounddiff = findall(abs.(diff(ground)).>0)
    segment_length = Vector{Float64}(undef, length(grounddiff)+1)
    segment_ground = Vector{Float64}(undef, length(grounddiff)+1)
    segment_ground[1] = ground[1]
    segment_length[1] = haversine(reverse(path[1,:]), reverse(path[grounddiff[1],:]), R_KM)
    for i in 2:length(segment_length)-1
        segment_ground[i] = ground[grounddiff[i]-1]
        segment_start = path[grounddiff[i-1],:]
        segment_end = path[grounddiff[i],:]
        segment_length[i] = haversine(reverse(segment_start), reverse(segment_end), R_KM)
    end
    segment_ground[end] = ground[end]
    segment_length[end] = haversine(reverse(path[grounddiff[end],:]), reverse(path[end,:]), R_KM)

    return segment_length, segment_ground
end

segment_length, segment_ground = segments(sspath, pathgnd)

# plot propagation path with waveguide segments
let fig = Figure(resolution = (1200,1200))
    ga = GeoAxis(fig[1,1]; coastlines = true, title = "Land-sea-ice mask",
        dest = "+proj=natearth", latlims = (40,90), lonlims = (-140, 30))

    surface!(ga, lonmesh, latmesh, LSI; 
        colormap="broc", colorrange=(-5,5), shading=false);
    scatter!(ga, reverse(Tx); color=:green)
    scatter!(ga, reverse(Rx); color=:red)
    lines!(ga, sspath; color=pathgnd, colormap="berlin", colorrange=(-2,2))
    

    fa = Axis(fig[2,1], title="waveguide segments")
    barplot!(fa, 1:length(segment_length), segment_length; 
        direction=:x, 
        color=segment_ground, colormap="berlin", colorrange=(-2,2))

    # make a legend!!!
    # labels = ["land", "sea", "ice"]
    # elements = [PolyElement(polycolor=)]

    fig

end

