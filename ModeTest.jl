# ModeTest.jl
# Todd Anderson
# 8 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 2. Mode finding

using Plots
using Plots.Measures
using RootsAndPoles

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, solvemodalequation, trianglemesh, defaultmesh

using Printf

## Part 1:
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/meshgrid2/

# define frequencies: 10 kHz, 20 kHz, 100 kHz
lowfrequency = Frequency(10e3)
midfrequency = Frequency(20e3)
highfrequency = Frequency(100e3)

# define day and night ionospheres
day = Species(QE, ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency)
night = Species(QE, ME, z->waitprofile(z, 85, 0.9), electroncollisionfrequency)

# waveguides
daywaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), day, GROUND[5])
nightwaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), night, GROUND[5])

# mode equations
day_low_me = PhysicalModeEquation(lowfrequency, daywaveguide)
day_mid_me = PhysicalModeEquation(midfrequency, daywaveguide)
day_high_me = PhysicalModeEquation(highfrequency, daywaveguide)
night_mid_me = PhysicalModeEquation(midfrequency, nightwaveguide)

# plot titles
day_low_title = "10 kHz\nh': 75 km, β: 0.35 km⁻¹"
day_mid_title = "20 kHz\nh': 75 km, β: 0.35 km⁻¹"
day_high_title = "100 kHz\nh': 75 km, β: 0.35 km⁻¹"
night_mid_title = "20 kHz\nh': 85 km, β: 0.9 km⁻¹"

# define rectangular mesh
Δr = 0.2
x = 0:Δr:90
y = -40:Δr:0
mesh = x .+ im*y';

function modeequationphase(me, mesh)
    phase = Vector{Float64}(undef, length(mesh))
    Threads.@threads for i in eachindex(mesh)
        f = solvemodalequation(deg2rad(mesh[i]), me)
        phase[i] = rad2deg(angle(f))
    end
    return phase
end

phase = modeequationphase(day_mid_me, mesh);

# plot mode equation phase θ for daytime ionosphere
heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(0, 90), ylims=(-40, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=day_mid_title,right_margin=2mm)

# zoom in on lowest-order modes
heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(30, 90), ylims=(-10, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=day_mid_title,right_margin=2mm)

# repeat for nighttime ionosphere
phase = modeequationphase(night_mid_me, mesh);
heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(0, 90), ylims=(-40, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=night_mid_title,right_margin=2mm)
        
# repeat for lower frequency
phase = modeequationphase(day_low_me, mesh);
heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(0, 90), ylims=(-40, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=day_low_title,right_margin=2mm)

# find roots and poles with the Global complex Roots and Poles Finding (GRPF) algorithm
# define mesh space
zbl = deg2rad(complex(30.0, -10.0)) # bottom left corner of complex domain
ztr = deg2rad(complex(89.9, 0.0)) # complex top right corner
Δr = deg2rad(0.5)

mesh = trianglemesh(zbl, ztr, Δr)

meshdeg = rad2deg.(mesh)

img = plot(real(meshdeg), imag(meshdeg); seriestype=:scatter,
        xlims=(80, 90), ylims=(-10, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        legend=false, size=(450,375));
plot!(img, [80, 90], [0, 0]; color="red");
plot!(img, [0, 90], [-90, 0]; color="red")

# daytime ionosphere, mid frequency
params = LMPParams().grpfparams
roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, day_mid_me),
                                                mesh, PlotData(), params);

z, edgecolors = getplotdata(tess, quads, phasediffs, g2f)

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
zdeg = rad2deg.(z)

twilightquads = [
        colorant"#9E3D36",
        colorant"#C99478",
        colorant"#7599C2",
        colorant"#5C389E",
        colorant"#404040",
        RGB(0.0, 0.0, 0.0)
]

img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth = 1.5,
            xlims=(30, 90), ylims=(-10, 0),
            xlabel="real(θ)", ylabel="imag(θ)", legend=false,
            title=day_mid_title);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
            seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
            seriestype=:scatter, markershape=:utriangle, markersize=5)

# daytime ionosphere, lower frequency
roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, day_low_me),
                                                mesh, PlotData(), params);

z, edgecolors = getplotdata(tess, quads, phasediffs, g2f)

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
zdeg = rad2deg.(z)

img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth = 1.5,
            xlims=(30, 90), ylims=(-10, 0),
            xlabel="real(θ)", ylabel="imag(θ)", legend=false,
            title=day_low_title);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
            seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
            seriestype=:scatter, markershape=:utriangle, markersize=5)

# daytime ionosphere, higher frequency
begin
    roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, day_high_me),
                                                    mesh, PlotData(), params);

    z, edgecolors = getplotdata(tess, quads, phasediffs, g2f)

    rootsdeg = rad2deg.(roots)
    polesdeg = rad2deg.(poles)
    zdeg = rad2deg.(z)

    img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth = 1.5,
                xlims=(30, 90), ylims=(-10, 0),
                xlabel="real(θ)", ylabel="imag(θ)", legend=false,
                title=day_high_title);
    plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
                seriestype=:scatter, markersize=5);
    plot!(img, real(polesdeg), imag(polesdeg); color="red",
                seriestype=:scatter, markershape=:utriangle, markersize=5)
end

# nighttime ionosphere, mid frequency
begin
    roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, night_mid_me),
                                                    mesh, PlotData(), params);

    z, edgecolors = getplotdata(tess, quads, phasediffs, g2f)

    rootsdeg = rad2deg.(roots)
    polesdeg = rad2deg.(poles)
    zdeg = rad2deg.(z)

    img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth = 1.5,
                xlims=(30, 90), ylims=(-10, 0),
                xlabel="real(θ)", ylabel="imag(θ)", legend=false,
                title=night_mid_title);
    plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
                seriestype=:scatter, markersize=5);
    plot!(img, real(polesdeg), imag(polesdeg); color="red",
                seriestype=:scatter, markershape=:utriangle, markersize=5)
end

## Part 2
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/meshgrid2/#Mesh-grid-for-mode-finding-Part-2
# Look at scenario were GRPF can fail to find all roots and poles

# segmented_scenario is known to miss roots if we use the same trianglemesh parameters
# (with Δr = deg2rad(0.5)) as in part 1 of this example.  Define the HomogeneousWaveguide
# from the second half of segmented_scenario known to have the missing roots
frequecy = Frequency(24e3)
electrons = Species(QE, ME, z->waitprofile(z, 80, 0.45), electroncollisionfrequency)
waveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), electrons, Ground(15, 0.001))
me = PhysicalModeEquation(frequecy, waveguide)
title = "24 kHz, h': 80 km, β: 0.45 km⁻¹"

# compute mode equation on a fine grid
Δr = 0.2
x = 30:Δr:90
y = -10:Δr:0
mesh = x .+ 1im*y';

# reuse modeequationphase() from part 1
phase = modeequationphase(me, mesh);

heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlim=(30, 90), ylims=(-10, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=title,
        right_margin=2mm)

# run grpf with Δr = 0.5
zbl = deg2rad(complex(30.0, -10.0))
ztr = deg2rad(complex(89.9, 0.0))
Δr = deg2rad(0.5)

mesh = trianglemesh(zbl, ztr, Δr)

params = LMPParams().grpfparams
roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, me), 
                                                mesh, PlotData(), params);

z, edgecolors = getplotdata(tess, quads, phasediffs, g2f)

zdeg = rad2deg.(z)
rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
rpstr = @sprintf("found %i roots and %i poles", length(roots), length(poles))
titlestr = string(title, "\n", rpstr)
# reuse twilightquads from Part 1

img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth = 1.5,
            xlims=(30, 90), ylims=(-10, 0),
            xlabel="real(θ)", ylabel="imag(θ)", legend=false,
            title=titlestr);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
            seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
            seriestype=:scatter, markershape=:utriangle, markersize=5)
# compare this plot with previous heatmap: root-pole pairs in upper right are missed by 
# grpf, because they are too closely-spaced.

# re-reun grpf with a finer mesh
Δr = deg2rad(0.2)

mesh = trianglemesh(zbl, ztr, Δr)

roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, me), 
                                                mesh, PlotData(), params);

z, edgecolors = getplotdata(tess, quads, phasediffs, g2f)

zdeg = rad2deg.(z)
rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
rpstr = @sprintf("found %i roots and %i poles", length(roots), length(poles))
titlestr = string(title, "\n", rpstr)

img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth = 1.5,
            xlims=(30, 90), ylims=(-10, 0),
            xlabel="real(θ)", ylabel="imag(θ)", legend=false,
            title=titlestr);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
            seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
            seriestype=:scatter, markershape=:utriangle, markersize=5)

# zoom in on upper right
img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth = 1.5,
            xlims=(80, 90), ylims=(-2, 0),
            xlabel="real(θ)", ylabel="imag(θ)", legend=false,
            title=titlestr);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
            seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
            seriestype=:scatter, markershape=:utriangle, markersize=5)

# plot LongwaveModePropagator.defaultmesh
# uses spacing of Δr = deg2rad(0.1) for upper right corner, Δr = deg2rad(0.5) otherwise

mesh = defaultmesh(frequecy)
meshdeg = rad2deg.(mesh)

img = plot(real(meshdeg), imag(meshdeg); seriestype=:scatter,
        xlabel="real(θ)", ylabel="imag(θ)",
        size=(450,375), legend=false);
plot!(img, [30, 90], [0, 0]; color="red");
plot!(img, [80, 90], [-10, 0]; color="red")

# solve modal equation on defaultmesh
roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, me), 
                                                mesh, PlotData(), params);

z, edgecolors = getplotdata(tess, quads, phasediffs, g2f)

zdeg = rad2deg.(z)
rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
rpstr = @sprintf("found %i roots and %i poles", length(roots), length(poles))
titlestr = string(title, "\n", rpstr)

img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth = 1.5,
            xlims=(30, 90), ylims=(-10, 0),
            xlabel="real(θ)", ylabel="imag(θ)", legend=false,
            title=titlestr);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
            seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
            seriestype=:scatter, markershape=:utriangle, markersize=5)
