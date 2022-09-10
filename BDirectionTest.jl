# BDirectionTest.jl
# Todd Anderson
# 10 September 2022
#
# Explore examples in LongwaveModePropagator.jl by @fgasdia
# 5. Magnetic field direction
# https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/magneticfield/
#
# Explore the influence of magnetic field direction on daytime and nighttime propagation
# path; compare results to LWPC

## Generate scenarios

using Dates, Printf
using HDF5
using Plots
using ProgressLogging, TerminalLoggers
using Logging: global_logger
# Progress and next! are not found in LongwaveModePropagator.  Guessing Progress 
# is from ProgessLogging, but I'm not sure where next! is.  Look into how to emulate 
# what this example is trying to do with ProgressLogging.  Also look into where buildrun 
# is defined in LongwaveModePropagator.

using LongwaveModePropagator
using LongwaveModePropagator: buildrun # Progress and next! not defined in LongwaveModePropagator

# define global OUTPUT_RANGES at which electric field will be computed, and magnetic field
# dip angles and azimuth angles B_DIPS and B_AZS

const OUTPUT_RANGES = 0:5e3:3000e3  # m
const B_DIPS = [90.0, 60, 60, 60, 60, 60, 60, 60, 60] # degrees
const B_AZS = [0.0, 0, 45, 90, 135, 180, 225, 270, 315] # degrees: vertical; N, E, S, W

# write function to generate BatchInput of scenarios

function generate(hp, β)
    batch = BatchInput{ExponentialInput}()
    batch.name = "Magnetic Field Tests"
    batch.description = "Varying magnetic field directions: vertical, N, E, S, W."
    batch.datetime = Dates.now()

    # Constants
    frequency = 24e3
    output_ranges = collect(OUTPUT_RANGES)
    segment_ranges = [0.0]
    hprime = [hp]
    beta = [β]

    # Ocean
    ground_epsr = [81]
    ground_sigma = [4.0]

    b_mag = [50e-6]
    b_dips = deg2rad.(B_DIPS)
    b_azs = deg2rad.(B_AZS)

    N = length(b_azs)
    inputs = Vector{ExponentialInput}(undef, N)
    for i in 1:N
        input = ExponentialInput()

        input.name = @sprintf("%d_%.0f_%.0f", i, b_dips[i], b_azs[i])
        input.description = "Wait ionosphere with ocean ground at 24 kHz."
        input.datetime = Dates.now()
        input.segment_ranges = segment_ranges
        input.hprimes = hprime
        input.betas = beta
        input.b_mags = b_mag
        input.b_dips = [b_dips[i]]
        input.b_azs = [b_azs[i]]
        input.ground_sigmas = ground_sigma
        input.ground_epsrs = ground_epsr
        input.frequency = frequency
        input.output_ranges = output_ranges

        inputs[i] = input
    end
    batch.inputs = inputs

    # json_str = JSON3.write(batch)
    # 
    # open("bfields.json","w") do f
    #    write(f, json_str)
    # end

    return batch
end

batch = generate(82.0, 0.6)

# # test progress log
# N = length(B_AZS)
# a = Vector{Float64}(undef, N)

# global_logger(TerminalLogger(right_justify=120))
# @progress for i in 1:length(B_AZS)
#     a[i] = B_AZS[i] + rand()
#     sleep(1)
#     @info "Middle of computation"
# end
# @info "Done"

## Run the model

function runlmp(inputs, outfile)

    h5open(outfile, "cw") do fid
        # Create batch attributes if they don't already exist
        fid_attrs = attributes(fid)
        haskey(fid_attrs, "name") || (fid_attrs["name"] = inputs.name)
        haskey(fid_attrs, "description") || (fid_attrs["description"] = inputs.description)
        haskey(fid_attrs, "datetime") || (fid_attrs["datetime"] = string(Dates.now()))

        if haskey(fid, "outputs")
            g = fid["outputs"]
        else
            g = create_group(fid, "outputs")
        end
    end

    global_logger(TerminalLogger(right_justify=120))

    @progress for i in eachindex(inputs.inputs)
        # If we've already run this, we can skip ahead
        complete = h5open(outfile, "r") do fid
            g = open_group(fid, "outputs")
            haskey(g, string(i)) ? true : false
        end
        if complete 
            @info @sprintf("Already run %i!", i)
            continue; 
        end

        output = buildrun(inputs.inputs[i])

        h5open(outfile, "r+") do fid
            g = open_group(fid, "outputs")
            o = create_group(g, string(i))
            attributes(o)["name"] = output.name
            attributes(o)["description"] = output.description
            attributes(o)["datetime"] = string(Dates.now())
            o["output_ranges"] = output.output_ranges
            o["amplitude"] = output.amplitude
            o["phase"] = output.phase
        end
        @info @sprintf("Done with %i, next!", i)
    end
end

# file path from example:
# root_dir = dirname(dirname(pathof(LongwaveModePropagator)))
# examples_dir = joinpath(root_dir, "examples")
# lmpfile = joinpath(examples_dir, "bfields_lmp.h5")

out_dir = joinpath(pwd(), "output_files/")
lmpfile = joinpath(out_dir, "bfields_lmp.h5")

# run!
runlmp(batch, lmpfile)

## Plots
# process HDF5 files into arrays of amplitude and phase v. range and magnetic field azimuth
function process(outputs)
    dist = read(outputs["outputs"]["1"]["output_ranges"])
    mask = indexin(OUTPUT_RANGES, dist)

    agrid = Array{Float64}(undef, length(OUTPUT_RANGES), length(B_AZS))
    pgrid = similar(agrid)

    for i in eachindex(B_AZS)
        o = outputs["outputs"][string(i)]

        amp = read(o["amplitude"])
        phase = read(o["phase"])

        agrid[:,i] = amp[mask]
        pgrid[:,i] = rad2deg.(phase[mask])
    end

    return agrid, pgrid
end


agrid, pgrid = h5open(lmpfile, "r") do o
    agrid, pgrid = process(o)
end


labels = string.(trunc.(Int, B_AZS), "°")
labels[1] = "90°, "*labels[1]
labels[2] = "60°, "*labels[2]
labels[3:end] .= "      ".*labels[3:end]

colors = [palette(:phase, length(B_AZS))...]
pushfirst!(colors, RGB(0.0, 0, 0))

plot(OUTPUT_RANGES/1000, agrid;
    linewidth=1.5, palette=colors, colorbar=false, 
    xlabel="range (km)", ylabel="amplitude (dB)", 
    labels=permutedims(labels), legendtitle="  dip, az", legend=true)


# compare with LWPC

root_dir = dirname(dirname(pathof(LongwaveModePropagator)))
examples_dir = joinpath(root_dir, "examples")
lwpcfile = joinpath(examples_dir, "bfields_lwpc.h5")

lagrid, lpgrid = h5open(lwpcfile, "r") do o
    lagrid, lpgrid = process(o)
end

plot(OUTPUT_RANGES/1000, lagrid;
    linewidth=1.5, palette=colors, colorbar=false, 
    xlabel="range (km)", ylabel="amplitude (dB)", 
    labels=permutedims(labels), legendtitle="  dip, az", legend=true)

adifference = agrid - lagrid

plot(OUTPUT_RANGES/1000, adifference;
     linewidth=1.5, palette=colors, colorbar=false,
     xlabel="range (km)", ylabel="amplitude difference (dB)",
     labels=permutedims(labels), legendtitle="  dip, az", legend=true)

# compare with example LMP file

ex_lmpfile = joinpath(examples_dir, "bfields_lmp.h5")
eagrid, epgrid = h5open(ex_lmpfile, "r") do o
    eagrid, epgrid = process(o)
end

adifference = agrid - eagrid

plot(OUTPUT_RANGES/1000, adifference;
     linewidth=1.5, palette=colors, colorbar=false,
     xlabel="range (km)", ylabel="amplitude difference (dB)",
     labels=permutedims(labels), legendtitle="  dip, az", legend=true)

## Daytime ionosphere
# Does magnetic field direction have a different level of influence on daytime ionospheres?

lmpfile = joinpath(out_dir, "bfields_daytime_lmp.h5")
runlmp(generate(72.0, 0.3), lmpfile)

agrid, pgrid = h5open(lmpfile, "r") do o
    agrid, pgrid = process(o)
end

plot(OUTPUT_RANGES/1000, agrid;
     linewidth=1.5, palette=colors, colorbar=false,
     xlabel="range (km)", ylabel="amplitude (dB)",
     labels=permutedims(labels), legendtitle="  dip, az", legend=true)
#savefig("magneticfield_daytime_amplitude.png"); nothing # hide
