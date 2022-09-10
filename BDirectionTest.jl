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
using ProgressLogging: Progress, next!
# Progress and next! are not found in LongwaveModePropagator.  Guessing Progress 
# is from ProgessLogging, but I'm not sure where next! is.  Look into how to emulate 
# what this example is trying to do with ProgressLogging.  Also look into where buildrun 
# is defined in LongwaveModePropagator.

using LongwaveModePropagator
using LongwaveModePropagator: buildrun, next! # Progress not defined in LongwaveModePropagator

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

        input.name = @sprintf("%d_.0f%_.0f", i, b_dips[i], b_azs[i])
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


