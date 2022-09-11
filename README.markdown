# longwave

Testing out @fgasdia's [LongwaveModePropagator.jl](https://github.com/fgasdia/LongwaveModePropagator.jl).

## Error log
### [Introduction to defining scenarios](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/basic/)
No issues.
### [File-based I/O](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/io/)
Haven't done this one yet!

### [Mesh grid for mode finding part 1](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/meshgrid/), [part 2](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/meshgrid2/)
No issues.

### [Solvers for ionosphere reflection coefficient](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/integratedreflection/)
 1. Passing ````LMPParams()```` to ````LMP.dRdz```` doesn't work as written in the example:
    ```julia
      prob = ODEProblem{false}(LMP.dRdz, Rtop, (topheight, 0.0), (me, LMPParams()))
    ```
    Bypassed this by passing ````me```` as the only parameter:
    ```julia
      prob = ODEProblem{false}(LMP.dRdz, Rtop, (topheight, 0.0), me)
    ```
    Not sure what the impact on the ODE problem/solution is!
    
 2. ````solverstrings```` includes solver optional arguments, which makes some of them very long, which is difficult to plot.
    ```julia
      solverstrings = replace.(string.(solvers), "OrdinaryDiffEq."=>"")
    ```
    Solution: only use the substring before the first "(":
    ```julia
      solverstrings = string.(first.(split.(string.(solvers), "(")))
    ```
    This change improves the solver evaluation heatmap formatting.

### [Wavefield integration](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/wavefieldintegration/)
 1. Same problem and solution with ````solverstrings```` as above.
 2. My ````hx2_unscaled```` differs slightly from the example; this may be a precision error, so not really a problem.

### [Magnetic field direction](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/magneticfield/)
 1. ````Progress```` and ````next!```` are not found in ````LongwaveModePropagator````.  This is only important for progress logging in ````runlmp````.  To get around this, use:
    ```julia
      using ProgressLogging, TerminalLoggers
      using Logging: global_logger
    ```
    and make the following changes to ````runlmp````:
    ```julia
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
    ```
 2. The computed nighttime phase differs by a constant -20 dB offset from both the LWPC and example LMP runs.  Haven't been able to find the error yet.
 
 ## Versions:
 
 ```terminal
           Status ~\Project.toml
       [6e4b80f9] BenchmarkTools v1.3.1
       [336ed68f] CSV v0.10.4
       [13f3f980] CairoMakie v0.8.13
       [a216cea6] CompoundPeriods v0.5.1
       [0c46a032] DifferentialEquations v7.3.0
       [5789e2e9] FileIO v1.15.0
       [f67ccb44] HDF5 v0.16.11
       [a98d9a8b] Interpolations v0.14.5
       [033835bb] JLD2 v0.4.22
       [7f56f5a3] LSODA v0.7.0
       [38c6559a] LongwaveModePropagator v0.3.4
       [23992714] MAT v0.10.3
       [7eb4fadd] Match v1.2.0
       [1dea7af3] OrdinaryDiffEq v6.26.4
       [91a5bcdd] Plots v1.32.0
       [c3e4b0f8] Pluto v0.19.11
       [c46f51b8] ProfileView v1.5.1
       [33c8b6b6] ProgressLogging v0.1.4
       [a0859a10] RootsAndPoles v1.4.0
       [c3572dad] Sundials v4.10.1
       [5d786b92] TerminalLoggers v0.1.6
       [a759f4b9] TimerOutputs v0.5.21
       [ade2ca70] Dates
       [37e2e46d] LinearAlgebra
       [de0858da] Printf
 ```
