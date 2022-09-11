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
    ```
      prob = ODEProblem{false}(LMP.dRdz, Rtop, (topheight, 0.0), (me, LMPParams()))
    ```
    Bypassed this by passing ````me```` as the only parameter:
    ```
      prob = ODEProblem{false}(LMP.dRdz, Rtop, (topheight, 0.0), me)
    ```
    Not sure what the impact on the ODE problem/solution is!
    
 2. ````solverstrings```` includes solver optional arguments, which makes some of them very long, which is difficult to plot.
    ```
      solverstrings = replace.(string.(solvers), "OrdinaryDiffEq."=>"")
    ```
    Solution: only use the substring before the first "(":
    ```
      solverstrings = string.(first.(split.(string.(solvers), "(")))
    ```
    This change improves the solver evaluation heatmap formatting.

### [Wavefield integration](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/wavefieldintegration/)
 1. Same problem and solution with ````solverstrings```` as above.
 2. My ````hx2_unscaled```` differs slightly from the example; this may be a precision error, so not really a problem.

### [Magnetic field direction](https://fgasdia.github.io/LongwaveModePropagator.jl/dev/generated/magneticfield/)
 1. ````Progress```` and ````next!```` are not found in ````LongwaveModePropagator````.  This is only important for progress logging in ````runlmp````.  To get around this, use:
    ```
      using ProgressLogging, TerminalLoggers
      using Logging: global_logger
    ```
    and make the following changes to ````runlmp````:
    ```
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
    
