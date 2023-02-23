# iteratively exclude data from fit until fit error is below threshold
using LsqFit

# data from longwave/DispersionTest.jl
ωf = [2*pi*(6e3:1e3:18e3);];
ϕ = [   5.840479829632033
            3.7609869647037457
            8.513717721555789
            7.0863922989541015
           -0.35831763785328985
            5.149585430769938
            4.344839597842935
            9.638081754389809
            8.942412906068212
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

xyfit, xin, yin, xout, yout = iterfit(ωf, ϕ, thres);
fitcurve = xyfit.param[1].*ωf .+ xyfit.param[2] .+ xyfit.param[3]./(ωf);
fit_f₀ = (xyfit.param[3]*2*c/r)^(1/2)/(2*pi);

begin fig = Figure()
    fa = Axis(fig[1, 1], title=@sprintf("phase fit at r = %d km\nf₀ = %.2f kHz", r/1e3, fit_f₀/1e3),
        xlabel="frequency (kHz)",
        ylabel="phase (∘)")

    scatter!(fa, xin./(2*pi*1000), rad2deg.(yin);
        linewidth=2, color="black",
        label="in")

    scatter!(fa, xout./(2*pi*1000), rad2deg.(yout);
        linewidth=2, color="red",
        label="out")

    lines!(fa, ωf/(2*pi*1000), rad2deg.(fitcurve);
        linewidth=2, linestyle="-", color="black",
        label = "phase fit")

    axislegend()
    
    fig
end