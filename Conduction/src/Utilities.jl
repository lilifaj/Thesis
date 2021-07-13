q = 1.6e-19; # Electron's charge (C)
k = 1.38e-23; # Boltzman constant (J.K^-1)
hbar = 1.054571817*10^(-34) # Planck constant in J.s
ev = 1.602e-19; # Electron volt (J)
mutable struct Semiconductor
    alpha::Float64 # decay constant of the assumed hydrogen-like localized state wave functions (cm^-1)
    ModeEffect::Float64 # Mode effect of the phonons (J)
    Ni::Float64 # intrinsic semiconductor's density (cm^-3)
    Nd::Float64 # Doping states' density (cm^-3)
    Ed::Float64 # Energy to a vacant target site (J)
    nu::Float64 # Base electron jump rate
    Uf::Float64 # Fermi level (J)
    SigmaI::Float64 # Intrinsic semiconductor's gaussian width (J)
    SigmaD::Float64 # Doping states' gaussian width (J)
    function Semiconductor(alpha, ModeEffect, Ni, Nd, Ed, nu, Uf::Real, SigmaI::Real, SigmaD::Real)
        Uf_2 = Uf * ev;
        SigmaI_2 = SigmaI * ev;
        SigmaD_2 = SigmaD * ev;
        new(alpha, ModeEffect, Ni, Nd, Ed, nu, Uf_2, SigmaI_2, SigmaD_2)
    end
end

struct ConvergenceError <: Exception
end

Base.showerror(io::IO, e::ConvergenceError) = print(io, "Root function starting point does not lead to a convergence")

average_density_integral(f::Function, x::Real) = quadgk(
    x -> f(x),
    -x,
    x,
    rtol=1e-2,
)[1]

average_density_integral(f::Function, lower_value::Real, higher_value::Real) = quadgk(
    x -> f(x),
    lower_value,
    higher_value,
    rtol=1e-2,
)[1]

average_density(fn::Function, fd::Function, x::Real) =
    return average_density_integral(fn, x) / average_density_integral(fd, x)

average_density(fn::Function, fd::Function, lower_value::Real, higher_value::Real) =
    return average_density_integral(fn, lower_value, higher_value) / average_density_integral(fd, lower_value, higher_value)

beta(semiconductor, T, F) = (F * q) / (2 * semiconductor.alpha * k * T)

function asymptoteRnnPerco(semiconductor::Semiconductor, U::Real, T::Real)::Float64

    positiveAsymptote = (2 * semiconductor.alpha) * (k * T * 4 * pi * quadgk(r-> DOS(semiconductor, r, T) * (1 - F(semiconductor, r, T)), -Inf, +Inf)[1] / (3 * 2.8))^(-1/3)
    x = -U - semiconductor.ModeEffect / (k * T)

    A1 = semiconductor.SigmaI^2 / (x * (k * T)^2 + semiconductor.ModeEffect * k * T + semiconductor.SigmaI^2) * semiconductor.Ni * exp(-semiconductor.ModeEffect^2 / (2 * semiconductor.SigmaI^2) - (semiconductor.ModeEffect + semiconductor.Uf)/ (k * T)) / (sqrt(2pi) * semiconductor.SigmaI)

    A2 = semiconductor.SigmaD^2 / (x * (k * T)^2 + (semiconductor.ModeEffect - semiconductor.Ed) * k * T + semiconductor.SigmaD^2) * semiconductor.Nd * exp(-(semiconductor.ModeEffect - semiconductor.Ed)^2 / (2 * semiconductor.SigmaD^2) - (semiconductor.ModeEffect + semiconductor.Uf)/ (k * T)) / (sqrt(2pi) * semiconductor.SigmaD)

    J1(x) = exp(-x^2 * (k * T)^2 / (2 * semiconductor.SigmaI^2) - x * (1 + (semiconductor.ModeEffect * k * T) / (semiconductor.SigmaI^2)))

    J2(x) = exp(-x^2 * (k * T)^2 / (2 * semiconductor.SigmaD^2) - x * (1 + ((semiconductor.ModeEffect - semiconductor.Ed) * k * T) / (semiconductor.SigmaD^2)))

    negativeAsymptote = 0
    try
        negativeAsymptote = (2 * semiconductor.alpha) * (k * T * 4 * pi * (A1 * J1(x) + A2 * J2(x)) / (3 * 2.8))^(-1/3)
    catch
    end

    (negativeAsymptote < positiveAsymptote ? res = positiveAsymptote : res = negativeAsymptote)

    return res
end