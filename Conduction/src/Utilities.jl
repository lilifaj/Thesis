q = 1.6e-19; # Electron's charge (C)
k = 1.38e-23; # Boltzman constant (J.K^-1)
hbar = 1.054571817*10^(-34) # Planck constant in J.s

mutable struct Semiconductor
    alpha::Float64 # decay constant of the assumed hydrogen-like localized state wave functions (cm^-1)
    ModeEffect::Float64 # Mode effect of the phonons (J)
    Ni::Float64 # intrinsic semiconductor's density (cm^-3)
    Nd::Float64 # Doping states' density (cm^-3)
    Ed::Float64 # Energy to a vacant target site (J)
    F::Float64 # Field (V.cm^-1)
    nu::Float64 # Base electron jump rate
    Uf::Function # Fermi level (no unit)
    SigmaI::Function # Intrinsic semiconductor's gaussian width (J)
    SigmaD::Function # Doping states' gaussian width (J)
    beta::Function # Field Effect(No Unit)
    gamma::Function # Amount of disorder (J)
    omega_min::Real # Lower phonon frequency limit
    function Semiconductor(alpha, ModeEffect, Ni, Nd, Ed, F, nu, Uf::Real, SigmaI::Real, SigmaD::Real, Gamma::Real, omega_min::Real)
        FUf(T) = Uf * k * T
        FSigmaI(T) = SigmaI * k * T
        FSigmaD(T) = SigmaD * k * T
        Fbeta(T) = (F * q) / (2 * alpha * k * T)
        FGamma(T) = Gamma * k * T
        new(alpha, ModeEffect, Ni, Nd, Ed, F, nu, FUf, FSigmaI, FSigmaD, Fbeta, FGamma, omega_min)
    end
end

struct ConvergenceError <: Exception
end

Base.showerror(io::IO, e::ConvergenceError) = print(io, "Root function starting point does not lead to a convergence")

average_density_integral(f::Function, x::Real) = quadgk(
    x -> f(x),
    -x,
    x,
    rtol=1e-5,
)[1]

average_density_integral(f::Function, lower_value::Real, higher_value::Real) = quadgk(
    x -> f(x),
    lower_value,
    higher_value,
    rtol=1e-5,
)[1]

average_density(fn::Function, fd::Function, x::Real) =
    return average_density_integral(fn, x) / average_density_integral(fd, x)

average_density(fn::Function, fd::Function, lower_value::Real, higher_value::Real) =
    return average_density_integral(fn, lower_value, higher_value) / average_density_integral(fd, lower_value, higher_value)
