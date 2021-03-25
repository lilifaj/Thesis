struct Semiconductor
    k::Float64 # Boltzman constant (J.K^-1)
    q::Float64 # Electron's charge (C)
    alpha::Float64 # decay constant of the assumed hydrogen-like localized state wave functions (cm^-1)
    ModeEffect::Float64 # Mode effect of the phonons (J)
    Ni::Float64 # intrinsic semiconductor's density (cm^-3)
    Nd::Float64 # Doping states' density (cm^-3)
    Ed::Float64 # Energy to a vacant target site (J)
    F::Float64 # Field (V.cm^-1)
    nu::Float64 # Base electron jump rate
    Uf::Function # Fermi level (J)
    SigmaI::Function # Intrinsic semiconductor's gaussian width (J)
    SigmaD::Function # Doping states' gaussian width (J)
    beta::Function # Field Effect(No Unit)
    function Semiconductor(k, q, alpha, ModeEffect, Ni, Nd, Ed, F, nu, Uf::Float64, SigmaI::Float64, SigmaD::Float64)
        FUf(T) = Uf * k * T
        FSigmaI(T) = SigmaI * k * T
        FSigmaD(T) = SigmaD * k * T
        Fbeta(T) = (F * q) / (2 * alpha * k * T)
        new(k, q, alpha, ModeEffect, Ni, Nd, Ed, F, nu, FUf, FSigmaI, FSigmaD, Fbeta)
    end
end

struct ConvergenceError <: Exception
end

Base.showerror(io::IO, e::ConvergenceError) = print(io, "Root function starting point does not lead to a convergence")

average_density_integral(f::Function, x) = quadgk(
    x -> f(x),
    -x,
    x,
    rtol=1e-5,
)[1]

average_density(fn::Function, fd::Function, x::Real) =
    return average_density(fn, x) / average_density(fd, x)
