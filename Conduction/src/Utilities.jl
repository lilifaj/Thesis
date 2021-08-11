q = 1.6e-19; # Electron's charge (C)
k = 1.38e-23; # Boltzman constant (J.K^-1)
hbar = 1.054571817*10^(-34) # Planck constant in J.s
ev = 1.602e-19; # Electron volt (J)

# Paraneters for a semiconductor
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
    function Semiconductor(alpha::Real, ModeEffect::Real, Ni::Real, Nd::Real, Ed::Real, nu::Real, Uf::Real, SigmaI::Real, SigmaD::Real) # This is the parameter to input when callon Semiconductor constructor
        Uf_2 = Uf * ev;
        Ed_2 = Ed * ev;
        ModeEffect_2 = ModeEffect * ev;
        SigmaI_2 = SigmaI * ev;
        SigmaD_2 = SigmaD * ev;
        new(alpha, ModeEffect_2, Ni, Nd, Ed_2, nu, Uf_2, SigmaI_2, SigmaD_2)
    end
end

# Created a new error when the convergence is not attained when findind a root
struct ConvergenceError <: Exception
end

Base.showerror(io::IO, e::ConvergenceError) = print(io, "Root function starting point does not lead to a convergence")

# Compute the integral of a function f between -x and x
average_density_integral(f::Function, x::Real)::Float64 = quadgk(
    x -> f(x),
    -x,
    x,
    rtol=1e-2,
)[1]

# Compute the integral of a function f between lower_value and higher_value
average_density_integral(f::Function, lower_value::Real, higher_value::Real)::Float64 = quadgk(
    x -> f(x),
    lower_value,
    higher_value,
    rtol=1e-2,
)[1]

# Compute the weighted average with an integral covering -x to x
average_density(fn::Function, fd::Function, x::Real)::Float64 =
    return average_density_integral(fn, x) / average_density_integral(fd, x)

# Compute the weighted average with an integral covering lower_value to higher_value
average_density(fn::Function, fd::Function, lower_value::Real, higher_value::Real)::Float64 =
    return average_density_integral(fn, lower_value, higher_value) / average_density_integral(fd, lower_value, higher_value)

# Express the field by the formula F * q / (2 alpha * k_b * T)
beta(semiconductor, T, F) = (F * q) / (2 * semiconductor.alpha * k * T)

# Compute the number of occupied state for a certain energy level U
function occupiedStates(semiconductor::Semiconductor, U::Real, T::Real)::Float64
    return DOS(semiconductor, U, T) * Conduction.F(semiconductor, U, T)
end

# Compute the global diffusion for a certain energy level U, field intensity F and temperature T between -x and x
function overallDiffusion(semiconductor::Semiconductor, T::Real, F::Real, x_limit::Real)::Float64
    return Conduction.overall(semiconductor, Conduction.D, T, F, x_limit)
end

# Compute the global mobility for a certain energy level U, field intensity F and temperature T between -x and x
function overallMobility(semiconductor::Semiconductor, T::Real, F::Real, x_limit::Real)::Float64
    return Conduction.overall(semiconductor, Conduction.mobility, T, F, x_limit)
end

# Compute the global einstein ratio for a certain energy level U, field intensity F and temperature T between -x and x
function overallEin(semiconductor::Semiconductor, T::Real, F::Real,  x_limit::Real)::Float64
    return overall(semiconductor, Conduction.ein, T, F, x_limit)
end

function overall(semiconductor::Semiconductor, f::Function, T, F::Real, x_limit::Real)::Float64
    fd(x::Real) = occupiedStates(semiconductor, x, T);
    fn(x::Real) = occupiedStates(semiconductor, x, T) * f(semiconductor, x, T, F)

    return average_density(fn, fd, x_limit)
end

end