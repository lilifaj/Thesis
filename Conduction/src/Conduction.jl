module Conduction

using QuadGK, HCubature, Roots, Base.Threads
import Base.Threads.@spawn
include("Variables.jl")
include("Utilities.jl")

# DOS of a doped semiconductor (J^-1.cm^-3)
DOS(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = (semiconductor.Ni / semiconductor.SigmaI(T) * exp(-(U * k * T - semiconductor.ModeEffect)^2 / (2 * semiconductor.SigmaI(T)^2)) + semiconductor.Nd / semiconductor.SigmaD(T) * exp(-(U * k * T - semiconductor.ModeEffect + semiconductor.Ed)^2 / (2 * semiconductor.SigmaD(T)^2))) / sqrt(2 * pi)

DOSp2(semiconductor::Semiconductor, U::Real, T::Real) = 100^3 * (semiconductor.Ni) / semiconductor.SigmaI(T) * exp(-(log(U * k * T) - log(semiconductor.ModeEffect))^2 / (2 * (semiconductor.SigmaI(T) / hbar)^2)) / sqrt(2 * pi)

DOSp(semiconductor::Semiconductor, U::Real, T::Real) = 100^3 * (semiconductor.Ni) / semiconductor.SigmaI(T) * exp(-(U * k * T - semiconductor.ModeEffect)^2 / (2 * semiconductor.SigmaI(T)^2)) / sqrt(2 * pi)

# Fermi-Dirac distribution
F(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = 1 / (1 + exp(U - (semiconductor.ModeEffect + semiconductor.Uf(T)) / (k * T)))

# Number of free state within a sphere of radius R
N(semiconductor::Semiconductor, U::Real, T::Real, R::Real)::Float64 = (k * T) / (8 * semiconductor.alpha^3) * 2 * pi * hcubature(
    x -> DOS(semiconductor, var1(U, semiconductor.beta(T), R, x[1], x[2], x[3]), T) * (1 - F(semiconductor, var1(U, semiconductor.beta(T), R, x[1], x[2], x[3]), T)) * 1 / (1 - x[1])^2 * x[2]^2 * sin(x[3]),
    [0, 0, 0],
    [1, R, pi],
    rtol=1e-6)[1]

# Range to the nearest neighbour using a VRH approach
# function RnnVRH(semiconductor::Semiconductor, U::Real, T::Real)::Float64
#     i = 0;
#     while i < 10
#         try
#             return find_zero(r -> N(semiconductor, U, T, r) - 1, 5 + i * 10, Order0())
#             break
#         catch
#         end

#         i += 1;
#     end

#     throw(ConvergenceError())
# end

# Range to the nearest neighbour using a VRH approach
function RnnVRH(semiconductor::Semiconductor, U::Real, T::Real)::Float64
    Rnn(semiconductor, U, T, 1)
end

function Rnn(semiconductor::Semiconductor, U::Real, T::Real, x)::Float64
    i = 0;
    while i < 100
        try
            return find_zero(r -> N(semiconductor, U, T, r) - x, 5 + i * 10, Order0())
            break
        catch
        end

        i += 1;
    end

    throw(ConvergenceError())
end

# Range to the nearest neighbour using a percolation approach
RnnPerco(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = ( (4pi) / (3 * 2.8) * quadgk(
    r -> DOS(semiconductor, r, T) * (1 - F(semiconductor, r, T)),
    -Inf,
    U + semiconductor.ModeEffect / (k * T),
    rtol=1e-5
    )[1])^(-1/3) * 2 * semiconductor.alpha * (k * T)^(-1/3)

# Range to the nearest neighbour using a percolation approach taking into account the field
function RnnPercoField(semiconductor::Semiconductor, U::Real, T::Real)::Float64
    Rnn(semiconductor, U, T, 2.8)
end

# Effective distance of jump of an electron
function xf(semiconductor::Semiconductor, Rnn::Function, U::Real, T::Real)
    R = Rnn(semiconductor, U, T);
    functionI = [I1, I2, I3, I4]
    resultI = Array{Float64}(undef, 4)

    for i in 1:4
        resultI[i] = functionI[i](U, T, semiconductor, R)
    end

    return (resultI[1] + resultI[2]) / (resultI[3] + resultI[4])
end

function asymptoteRnnPerco(semiconductor::Semiconductor, U::Real, T::Real)::Float64

    positiveAsymptote = (2 * semiconductor.alpha) * (k * T * 4 * pi * quadgk(r-> DOS(semiconductor, r, T) * (1 - F(semiconductor, r, T)), -Inf, +Inf)[1] / (3 * 2.8))^(-1/3)
    x = -U - semiconductor.ModeEffect / (k * T)

    A1 = semiconductor.SigmaI(T)^2 / (x * (k * T)^2 + semiconductor.ModeEffect * k * T + semiconductor.SigmaI(T)^2) * semiconductor.Ni * exp(-semiconductor.ModeEffect^2 / (2 * semiconductor.SigmaI(T)^2) - (semiconductor.ModeEffect + semiconductor.Uf(T))/ (k * T)) / (sqrt(2pi) * semiconductor.SigmaI(T))

    A2 = semiconductor.SigmaD(T)^2 / (x * (k * T)^2 + (semiconductor.ModeEffect - semiconductor.Ed) * k * T + semiconductor.SigmaD(T)^2) * semiconductor.Nd * exp(-(semiconductor.ModeEffect - semiconductor.Ed)^2 / (2 * semiconductor.SigmaD(T)^2) - (semiconductor.ModeEffect + semiconductor.Uf(T))/ (k * T)) / (sqrt(2pi) * semiconductor.SigmaD(T))

    J1(x) = exp(-x^2 * (k * T)^2 / (2 * semiconductor.SigmaI(T)^2) - x * (1 + (semiconductor.ModeEffect * k * T) / (semiconductor.SigmaI(T)^2)))

    J2(x) = exp(-x^2 * (k * T)^2 / (2 * semiconductor.SigmaD(T)^2) - x * (1 + ((semiconductor.ModeEffect - semiconductor.Ed) * k * T) / (semiconductor.SigmaD(T)^2)))

    negativeAsymptote = 0
    try
        negativeAsymptote = (2 * semiconductor.alpha) * (k * T * 4 * pi * (A1 * J1(x) + A2 * J2(x)) / (3 * 2.8))^(-1/3)
    catch
    end

    (negativeAsymptote < positiveAsymptote ? res = positiveAsymptote : res = negativeAsymptote)

    return res
end

function electronMobility(semiconductor::Semiconductor, Rnn::Float64, xf::Float64)::Float64
    return semiconductor.nu * xf * exp(-Rnn) / (-semiconductor.F * 2 * semiconductor.alpha)
end

function t(semiconductor, Rnn::Float64, U::Real, T::Real)::Float64
    beta = semiconductor.beta(T)
    functionI = [It1, It2, It3, It4]
    resultI = Array{Float64}(undef, 4)

    for i in 1:4
        resultI[i] = functionI[i](U, T, semiconductor, Rnn)
    end

    return (resultI[1] + resultI[2]) / (resultI[3] + resultI[4])
end

# Electric diffusion (cm^2/s), first hypothesis
function D(semiconductor::Semiconductor, Rnn::Float64, xf::Float64, t::Float64)::Float64
    return xf^2 / (24 * semiconductor.alpha^2) * ((1 + semiconductor.nu * t * exp(-Rnn))^2 - 1) * semiconductor.nu * exp(-Rnn)
end

# Electric diffusion (cm^2/s), second hypothesis
function D_ter(semiconductor::Semiconductor, Rnn::Float64, xf::Float64, t::Float64)::Float64
    return 1 / 12 * (xf * Rnn * semiconductor.nu * exp(-Rnn) * t /  + 0.5 * Rnn^2) / semiconductor.alpha^2 * semiconductor.nu * exp(-Rnn)
end

function ein(semiconductor::Semiconductor, D::Function, Rnn::Float64, xf::Float64, t::Float64)::Float64
    return D(semiconductor, Rnn, xf, t) / electronMobility(semiconductor, Rnn, xf)
end

# function Dp(semiconductor::Semiconductor, T, eta, lower_value, higher_value)
#     fn(omega) = 1 / (eta * omega^2 * (exp(hbar * omega / (k * T)) - 1))
#     fd(omega) = omega^2 / (eta * semiconductor.gamma(T)^2 * (exp(hbar * omega / (k * T)) - 1))
#     return average_density(fn, fd, lower_value, higher_value)
# end

Dp(semiconductor, U, T) = semiconductor.gamma(T)^(-2) * (U * k * T/ hbar)^(-4)

# C(U, T) = exp(U / k / T) * (U / T)^2 / (k * (exp(U / k / T) - 1)^2)
C(U, T)::Float64 = U^2 * k / (exp(U/2) - exp(-U/2))^2

kp(semiconductor, T) = k * T * quadgk(
    r -> DOSp(semiconductor, r, T) * C(r, T) * Dp(semiconductor, r, T),
    semiconductor.omega_min * hbar / (k * T),
    +Inf
)[1]

function ke(semiconductor, T)
    function f(r)
        Rnn = Conduction.RnnVRH(semiconductor, r, T);
        xf = Conduction.xf(semiconductor, Rnn, r, T);
        t = Conduction.t(semiconductor, Rnn, r, T);
        return DOS(semiconductor, r, T) * C(r, t) * D_ter(semiconductor, Rnn, xf, t)
    end
    return k * T * quadgk(
    r -> f(r),
    -15,
    15)[1]
end

function occupiedStates(semiconductor::Semiconductor, U, T)
    return DOS(semiconductor, U, T) * F(semiconductor, U, T)
end

function overallDiffusion(semiconductor::Semiconductor, Rnn::Function, T, x_limit::Real)
    return Conduction.overallEinD(semiconductor, Conduction.D, Rnn, T, x_limit)
end

function overallMobility(semiconductor::Semiconductor, Rnn::Function, T, x_limit::Real)
    fd(x) = occupiedStates(semiconductor, x, T)

    function fn(x, Rnn)
        Rnn = Rnn(semiconductor, x, T);
        xf = Conduction.xf(semiconductor, Rnn, x, T);
        return occupiedStates(semiconductor, x, T) * electronMobility(semiconductor, Rnn, xf)
    end

    fn_final(x) = fn(x, Rnn);

    return average_density(fn_final, fd, x_limit);
end

function overallEinD(semiconductor::Semiconductor, f::Function, Rnn::Function, T, x_limit::Real)
    fd(x) = occupiedStates(semiconductor, x, T);

    function fn(x::Real, Rnn::Function)
        Rnn = Rnn(semiconductor, x, T);
        xf = Conduction.xf(semiconductor, Rnn, x, T);
        t = Conduction.t(semiconductor, Rnn, x, T);
        return occupiedStates(semiconductor, x, T) * f(semiconductor, Rnn, xf, t)
    end

    fn_final(x) = fn(x, Rnn);

    return average_density(fn_final, fd, x_limit)
end

function overallEin(semiconductor::Semiconductor, Rnn::Function, T, x_limit::Real)
    return overallEinD(semiconductor, (x, y, z, v) -> Conduction.ein(x, Conduction.D_ter, y, z, v), Rnn, T, x_limit)
end

end # module