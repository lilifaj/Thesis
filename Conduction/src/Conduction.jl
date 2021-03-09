module Conduction

using QuadGK, HCubature, Roots, Base.Threads
import Base.Threads.@spawn
include("Variables.jl")
include("Utilities.jl")

# DOS of a doped semiconductor (J^-1.cm^-3)
DOS(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = (semiconductor.Ni / semiconductor.SigmaI(T) * exp(-(U * semiconductor.k * T - semiconductor.ModeEffect)^2 / (2 * semiconductor.SigmaI(T)^2)) + semiconductor.Nd / semiconductor.SigmaD(T) * exp(-(U * semiconductor.k * T - semiconductor.ModeEffect + semiconductor.Ed)^2 / (2 * semiconductor.SigmaD(T)^2))) / sqrt(2 * pi)

# Fermi-Dirac distribution
F(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = 1 / (1 + exp(U - (semiconductor.ModeEffect + semiconductor.Uf(T)) / (semiconductor.k * T)))

# Number of free state within a sphere of radius R
N(semiconductor::Semiconductor, U::Real, T::Real, R::Real)::Float64 = (semiconductor.k * T) / (8 * semiconductor.alpha^3) * 2 * pi * hcubature(
    x -> DOS(semiconductor, var1(U, semiconductor.beta(T), R, x[1], x[2], x[3]), T) * (1 - F(semiconductor, var1(U, semiconductor.beta(T), R, x[1], x[2], x[3]), T)) * 1 / (1 - x[1])^2 * x[2]^2 * sin(x[3]),
    [0, 0, 0],
    [1, R, pi],
    rtol=1e-8)[1]

# Range to the nearest neighbour using a VRH approach
RnnVRH(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = find_zero(r -> N(semiconductor, U, T, r) - 1, 5, Order0())

# Range to the nearest neighbour using a percolation approach
RnnPerco(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = ( (4pi) / (3 * 2.8) * quadgk(
    r -> DOS(semiconductor, r, T) * (1 - F(semiconductor, r, T)),
    -Inf,
    U + semiconductor.ModeEffect / (semiconductor.k * T),
    rtol=1e-8
    )[1])^(-1/3) * 2 * semiconductor.alpha * (semiconductor.k * T)^(-1/3)

# Range to the nearest neighbour using a percolation approach taking into account the field
RnnPercoField(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = find_zero(r -> N(semiconductor, U, T, r) - 2.8, 5, Order0())

# Effective distance of jump of an electron
function xf(semiconductor::Semiconductor, Rnn::Float64, U::Real, T::Real)
    beta = semiconductor.beta(T)
    functionI = [I1, I2, I3, I4]
    resultI = Array{Float64}(undef, 4)

    for i in 1:4
        resultI[i] = functionI[i](U, T, semiconductor, Rnn)
    end

    return (resultI[1] + resultI[2]) / (resultI[3] + resultI[4])
end

function asymptoteRnnPerco(semiconductor::Semiconductor, U::Real, T::Real)::Float64

    positiveAsymptote = (2 * semiconductor.alpha) * (semiconductor.k * T * 4 * pi * quadgk(r-> DOS(semiconductor, r, T) * (1 - F(semiconductor, r, T)), -Inf, +Inf)[1] / (3 * 2.8))^(-1/3)
    x = -U - semiconductor.ModeEffect / (semiconductor.k * T)

    A1 = semiconductor.SigmaI(T)^2 / (x * (semiconductor.k * T)^2 + semiconductor.ModeEffect * semiconductor.k * T + semiconductor.SigmaI(T)^2) * semiconductor.Ni * exp(-semiconductor.ModeEffect^2 / (2 * semiconductor.SigmaI(T)^2) - (semiconductor.ModeEffect + semiconductor.Uf(T))/ (semiconductor.k * T)) / (sqrt(2pi) * semiconductor.SigmaI(T))

    A2 = semiconductor.SigmaD(T)^2 / (x * (semiconductor.k * T)^2 + (semiconductor.ModeEffect - semiconductor.Ed) * semiconductor.k * T + semiconductor.SigmaD(T)^2) * semiconductor.Nd * exp(-(semiconductor.ModeEffect - semiconductor.Ed)^2 / (2 * semiconductor.SigmaD(T)^2) - (semiconductor.ModeEffect + semiconductor.Uf(T))/ (semiconductor.k * T)) / (sqrt(2pi) * semiconductor.SigmaD(T))

    J1(x) = exp(-x^2 * (semiconductor.k * T)^2 / (2 * semiconductor.SigmaI(T)^2) - x * (1 + (semiconductor.ModeEffect * semiconductor.k * T) / (semiconductor.SigmaI(T)^2)))

    J2(x) = exp(-x^2 * (semiconductor.k * T)^2 / (2 * semiconductor.SigmaD(T)^2) - x * (1 + ((semiconductor.ModeEffect - semiconductor.Ed) * semiconductor.k * T) / (semiconductor.SigmaD(T)^2)))

    negativeAsymptote = 0
    try
        negativeAsymptote = (2 * semiconductor.alpha) * (semiconductor.k * T * 4 * pi * (A1 * J1(x) + A2 * J2(x)) / (3 * 2.8))^(-1/3)
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

function D(semiconductor::Semiconductor, Rnn::Float64, xf::Float64, t::Float64)::Float64
    return xf^2 / 6 * ((1 + semiconductor.nu * t * exp(-Rnn))^2 - 1) * semiconductor.nu * exp(-Rnn)
end

function D_no_t(semiconductor::Semiconductor, Rnn::Float64)
    return Rnn^2 / (6 * (2 * semiconductor.alpha)^2) * semiconductor.nu * exp(-Rnn)
end

end # module