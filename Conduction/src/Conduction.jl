module Conduction

using Base: Real, Float64
using QuadGK, HCubature, Roots, Base.Threads
import Base.Threads.@spawn
include("Variables.jl")
include("Utilities.jl")

# Please note that in this program, every energy given as a parameter of a function should be given in reduced units: u = U / k * T

# DOS of a doped semiconductor (J^-1.cm^-3)
DOS(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = (semiconductor.Ni / semiconductor.SigmaI * exp(-(U * k * T - semiconductor.ModeEffect)^2 / (2 * semiconductor.SigmaI^2)) + semiconductor.Nd / semiconductor.SigmaD * exp(-(U * k * T - semiconductor.ModeEffect + semiconductor.Ed)^2 / (2 * semiconductor.SigmaD^2))) / sqrt(2 * pi)

# Phonon DOS (J^-1 m^-1)
DOSp(semiconductor::Semiconductor, U::Real, T::Real) = 100^3 * (semiconductor.Ni) / semiconductor.SigmaI * exp(-(U * k * T - semiconductor.ModeEffect)^2 / (2 * semiconductor.SigmaI^2)) / sqrt(2 * pi)

# Fermi-Dirac distribution
F(semiconductor::Semiconductor, U::Real, T::Real)::Float64 = 1 / (1 + exp(U - (semiconductor.ModeEffect + semiconductor.Uf) / (k * T)))

# Number of free state within a sphere of radius R
N(semiconductor::Semiconductor, U::Real, T::Real, R::Real, F::Real)::Float64 = (k * T) / (8 * semiconductor.alpha^3) * 2 * pi * hcubature(
    x -> DOS(semiconductor, var1(U, Conduction.beta(semiconductor, T, F), R, x[1], x[2], x[3]), T) * (1 - Conduction.F(semiconductor, var1(U, Conduction.beta(semiconductor, T, F), R, x[1], x[2], x[3]), T)) * 1 / (1 - x[1])^2 * x[2]^2 * sin(x[3]),
    [0, 0, 0],
    [1, R, pi],
    rtol=1e-5)[1]

# Range to the nearest neighbour using a VRH approach (reduced unit)
function RnnVRH(semiconductor::Semiconductor, U::Real, T::Real, F::Real)::Float64
    Rnn(semiconductor, U, T, 1, F)
end

# Range to the nearest neighbour using a percolation approach taking into account the field (reduced unit)
function RnnPercoField(semiconductor::Semiconductor, U::Real, T::Real, F::Real)::Float64
    Rnn(semiconductor, U, T, 2.8, F)
end

function Rnn(semiconductor::Semiconductor, U::Real, T::Real, x::Real, F::Real)::Float64
    i = 0;
    while i < 100
        try
            return find_zero(r -> N(semiconductor, U, T, r, F) - x, 5 + i * 10, Order0())
            break
        catch
        end

        i += 1;
    end

    throw(ConvergenceError())
end

# Effective distance of jump of an electron (reduced unit)
function xf(semiconductor::Semiconductor, Rnn::Real, U::Real, T::Real, F::Real)
    functionI = [I1, I2, I3, I4]
    resultI = Array{Float64}(undef, 4)

    for i in 1:4
        resultI[i] = functionI[i](U, T, semiconductor, Rnn, F)
    end

    return (resultI[1] + resultI[2]) / (resultI[3] + resultI[4])
end

function xf(semiconductor::Semiconductor, U::Real, T::Real, F::Real)::Float64
    R = Conduction.RnnVRH(semiconductor, U, T, F);

    return xf(semiconductor, R, U, T, F)
end

# Return the mobility for the charge carriers (cm^2 V^-1 s^-1)
function mobility(semiconductor::Semiconductor, Rnn::Float64, xf::Float64, F::Real)::Float64
    return semiconductor.nu * xf * exp(-Rnn) / (-F * 2 * semiconductor.alpha)
end

function mobility(semiconductor::Semiconductor, U, T, F)::Float64
    R = Conduction.RnnVRH(semiconductor, U, T, F);
    xf = Conduction.xf(semiconductor, R, U, T, F);

    return mobility(semiconductor, R, xf, F)
end

# Stochastic time of trapping (s)
function t(semiconductor, Rnn::Float64, U::Real, T::Real, F::Real)::Float64
    functionI = [It1, It2, It3, It4]
    resultI = Array{Float64}(undef, 4)

    for i in 1:4
        resultI[i] = functionI[i](U, T, semiconductor, Rnn, F)
    end

    return (resultI[1] + resultI[2]) / (resultI[3] + resultI[4])
end

function t(semiconductor::Semiconductor, U, T, F)::Float64
    R = Conduction.RnnVRH(semiconductor, U, T, F);
    return t(semiconductor, R, U, T, F)
end

# Electric diffusion (cm^2 s^-!)
function D(semiconductor::Semiconductor, Rnn::Float64, xf::Float64, t::Float64)::Float64
    return 1 / 6 * (2 * xf * Rnn *  semiconductor.nu * exp(-Rnn) * t + Rnn^2) * semiconductor.nu * exp(-Rnn) / (4 * semiconductor.alpha^2)
end

function D(semiconductor::Semiconductor, U::Real, T::Real, F::Real)::Float64
    R = Conduction.RnnVRH(semiconductor, U, T, F);
    xf = Conduction.xf(semiconductor, R, U, T, F);
    t = Conduction.t(semiconductor, R, U, T, F);

    return D(semiconductor, R, xf, t)
end

# Einstein ratio (k_b T / q units)
function ein(semiconductor::Semiconductor, Rnn::Float64, xf::Float64, t::Float64, F::Real)::Float64
    return Conduction.D(semiconductor, Rnn, xf, t) / Conduction.mobility(semiconductor, Rnn, xf, F)
end

function ein(semiconductor::Semiconductor, U::Real, T::Real, F::Real)
    R = Conduction.RnnVRH(semiconductor, U, T, F);
    xf = Conduction.xf(semiconductor, R, U, T, F);
    t = Conduction.t(semiconductor, R, U, T, F);

    return ein(semiconductor, R, xf, t, F)
end

# Heat conduction by frequency (J K^-1)
C(U)::Float64 = U^2 * k / (exp(U/2) - exp(-U/2))^2

# Heat conduction for phonons (W m^-1 K^-1)
kp(semiconductor, T) = quadgk(
    r -> Conduction.k * T * Conduction.DOSp(semiconductor, r, T) * Conduction.C(r) * 4.10e-6,
    7.95e-21 / (Conduction.k * T),
    7.95e-20 / (Conduction.k * T)
)[1];

# Heat conduction for charge carriers (W m^-1 K^-1)
function ke(semiconductor, T, F)
    return 100 * k * T * (quadgk(
        r ->  DOS(semiconductor, r, T) * C(r) * D(semiconductor, r, T, F),
        0,
        15)[1] +
        quadgk(
        r ->  DOS(semiconductor, r, T) * C(r) * D(semiconductor, r, T, F),
        0,
        -15)[1])
end
end # module