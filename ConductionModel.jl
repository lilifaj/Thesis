using QuadGK, HCubature, Roots, Base.Threads
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
    SigmaI::Function #  Intrinsic semiconductor's gaussian width (J)
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

# DOS of a doped semiconductor (J^-1.cm^-3)
DOS(semiconductor::Semiconductor, U, T) = (semiconductor.Ni / semiconductor.SigmaI(T) * exp(-(U * semiconductor.k * T - semiconductor.ModeEffect)^2 / (2 * semiconductor.SigmaI(T)^2)) + semiconductor.Nd / semiconductor.SigmaD(T) * exp(-(U * semiconductor.k * T - semiconductor.ModeEffect + semiconductor.Ed)^2 / (2 * semiconductor.SigmaD(T)^2))) / sqrt(2 * pi)

# Fermi-Dirac distribution
F(semiconductor::Semiconductor, U, T) = 1 / (1 + exp(U - (semiconductor.ModeEffect + semiconductor.Uf(T)) / (semiconductor.k * T)))

# Number of enclosed state by a radius of R
var1(U, beta, R, t, Rp, theta) = R + U - Rp * (1 + beta * cos(theta)) - t / (1 - t)

N(semiconductor::Semiconductor, U, T, R) = (semiconductor.k * T) / (8 * semiconductor.alpha^3) * 2 * pi * hcubature(
    x -> DOS(semiconductor, var1(U, semiconductor.beta(T), R, x[1], x[2], x[3]), T) * (1 - F(semiconductor, var1(U, semiconductor.beta(T), R, x[1], x[2], x[3]), T)) * 1 / (1 - x[1])^2 * x[2]^2 * sin(x[3]),
    [0, 0, 0],
    [1, R, pi],
    rtol=1e-5)[1]

RnnVRH(semiconductor::Semiconductor, U, T) = find_zero(r -> N(semiconductor, U, T, r) - 1, 5, Order0())

RnnPerco(semiconductor::Semiconductor, U, T) = ( (4pi) / (3 * 2.8) * quadgk(
    r -> DOS(semiconductor, r, T) * (1 - F(semiconductor, r, T)),
    -Inf,
    U + semiconductor.ModeEffect / (semiconductor.k * T),
    rtol=1e-8
    )[1])^(-1/3) * 2 * semiconductor.alpha * (semiconductor.k * T)^(-1/3)


RnnPercoField(semiconductor::Semiconductor, U, T) = find_zero(r -> N(semiconductor, U, T, r) - 2.8, 5, Order0())

var2(r, U, T, semiconductor, theta, Rnn) = Rnn * (r * (1 + semiconductor.beta(T) * cos(theta)) - semiconductor.beta(T) * cos(theta)) + U

var3(r, U, T, semiconductor, theta, Rnn) = U - Rnn * semiconductor.beta(T) * cos(theta) - r / (1 - r)

I1(U, T, semiconductor, Rnn::Float64) = 0.5 * Rnn * hcubature(
    x -> DOS(semiconductor, var2(x[1], U, T, semiconductor, x[2], Rnn), T) * (1 - F(semiconductor, var2(x[1], U, T, semiconductor, x[2], Rnn), T)) * sin(2 * x[2]) * (Rnn - var2(x[1], U, T, semiconductor, x[2], Rnn) + U)^3 / (1 + semiconductor.beta(T) * cos(x[2]))^2,
    [0, 0],
    [1, pi],
    rtol=1e-5,
    atol=1
)[1]

I2(U, T, semiconductor, Rnn::Float64) = 0.5 * hcubature(
    x -> DOS(semiconductor, var3(x[1], U, T, semiconductor, x[2], Rnn), T) * (1 - F(semiconductor, var3(x[1], U, T, semiconductor,x[2], Rnn), T)) * Rnn^3 * sin(2 * x[2]) / (1 - x[1])^2,
    [0, 0],
    [1, pi],
    rtol=1e-5,
    atol=1
)[1]

I3(U, T, semiconductor, Rnn::Float64) =Rnn * hcubature(
    x -> DOS(semiconductor, var2(x[1], U, T, semiconductor, x[2], Rnn), T) * (1 - F(semiconductor, var2(x[1], U, T, semiconductor, x[2], Rnn), T)) * sin(x[2]) * (Rnn - var2(x[1], U, T, semiconductor, x[2], Rnn) + U)^2 / (1 + semiconductor.beta(T) * cos(x[2])),
    [0, 0],
    [1, pi],
    rtol=1e-5,
    atol=1
)[1]

I4(U, T, semiconductor, Rnn::Float64) = hcubature(
    x -> DOS(semiconductor, var3(x[1], U, T, semiconductor, x[2], Rnn), T) * (1 - F(semiconductor, var3(x[1], U, T, semiconductor, x[2], Rnn), T)) * Rnn^2 * sin(x[2]) / (1 - x[1])^2,
    [0, 0],
    [1, pi],
    rtol=1e-5,
    atol=1
)[1]

function xf(semiconductor::Semiconductor, Rnn::Function, U, T)
    Rnn = Rnn(semiconductor, U, T)
    beta = semiconductor.beta(T)
    functionI = [I1, I2, I3, I4]
    resultI = Array{Float64}(undef, 4)

    Threads.@threads for i in 1:4
        resultI[i] = functionI[i](U, T, semiconductor, Rnn)
    end

    return (resultI[1] + resultI[2]) / (resultI[3] + resultI[4])
end

function asymptoteRnnPerco(semiconductor::Semiconductor, U, T)

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

function electronMobility(semiconductor::Semiconductor, Rnn::Function, U, T)
    return semiconductor.nu * xf(semiconductor, Rnn, U, T) * exp(-Rnn(semiconductor, U, T)) / (-semiconductor.F * 2 * semiconductor.alpha)
end