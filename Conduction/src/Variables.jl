# Number of enclosed state by a radius of R
var1(U, beta, R, t, Rp, theta) = R + U - Rp * (1 + beta * cos(theta)) - t / (1 - t)

var2(semiconductor, r, U, T, F, theta, Rnn) = Rnn * (r * (1 + Conduction.beta(semiconductor, T, F) * cos(theta)) - Conduction.beta(semiconductor, T, F) * cos(theta)) + U

var3(semiconductor, r, U, T, F, theta, Rnn) = U - Rnn * Conduction.beta(semiconductor, T, F) * cos(theta) - r / (1 - r)

I1(U, T, semiconductor, Rnn::Float64, F::Real) = 0.5 * Rnn * hcubature(
    x -> DOS(semiconductor, var2(semiconductor, x[1], U, T, F, x[2], Rnn), T) * (1 - Conduction.F(semiconductor, var2(semiconductor, x[1], U, T, F, x[2], Rnn), T)) * sin(2 * x[2]) * (Rnn - var2(semiconductor, x[1], U, T, F, x[2], Rnn) + U)^3 / (1 + Conduction.beta(semiconductor, T, F) * cos(x[2]))^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I2(U, T, semiconductor, Rnn::Float64, F::Real) = 0.5 * hcubature(
    x -> DOS(semiconductor, var3(semiconductor, x[1], U, T, F, x[2], Rnn), T) * (1 - Conduction.F(semiconductor, var3(semiconductor, x[1], U, T, F,x[2], Rnn), T)) * Rnn^3 * sin(2 * x[2]) / (1 - x[1])^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I3(U, T, semiconductor, Rnn::Float64, F::Real) = Rnn * hcubature(
    x -> DOS(semiconductor, var2(semiconductor, x[1], U, T, F, x[2], Rnn), T) * (1 - Conduction.F(semiconductor, var2(semiconductor, x[1], U, T, F, x[2], Rnn), T)) * sin(x[2]) * (Rnn - var2(semiconductor, x[1], U, T, F, x[2], Rnn) + U)^2 / (1 + Conduction.beta(semiconductor, T, F) * cos(x[2])),
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I4(U, T, semiconductor, Rnn::Float64, F::Real) = hcubature(
    x -> DOS(semiconductor, var3(semiconductor, x[1], U, T, F, x[2], Rnn), T) * (1 - Conduction.F(semiconductor, var3(semiconductor, x[1], U, T, F, x[2], Rnn), T)) * Rnn^2 * sin(x[2]) / (1 - x[1])^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

var4(V, r, theta, U, T, Rnn, F::Real, semiconductor) = V * (Rnn - r) + U - r * Conduction.beta(semiconductor, T, F) * cos(theta)

var5(V, r, theta, U, T, Rnn, F::Real, semiconductor) = V / (V - 1) + U - r * Conduction.beta(semiconductor, T, F) * cos(theta)

It1(U, T, semiconductor, Rnn::Float64, F::Real) = hcubature(
    x -> (DOS(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T) * (1 - Conduction.F(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T)) * exp((1 + Conduction.beta(semiconductor, T, F) * cos(x[3])) * x[2] + var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor) - U) / semiconductor.nu * 2 * pi * x[2]^2 * sin(x[3])) * (Rnn - x[2]),
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It2(U, T, semiconductor, Rnn::Float64, F::Real) = hcubature(
    x -> (DOS(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T) * (1 - Conduction.F(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T)) * exp((1 + Conduction.beta(semiconductor, T, F) * cos(x[3])) * x[2]) / semiconductor.nu * 2 * pi * x[2]^2 * sin(x[3])) / (1 - x[1])^2,
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It3(U, T, semiconductor, Rnn::Float64, F::Real) = hcubature(
    x -> (DOS(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T) * (1 - Conduction.F(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T))) * (Rnn - x[2]) * 2 * pi * x[2]^2 * sin(x[3]),
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It4(U, T, semiconductor, Rnn::Float64, F::Real) = hcubature(
    x -> (DOS(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T) * (1 - Conduction.F(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T))) * 2 * pi * x[2]^2 * sin(x[3])/ (1 - x[1])^2,
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]
