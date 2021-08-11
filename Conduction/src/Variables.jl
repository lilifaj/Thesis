# Variables for the various change of variable in integrals

# Change of variable for the number of free states N
var1(U::Real, beta::Real, R::Real, t::Real, Rp::Real, theta::Real)::Float64 = R + U - Rp * (1 + beta * cos(theta)) - t / (1 - t)

# Change of variable for real hopped distance xf
var2(semiconductor::Semiconductor, r::Real, U::Real, T::Real, F::Real, theta::Real, Rnn::Real)::Float64 = Rnn * (r * (1 + Conduction.beta(semiconductor, T, F) * cos(theta)) - Conduction.beta(semiconductor, T, F) * cos(theta)) + U

# Change of variable for real hopped distance xf
var3(semiconductor::Semiconductor, r::Real, U::Real, T::Real, F::Real, theta::Real, Rnn::Real)::Float64 = U - Rnn * Conduction.beta(semiconductor, T, F) * cos(theta) - r / (1 - r)

# Change of variable for the stochastic time of trap t
var4(V::Real, r::Real, theta::Real, U::Real, T::Real, Rnn::Real, F::Real, semiconductor::Semiconductor)::Float64 = V * (Rnn - r) + U - r * Conduction.beta(semiconductor, T, F) * cos(theta)

# Change of variable for the stochastic time of trap t
var5(V::Real, r::Real, theta::Real, U::Real, T::Real, Rnn, F::Real, semiconductor::Semiconductor)::Float64 = V / (V - 1) + U - r * Conduction.beta(semiconductor, T, F) * cos(theta)


# Integrals for the real hopped distabce xf
I1(U::Real, T::Real, semiconductor::Semiconductor, Rnn::Float64, F::Real)::Float64 = 0.5 * Rnn * hcubature(
    x -> DOS(semiconductor, var2(semiconductor, x[1], U, T, F, x[2], Rnn), T) * (1 - Conduction.F(semiconductor, var2(semiconductor, x[1], U, T, F, x[2], Rnn), T)) * sin(2 * x[2]) * (Rnn - var2(semiconductor, x[1], U, T, F, x[2], Rnn) + U)^3 / (1 + Conduction.beta(semiconductor, T, F) * cos(x[2]))^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I2(U::Real, T::Real, semiconductor::Semiconductor, Rnn::Float64, F::Real)::Float64 = 0.5 * hcubature(
    x -> DOS(semiconductor, var3(semiconductor, x[1], U, T, F, x[2], Rnn), T) * (1 - Conduction.F(semiconductor, var3(semiconductor, x[1], U, T, F,x[2], Rnn), T)) * Rnn^3 * sin(2 * x[2]) / (1 - x[1])^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I3(U::Real, T::Real, semiconductor::Semiconductor, Rnn::Float64, F::Real)::Float64 = Rnn * hcubature(
    x -> DOS(semiconductor, var2(semiconductor, x[1], U, T, F, x[2], Rnn), T) * (1 - Conduction.F(semiconductor, var2(semiconductor, x[1], U, T, F, x[2], Rnn), T)) * sin(x[2]) * (Rnn - var2(semiconductor, x[1], U, T, F, x[2], Rnn) + U)^2 / (1 + Conduction.beta(semiconductor, T, F) * cos(x[2])),
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I4(U::Real, T::Real, semiconductor::Semiconductor, Rnn::Float64, F::Real)::Float64 = hcubature(
    x -> DOS(semiconductor, var3(semiconductor, x[1], U, T, F, x[2], Rnn), T) * (1 - Conduction.F(semiconductor, var3(semiconductor, x[1], U, T, F, x[2], Rnn), T)) * Rnn^2 * sin(x[2]) / (1 - x[1])^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

# Integrals for the stochastic time of trap t
It1(U::Real, T::Real, semiconductor::Semiconductor, Rnn::Float64, F::Real)::Float64 = hcubature(
    x -> (DOS(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T) * (1 - Conduction.F(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T)) * exp((1 + Conduction.beta(semiconductor, T, F) * cos(x[3])) * x[2] + var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor) - U) / semiconductor.nu * 2 * pi * x[2]^2 * sin(x[3])) * (Rnn - x[2]),
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It2(U::Real, T::Real, semiconductor::Semiconductor, Rnn::Float64, F::Real)::Float64 = hcubature(
    x -> (DOS(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T) * (1 - Conduction.F(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T)) * exp((1 + Conduction.beta(semiconductor, T, F) * cos(x[3])) * x[2]) / semiconductor.nu * 2 * pi * x[2]^2 * sin(x[3])) / (1 - x[1])^2,
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It3(U::Real, T::Real, semiconductor::Semiconductor, Rnn::Float64, F::Real)::Float64 = hcubature(
    x -> (DOS(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T) * (1 - Conduction.F(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T))) * (Rnn - x[2]) * 2 * pi * x[2]^2 * sin(x[3]),
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It4(U::Real, T::Real, semiconductor::Semiconductor, Rnn::Float64, F::Real)::Float64 = hcubature(
    x -> (DOS(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T) * (1 - Conduction.F(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, F, semiconductor), T))) * 2 * pi * x[2]^2 * sin(x[3])/ (1 - x[1])^2,
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]
