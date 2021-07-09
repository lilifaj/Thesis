# Number of enclosed state by a radius of R
var1(U, beta, R, t, Rp, theta) = R + U - Rp * (1 + beta * cos(theta)) - t / (1 - t)

var2(r, U, T, semiconductor, theta, Rnn) = Rnn * (r * (1 + semiconductor.beta(T) * cos(theta)) - semiconductor.beta(T) * cos(theta)) + U

var3(r, U, T, semiconductor, theta, Rnn) = U - Rnn * semiconductor.beta(T) * cos(theta) - r / (1 - r)

I1(U, T, semiconductor, Rnn::Float64) = 0.5 * Rnn * hcubature(
    x -> DOS(semiconductor, var2(x[1], U, T, semiconductor, x[2], Rnn), T) * (1 - F(semiconductor, var2(x[1], U, T, semiconductor, x[2], Rnn), T)) * sin(2 * x[2]) * (Rnn - var2(x[1], U, T, semiconductor, x[2], Rnn) + U)^3 / (1 + semiconductor.beta(T) * cos(x[2]))^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I2(U, T, semiconductor, Rnn::Float64) = 0.5 * hcubature(
    x -> DOS(semiconductor, var3(x[1], U, T, semiconductor, x[2], Rnn), T) * (1 - F(semiconductor, var3(x[1], U, T, semiconductor,x[2], Rnn), T)) * Rnn^3 * sin(2 * x[2]) / (1 - x[1])^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I3(U, T, semiconductor, Rnn::Float64) =Rnn * hcubature(
    x -> DOS(semiconductor, var2(x[1], U, T, semiconductor, x[2], Rnn), T) * (1 - F(semiconductor, var2(x[1], U, T, semiconductor, x[2], Rnn), T)) * sin(x[2]) * (Rnn - var2(x[1], U, T, semiconductor, x[2], Rnn) + U)^2 / (1 + semiconductor.beta(T) * cos(x[2])),
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

I4(U, T, semiconductor, Rnn::Float64) = hcubature(
    x -> DOS(semiconductor, var3(x[1], U, T, semiconductor, x[2], Rnn), T) * (1 - F(semiconductor, var3(x[1], U, T, semiconductor, x[2], Rnn), T)) * Rnn^2 * sin(x[2]) / (1 - x[1])^2,
    [0, 0],
    [1, pi],
    rtol=1e-3,
    atol=1
)[1]

var4(V, r, theta, U, T, Rnn, semiconductor) = V * (Rnn - r) + U - r * semiconductor.beta(T) * cos(theta)

var5(V, r, theta, U, T, Rnn, semiconductor) = V / (V - 1) + U - r * semiconductor.beta(T) * cos(theta)

It1(U, T, semiconductor, Rnn::Float64) = hcubature(
    x -> (DOS(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, semiconductor), T) * (1 - F(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, semiconductor), T)) * exp((1 + semiconductor.beta(T) * cos(x[3])) * x[2] + var4(x[1], x[2], x[3], U, T, Rnn, semiconductor) - U) / semiconductor.nu * 2 * pi * x[2]^2 * sin(x[3])) * (Rnn - x[2]),
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It2(U, T, semiconductor, Rnn::Float64) = hcubature(
    x -> (DOS(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, semiconductor), T) * (1 - F(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, semiconductor), T)) * exp((1 + semiconductor.beta(T) * cos(x[3])) * x[2]) / semiconductor.nu * 2 * pi * x[2]^2 * sin(x[3])) / (1 - x[1])^2,
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It3(U, T, semiconductor, Rnn::Float64) = hcubature(
    x -> (DOS(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, semiconductor), T) * (1 - F(semiconductor, var4(x[1], x[2], x[3], U, T, Rnn, semiconductor), T))) * (Rnn - x[2]) * 2 * pi * x[2]^2 * sin(x[3]),
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

It4(U, T, semiconductor, Rnn::Float64) = hcubature(
    x -> (DOS(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, semiconductor), T) * (1 - F(semiconductor, var5(x[1], x[2], x[3], U, T, Rnn, semiconductor), T))) * 2 * pi * x[2]^2 * sin(x[3])/ (1 - x[1])^2,
    [0, 0, 0],
    [1, Rnn, pi],
    rtol=1e-5,
    atol=1
)[1]

ft11(U, T, semiconductor, Rnn::Float64, r, theta) = quadgk(
    x -> DOS(semiconductor, x, T) * (1 - F(semiconductor, x, T)) * exp((1 + semiconductor.beta(T) * cos(theta)) * r + x - U) / semiconductor.nu,
    U - r * semiconductor.beta(T) * cos(theta),
    Rnn + U - r * (1 + semiconductor.beta(T) * cos(theta)),
)[1]

ft12(U, T, semiconductor, Rnn::Float64, r, theta) = quadgk(
    x -> DOS(semiconductor, x, T) * (1 - F(semiconductor, x, T)) * exp((1 + semiconductor.beta(T) * cos(theta)) * r) / semiconductor.nu,
    -Inf,
    U - r * semiconductor.beta(T) * cos(theta),
)[1]

ft13(U, T, semiconductor, Rnn::Float64, r, theta) = quadgk(
    x -> DOS(semiconductor, x, T) * (1 - F(semiconductor, x, T)),
    U - r * semiconductor.beta(T) * cos(theta),
    Rnn + U - r * (1 + semiconductor.beta(T) * cos(theta)),
)[1]

ft14(U, T, semiconductor, Rnn::Float64, r, theta) = quadgk(
    x -> DOS(semiconductor, x, T) * (1 - F(semiconductor, x, T)),
    -Inf,
    U - r * semiconductor.beta(T) * cos(theta),
)[1]

ftt2(U, T, semiconductor, Rnn::Float64, theta, f) = quadgk(
    x -> f(U, T, semiconductor, Rnn, x, theta) * 2 * pi * x^2,
    0,
    Rnn,
    rtol=1e-5
)[1]

ftt3(U, T, semiconductor, f) = quadgk(
    x -> sin(x) * f(U, T, semiconductor, x),
    0,
    pi
)