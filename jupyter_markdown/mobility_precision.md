```julia
include("../src/Conduction.jl");
using Plots, LaTeXStrings;
```


```julia
semiconductor = Conduction.Semiconductor(
1.38 * 10^-23, # Boltzman constant (J.K^-1)
1.6*10^-19, # Electron's charge (C)
10^7, # decay constant of the assumed hydrogen-like localized state wave functions (cm^-1)
0.1 * 1.6*10^-19, # Mode effect of the phonons (J)
2.1 * 10^18, # intrinsic semiconductor's density (cm^-3)
2.1 * 10^18, # Doping states' density (cm^-3)
0.1 * 1.6 * 10^-19, # Energy to a vacant target site (J)
-2*10^5, # Field (V.cm^-1)
10^13, # Base electron jump rate
-10.0, # Fermi level (J)
2.7, # Intrinsic semiconductor's gaussian width (J)
2.7 # Doping states' gaussian width (J)
);

function fd(semiconductor, U, T)
    return Conduction.DOS(semiconductor, U, T) * (1 - Conduction.F(semiconductor, U, T))
end

function fn(semiconductor, U, T)
    Rnn = Conduction.RnnVRH(semiconductor, U, T)
    xf = Conduction.xf(semiconductor, Rnn, U, T)
    t = Conduction.t(semiconductor, Rnn, U, T)
    return Conduction.ein(semiconductor, Rnn, xf, t) * fd(semiconductor, U, T)
end
```




    fn (generic function with 1 method)




```julia
T = 300;

Conduction.average_density(x -> fn(semiconductor, x, T), x -> fd(semiconductor, x, T), 10)
```




    6.117974971777466e11




```julia
Conduction.average_density(x -> fn(semiconductor, x, T), x -> fd(semiconductor, x, T), 15)
```




    6.082848630281743e11




```julia
Conduction.average_density(x -> fn(semiconductor, x, T), x -> fd(semiconductor, x, T), 16)
```




    6.082801453968109e11




```julia
res = Conduction.average_density(x -> fn(semiconductor, x, T), x -> fd(semiconductor, x, T), 20)
```




    6.082804191878281e11




```julia
res * semiconductor.k * T / semiconductor.q
```




    1.5739255846485054e10


