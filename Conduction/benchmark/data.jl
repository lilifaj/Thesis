using Pkg
Pkg.activate(String(@__DIR__) * "/..")
Pkg.instantiate()

include("../src/Conduction.jl");
using Plots, VegaLite, LaTeXStrings, QuadGK, BenchmarkTools;

range = -10.:0.5:10.;

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

function einsteinRelation(semiconductor, U, T)
    r = Conduction.RnnVRH(semiconductor, U, T)
    x = Conduction.xf(semiconductor, r, U, T)
    t = Conduction.t(semiconductor, r, U, T)
    return Conduction.D(semiconductor, r, x, t) / Conduction.electronMobility(
        semiconductor, r, x)
end;

function einsteinRelation_no_t(semiconductor, U, T)
    r = Conduction.RnnVRH(semiconductor, U, T)
    x = Conduction.xf(semiconductor, r, U, T)
    return Conduction.D_no_t(semiconductor, r) / Conduction.electronMobility(
        semiconductor, r, x)
end

function criteria_D(semiconductor, U, T)
    r = Conduction.RnnVRH(semiconductor, U, T)
    return semiconductor.nu * exp(-r) * Conduction.t(semiconductor, r, U, T)
end;

function to_benchmark()
    pmap(x -> einsteinRelation(semiconductor, x, 300.), range);
end