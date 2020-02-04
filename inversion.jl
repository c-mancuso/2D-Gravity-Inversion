#####################################################
#Julia 1.10 Written by C Mancuso March 2019
#2D Constrained Gravity Inversion
#Based on equations from Vatankhah, S., E. Ardestani,
#and M. Ashtari Jafari. (2014) "A method for 2-dimens
#ional inversion of gravity data."
#####################################################

using LinearAlgebra;
using NPZ;
using Distributions;
using Statistics;
using Random;

include("dataKernal.jl");       #forward model
include("conjugateGrad.jl");    #cg solver

#point to a 1D npy file of gravity measurments
gravity_data=npzread("bougeir.npy"); #reduced, check units!

#SET MODEL PARAMETERS
DEPTH = 30; LENGTH = 200; #number of voxels
z = 300;x = 300; #voxel size
model = (ones(DEPTH,LENGTH))*2;
background=2.5; #constant starting background density gcm3

elevation = zeros(LENGTH) #assuming elevation at zero

#SET PARAMETERS
ζ = 1e-4;
β = 1.00; ε = 1e-12;
ITERATIONS=4;

#***SHOULD NOT NEED TO EDIT BELOW THIS LINE***

#fix for reproducability
Random.seed!(44); println("Random Seed No. '",randstring(),"'\n") 

#get dimensions
M = size(model)[1] * size(model)[2];
N = size(model)[2];

#generate a G with assuming no change in datum height (z)
G = generate_Kernal(M, zeros(N), x, z);

g_obs = vcat(ones(20)*4,gravity_data);  #zero pad for stability

#make an H matrix "flatness"
H = Matrix{Float64}(I,M,M) * -1; #smoothness matrix
for i = 1:M-1, j = 2:M
    H[i,i+1] = 1;
    #H[j,j-1] = 1;
end

#make A matrix [G]
############## [H]

A = hcat(G', ζ * H')';

#make a P matrix ("hard constraint")
P = Matrix{Float64}(I, M, M)

#make a Q matrix ("depth weighting matrix")
Q = Matrix{Float32}(I, M, M);


for j = 1:M
    Q_jj = 1 / (((j*z) + ε)^β)
    Q[j,:] = Q[j,:] * Q_jj
end

#make a V matrix ("compactness")
𝓥 = Matrix{Float64}(I,M,M); #initially set to Identity

#set any initial density estimates
ρ_back = ones(M) * background
ρ = Array{Float64}(undef,M) + ρ_back; ρ_min = 2; ρ_max = 3;
null_vector = zeros(M)

for i = 1:M #update hard constraint
    if ρ[i] != 2.5
        P[i,i] = 1e-2
    end
end

for i = 1:ITERATIONS
    W = inv(P) * Q * 𝓥
    g_cal = G * ρ
    b = vcat(g_obs - g_cal, null_vector);
    Θ = return_m((A * inv(W)) * transpose((A * inv(W))), b, 40); #cg

    global ρ = ρ + (inv(W) * transpose(A * inv(W)) * Θ)

    for j = 1:M
        𝓥[j,j] = 1 / ((ρ[j])^2 + 1e-11)
        if ρ[j] < ρ_min || ρ[j] > ρ_max
            P[j,j] = 1e-2
        else
            P[j,j] = 1
        end
    end
end

model_est = reshape(ρ, (LENGTH, DEPTH));
npzwrite("me_chibo_mest.npy", model_est) #write out to a numpy file

#END
