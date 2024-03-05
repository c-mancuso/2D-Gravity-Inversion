# Julia 1.10
# Generate the Data Kernel "G" for the 2D Gravity Inversion
# Based on Potential Theory in Gravity and Magnetic Applications
# By Richard J. Blakely
# using NPZ
function returnCellularResponse(x_stat, z_stat, xcorn, zcorn)
    γ = 6.67408e-11; # gravitational constant in m3 kg-1 s-2
    Σ = 0;

    ncorn = Int16(size(xcorn)[1]) #get number of polygon corners

    for n = 1:ncorn
        if n == ncorn
            n2 = 1;
        else
            n2 = n + 1;
        end
        x1 = xcorn[n] - x_stat;
        z1 = zcorn[n] - z_stat;
        x2 = xcorn[n2] - x_stat;
        z2 = zcorn[n2] - z_stat;

        r1sq = x1^2 + z1^2;
        r2sq = x2^2 + z2^2;

        if r1sq == 0
            println("field point on a corner")
            return
        end
        if r1sq == 0
            println("field point on a corner")
            return
        end
        denom = z2 - z1
        if denom == 0
            denom = 1e-9;
        end

        α = (x2 - x1) / denom
        β = (x1 - α * z1)
        factor = β / (1.0 + α^2);
        term1 = (log(r2sq)) - (log(r1sq));
        term2 = atan(z2,x2) - atan(z1,x1);

        Σ = (Σ + (factor * (term1 - α * term2)));
    end

return (2 * γ * Σ * 1e8)
end

function generateKernel(M, z_stations, x_mult, z_mult)
    # Assumes that every column has a station in the center of it
    N = Int16(size(z_stations)[1]);
    z_max = Int16(M / N);

    # Set up a G
    G = zeros(N, M);

    # Make N numbered x by z model spaces
    model_space = Array{Float64}(zeros(N, z_max, N));
    print("(kernel.jl) the model is ", z_max, " by ", N);

    # Go through and calculate g_i of each cell within each model_space
    for i in 1:N             #all models; at each column
        for j in 1:N         #all columns
            for k in 1:z_max #all rows
                xcorners = [((j * x_mult) - x_mult), (j * x_mult), (j * x_mult), ((j * x_mult) -x_mult)];
                zcorners = [((k * z_mult) - z_mult), ((k * z_mult) - z_mult), (k * z_mult), (k * z_mult)];
                model_space[i,k,j] = returnCellularResponse((i * x_mult) - (x_mult / 2), z_stations[i], xcorners, zcorners);
            end
        end
        G[i,:] .= Iterators.flatten(transpose(model_space[i,:,:]));
    end
return G
end
