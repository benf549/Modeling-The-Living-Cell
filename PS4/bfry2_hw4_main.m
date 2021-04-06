%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 03/15/21, coded on MATLAB _R2020b_ 
% 
%Problem Set 4 - Modeling The Living Cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 2 - Molecular Dynamics
positions = load("init_crds_boxl_3.5.dat");

T_star = 1;
sigma_star = 1;
epsilon_star = 0.25;
density_star = 0.7;
L = 3.5;
time_step = 0.005;
maxtime = 50;

mass_star = density_star * L^3 / 30;


%Equilibrate for 50 time units
[o, eq_positions] = MolecularDynamics(positions, sigma_star, epsilon_star, T_star, mass_star, L, time_step, maxtime);
%Sample for 50 additional time units
[output, f] = MolecularDynamics(eq_positions, sigma_star, epsilon_star, T_star, mass_star, L, time_step, maxtime);
%Visualize sampled energies and temperatures
PlotMDOutput(output, time_step, maxtime, 1, "T_{init} = 1, \Deltat =  0.005")


%Calculate average potential and temperature from output
avg_potential = sum(output(1, :)) / length(output(1,:));
fprintf("MD Avg. Potential Energy w/ initial T=%.2f as: %.4f\n", T_star, avg_potential);
avg_temp = sum(output(3, :)) / length(output(3, :));
fprintf("MD Avg. Temperature w/ initial T=%.2f as: %.4f\n\n", T_star, avg_temp)


%Repeat for T = 0.5 and T = 0.
T_star = 0.5;
[o, eq_positions] = MolecularDynamics(positions, sigma_star, epsilon_star, T_star, mass_star, L, time_step, maxtime);
[output, f] = MolecularDynamics(eq_positions, sigma_star, epsilon_star, T_star, mass_star, L, time_step, maxtime);
PlotMDOutput(output, time_step, maxtime, 3, "T_{init} = 0.5, \Deltat =  0.005")
avg_potential = sum(output(1, :)) / length(output(1,:));
fprintf("MD Avg. Potential Energy w/ initial T=%.2f as: %.4f\n", T_star, avg_potential);
avg_temp = sum(output(3, :)) / length(output(3, :));
fprintf("MD Avg. Temperature w/ initial T=%.2f as: %.4f\n\n", T_star, avg_temp)

T_star = 0.1;
[o, eq_positions] = MolecularDynamics(positions, sigma_star, epsilon_star, T_star, mass_star, L, time_step, maxtime);
[output, f] = MolecularDynamics(eq_positions, sigma_star, epsilon_star, T_star, mass_star, L, time_step, maxtime);
PlotMDOutput(output, time_step, maxtime, 5, "T_{init} = 0.1, \Deltat =  0.005")
avg_potential = sum(output(1, :)) / length(output(1,:));
fprintf("MD Avg. Potential Energy w/ initial T=%.2f as: %.4f\n", T_star, avg_potential);
avg_temp = sum(output(3, :)) / length(output(3, :));
fprintf("MD Avg. Temperature w/ initial T=%.2f as: %.4f\n\n", T_star, avg_temp)


% Energy is not constant when the time step is too large.
T_star = 1;
[output2, f] = MolecularDynamics(positions, sigma_star, epsilon_star, T_star, mass_star, L, .1, 10);
PlotMDOutput(output2, .1, 10, 7, "T_{init} = 1, \Deltat =  0.1")

%% Problem 3 - Brownian Dynamics
kT = 10;
dt = 1e-4;
D = 1;
num_loops = 50001;
lambda = (cot(70/180*pi))^2;

%Run brownian dynamics simulation for kT = 10, 5, 1, and 0.5 and plot
[coords_kT10, x] = BrownianDynamics(kT, D, dt, lambda, num_loops, true, "kT = 10", 1, false); %#ok<*ASGLU>
[coords_kT5, x] = BrownianDynamics(5, D, dt, lambda, num_loops, true, "kT = 5", 2, false);
[coords_kT1, x] = BrownianDynamics(1, D, dt, lambda, num_loops, true, "kT = 1", 3, false);
[coords_kT0p5, x] = BrownianDynamics(0.5, D, dt, lambda, num_loops, true, "kT = 0.5", 4, false);


%Run BD at kT = 1 for 1000 trials and count number with 0.2 of (1,0)
count = 0;
for j = 1:1000
    [out, isWithin02] = BrownianDynamics(1, D, dt, lambda, num_loops, false, "", 0, true);
    count = count + isWithin02;
end

%Prints number of runs within 0.2 to console.
fprintf("%d of 1000 BD trajectories with kT=1 were within 0.2 of (1,0) basin minima\n", count)














