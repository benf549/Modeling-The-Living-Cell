%% Problem 1 - Langevin Dynamics
rng(100)
%Parameters
filename = "5bead.inp";
nstep = 10000;
temp = 2;
drag_coeff = 0.05; %s^-1
mass = 1;
k = 1;
time_step = 0.003; %s

%Run LD Simulation with given parameters
[pes, kes] = LangevinDynamics(filename, nstep, temp, drag_coeff, mass, k, time_step);

%Plot output
figure(1)
plot(0:time_step:nstep*time_step, pes)
hold on
plot(0:time_step:nstep*time_step, kes)
legend("Potential Energy", "Kinetic Energy")
xlabel("Time (seconds)")
ylabel("Energy (reduced units)")
title("Plot of Energy as a Function of Time Over 10,000 Time Steps");
hold off

%Calculate averages for final 1000 steps
avg_pe_1000 = sum(pes(nstep+1-1000:nstep+1))/1000;
avg_ke_1000 = sum(kes(nstep+1-1000:nstep+1))/1000;
fprintf("The final 1000 time steps had an average\nPE = %f\nKE = %f\n", ...
    avg_pe_1000, avg_ke_1000);



