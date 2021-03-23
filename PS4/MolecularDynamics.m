

function [energies, positions] = MolecularDynamics(positions, sigma_star, epsilon_star, T_star, mass_star, L, dt, maxtime)    
    num_loops = maxtime/dt + 1; %simulate out to 20 time units

    %Initialize Velocities and Zero the Momentum
    velocities = randn(length(positions), 3) .* sqrt(T_star);
    adjustment = sum(velocities)/length(velocities);
    velocities = velocities - adjustment;

    %Initialize vectors to hold potential and kinetic energies and temperature calculations
    potentials = zeros(length(num_loops), 1);
    kinetics = zeros(length(num_loops), 1);
    temperatures = zeros(length(num_loops), 1);

    %Initialize Forces and Potential at t=0
    [potential, forces] = calculateForces(positions, sigma_star, epsilon_star, L);
    potentials(1) = potential;
    
    %Calcualte kinetic energy at t=0
    kinetic = sum(sum( (mass_star/2) .* (velocities .^ 2) ));    
    kinetics(1) = kinetic;
    
    %Calculate temperature at t=0
    deg_freedom = 3*length(positions) - 3;
    temperature = kinetic/deg_freedom;
    temperatures(1) = temperature;

    for k = 2:num_loops
        %Update positions
        positions = positions + (dt.*velocities) + (( (dt^2) / (2*mass_star) ) .* forces);
        %Partial Velocity Update
        velocities = velocities + ((dt/(2*mass_star)) .* forces);
        %Update Forces
        [potential, forces] = calculateForces(positions, sigma_star, epsilon_star, L);
        potentials(k) = potential;
        %Complete Velocity Update
        velocities = velocities + ((dt/(2*mass_star)) .* forces);
        %Calculate new KE and T
        kinetic = sum(sum( (mass_star/2) .* (velocities .^ 2) ));  
        kinetics(k) = kinetic;
        temperatures(k) = kinetic/deg_freedom;
    end
    
    energies = [potentials; kinetics; temperatures];
end


function [potential, forces] = calculateForces(positions, sigma, epsilon, L)
    potential = 0;
    forces = zeros(length(positions), 3);
    for i = 1:length(positions)
        for j = i+1:length(positions)
            %absolute positions of i and j
            ri = positions(i, :);
            rj = positions(j, :);

            %Calculate Distance Between Pair
            DeltX = ri(1)-rj(1);
            DeltY = ri(2)-rj(2);
            DeltZ = ri(3)-rj(3);

            %Enforce PBC
            DeltX = DeltX - L*round(DeltX/L);
            DeltY = DeltY - L*round(DeltY/L);
            DeltZ = DeltZ - L*round(DeltZ/L);
            Deltas = [DeltX, DeltY, DeltZ];

            %Calculate Pairwise Potentials and Forces
            rij_squared = (DeltX)^2 + (DeltY)^2 + (DeltZ)^2;
            [V_ij, F_ij] = calculateLJ(rij_squared, Deltas, sigma, epsilon);
            potential = potential + V_ij;
            forces(i, :) =  forces(i, :) + F_ij;
            forces(j, :) =  forces(j, :) - F_ij;
        end
    end
end

function [PE, F] = calculateLJ(r_ij_squared, deltas, sigma, epsilon)
    s12 = sigma^12;
    s6 = sigma^6;
    r12 = r_ij_squared^6;
    r6 = r_ij_squared^3;
    
    PE = 4*epsilon*((s12/r12) - (s6/r6));
    F = deltas.*( ((48*epsilon)/r_ij_squared) * ( (s12/r12) - (0.5*s6/r6) ));
end







