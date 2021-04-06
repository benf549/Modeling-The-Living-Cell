%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 4/05/21, coded on MATLAB _R2020b_ 
% 
% Performs Langevin Dynamics simulation for coarse-grained sequence of Polar and
% Hydrophobic beads. Reads data from .inp file containing bead polarity and
% initial coordinates. Implements Velocity Verlet Algorithm modified to
% work with low-friction langevin simulation by applying BBK algorithm.
%
% Outputs two vectors. First the vector of potential energies and second
% the vector of kinetic energies as calculated over the course of the
% simulation.
% 
% Input is the *.inp filename to read in, the number of steps to simulate
% out to, the temperature of the system, the drag coefficient, the mass of
% each bead, the reduced boltzmann constant, and the time step.
%
% Inputs are given in reduced units.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [potentials, kinetics] = LangevinDynamics(filename, nstep, temp, drag_coeff, mass, k, time_step)

    %Initialize Data
    data = importdata(filename, " ", 1);
    coords = data.data;
    labels = data.textdata(2:length(coords)+1);

    %Initialize Velocities
    numres = length(coords);
    velocities = initializeVelocities(numres, k, temp, mass);

    %Initialize Forces
    [config_forces, pe] = calcConfigForce(coords, labels);
    kick_forces = calcRandomKicks(numres, k, temp, mass, drag_coeff, time_step);

    %Initialize output vectors
    potentials = zeros(nstep + 1, 1);
    kinetics = zeros(nstep + 1, 1);
    potentials(1) = pe;
    kinetics(1) = mass*sum(sum(velocities.^2))/2;
    
    for i = 1:nstep
        %Calculate Drag Forces and update forces array
        drag_forces = calcDragForce(velocities, mass, drag_coeff);
        forces = config_forces + kick_forces + drag_forces;

        %Update Positions
        coords = coords + time_step*velocities + 0.5*(time_step^2)*forces/mass;

        %Partial Velocity Update
        prefix = 1/(1 + time_step*drag_coeff/2);
        velocities = prefix * (velocities + time_step*forces/(2*mass));

        %Force Update Based on New Positions
        [config_forces, pe] = calcConfigForce(coords, labels);
        kick_forces = calcRandomKicks(numres, k, temp, mass, drag_coeff, time_step);
        forces = config_forces + kick_forces;

        %Complete Velocity Update
        velocities = velocities + ( prefix * (time_step*forces/(2*mass)) );

        %Calculate New Kinetic Energy
        ke = mass*sum(sum(velocities.^2))/2;

        %Update energy vectors
        potentials(i+1) = pe;
        kinetics(i+1) = ke;
    
    end

end




