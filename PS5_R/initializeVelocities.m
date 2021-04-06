
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 4/05/21, coded on MATLAB _R2020b_ 
% 
% Initializes velocities for Langevin Dynamics simulation. Ensures momentum
% of system is zero.
% 
% Input is the number of residues to initialize 3-dimensional velocities
% for, the reduced boltzmann constant, the temperature of the system, and
% the mass of each bead.
%
% Output is the matrix of initialized velocities where each row represents
% a bead and each column is the X, Y, and Z component respectively.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [velocities] = initializeVelocities(numres, k, T, m)
    %Initialize random velocities as normal with sd sqrt(kt/m)
    velocities = randn(numres, 3) * sqrt(k*T/m);
    
    %Calculate average in each dimension
    avg_v = sum(velocities)/numres;
    
    %Zero momentum
    velocities = velocities - avg_v;
end
