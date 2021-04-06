%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 4/05/21, coded on MATLAB _R2020b_ 
% 
% Creates a matrix of random kicks to add to force vector in langevin
% dynamics simulation. 
% 
% Input is the number of residues to create 3-dimensional random kicks for,
% the reduced boltzmann constant, the initial temperature of the system, the 
% mass of each bead, the drag coefficient, and the time step.
%
% Output is a matrix representing the force on each bead (row) as a result
% of random kicks
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [force] = calcRandomKicks(numres, k, T, mass, drag_coeff, dt)
    force = randn(numres, 3) .* sqrt(2*k*T*mass*drag_coeff/dt);
end