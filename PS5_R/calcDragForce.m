%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 4/05/21, coded on MATLAB _R2020b_ 
% 
% Creates a matrix of drag forces to be added to total forces on beads in
% Langevin Dynamics simulation. The drag force opposes the velocity in each
% direction. 
% 
% Input is the matrix of velocities of each bead (row) in the polymer, the
% mass of each bead, and the drag coefficient
%
% Output is a matrix representing the force on each bead (row) as a result
% of the drag force.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [force] = calcDragForce(veloc, mass, drag_coeff)
    %calculate drag force as -m*drag_coeff*velocity
    force = -mass*drag_coeff.*veloc;
end