%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 4/05/21, coded on MATLAB _R2020b_ 
% 
% Creates a matrix of configurational forces acting on the beads of the
% polymer being simulated in our Langevin Dynamics simulation and
% calculates the total potential energy of a given configuaration. For
% adjacent beads, the force is due to a harmonic oscillator potential
% representing bonding while non-adjacent beads have the Lennard Jones
% potential applied with parameters varying by pairwise polarity interations.
% 
% Input is the position of each bead (row) and the polarity of the bead as
% a vector of labels.
%
% Output is a matrix representing the force on each bead (row) as a result
% of the bead configuration and a number representing the total potential of 
% the configuation.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [forces, enpot] = calcConfigForce(positions, labels)
    %LJ Constants
    sigma = 1;
    epsilon_HH = 1;
    epsilon_HP = 2/3;
    epsilon_PP = 2/3;
    
    %Bond Harmonic Oscillator Constants
    k = 20;
    l_e = 1;

    %initialize force mtx
    numres = length(positions);
    forces = zeros(numres, 3);
    
    %initialize potential counter
    enpot = 0;
    
    %loop over residues
    for i=1:numres
        for j=i+1:numres
            %Test if adjacent, if so apply bond potential otherwise apply
            %LJ potential by type of residue each bead is
            if abs(j-i) == 1
               [forces, pe] = BondPotential(i, j, positions, forces, k, l_e);
            elseif char(labels(i)) == 'H' && char(labels(j)) == 'H'
               [forces, pe] = LennardJonesPotential(i, j, positions, forces, sigma, epsilon_HH);
            elseif char(labels(i)) == 'P' && char(labels(j)) == 'P'
               [forces, pe] = LennardJonesPotential(i, j, positions, forces, sigma, epsilon_PP);
            else
               [forces, pe] = LennardJonesPotential(i, j, positions, forces, sigma, epsilon_HP);

            end
            enpot = enpot + pe;
        end
    end


end

function [forces, pe] = BondPotential(i, j, positions, forces, k, l_e)
    %calculate axis-wise distances
    xij = positions(i, 1) - positions(j, 1);
    yij = positions(i, 2) - positions(j, 2);
    zij = positions(i, 3) - positions(j, 3);
    
    %calculate distance between beads.
    rij = sqrt(xij^2 + yij^2 + zij^2);

    %total forces from F = -dV/dr
    f = -k*(rij - l_e)/rij;
    
    %update force on beads scaled by axis
    forces(i, :) = forces(i, :) + f*[xij yij zij];
    forces(j, :) = forces(j, :) - f*[xij yij zij];
    
    %calculate potential
    pe = 0.5*k*((rij - l_e)^2);
end

function [forces, pe] = LennardJonesPotential(i, j, positions, forces, sigma, eps)
    %calculate axis-wise distances
    xij = positions(i,1) - positions(j,1);
    yij = positions(i,2) - positions(j,2);
    zij = positions(i,3) - positions(j,3);    

    % Squared distance in 3 dimensions
    r2 = xij^2 + yij^2 + zij^2;
    sigma2 = sigma^2;
    r2i = sigma2/r2;
    r6i = r2i^3;
    
    %Force determined from F = -dV/dr
    f = 48 * eps * r2i * (r6i^2 - 0.5*r6i);

    %update force on beads scaled by axis
    forces(i,:) = forces(i,:) + f*[xij yij zij];
    forces(j,:) = forces(j,:) - f*[xij yij zij];
    
    %calculate pair potential
    pe = 4*eps*(r6i^2 - r6i);
end