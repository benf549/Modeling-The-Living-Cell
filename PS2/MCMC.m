%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/14/21, coded on MATLAB _R2020b_ 
% 
% Performs a Markov Chain Monte Carlo simulation on a matrix of R^3 vectors.
% Applies the Monte Carlo algorithm and generates a vector of Lennard-Jones
% potential energies representing the energies of the system for each step
% in random walk.
%
% Output is a vector of energies corresponding to the energy calculated for
% each step in the random walk. 
% 
% Input is the number of steps in random walk, the algorithm parameters
% delta, beta, L, epsilon, and sigma as explained in Frenkel and Smit ch 3,
% and the matrix of R^3 vectors on which the random walk will be performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [energy_array] = MCMC(num_cycles, delta, beta, L, epsilon, sigma, data_matrix)
    
    temp_out = zeros(1, num_cycles);
    
    for i = 1:num_cycles
    
        %select a random particle
        selection = randi([1,30]);
        old = [data_matrix(selection, 1) data_matrix(selection, 2) data_matrix(selection, 3)];
    
        %calculate energy of old state. Calculate it if first call,
        %otherwise use memoized value to avoid repeat calculation
        if i > 1 
            old_energy = temp_out(i-1);
        else
            old_energy = calc_LJ_Potential(L, epsilon, sigma, data_matrix);
        end
        
        %Randomly independently adjust x, y, z coordinates of old
        data_matrix(selection, 1) = old(1) + delta*(rand()-0.5);
        data_matrix(selection, 2) = old(2) + delta*(rand()-0.5);
        data_matrix(selection, 3) = old(3) + delta*(rand()-0.5);
    
        %Calculate energy of new state
        new_energy = calc_LJ_Potential(L, epsilon, sigma, data_matrix);
        
        %Define acceptance criterion with Boltzmann factor.
        acc = min([1, exp(-beta * (new_energy - old_energy))]);
        
        %If a randomly generated number is less than the acceptance criterion, we accept the old to new transition
        if rand() < acc
            %Append the calculated energy of the new state to energy array for calculation of average
            temp_out(i) = new_energy;

        else
            %If we do not accept the transition, we keep the old position and energy and append to energy array for average calculation
            data_matrix(selection, 1) = old(1);
            data_matrix(selection, 2) = old(2);
            data_matrix(selection, 3) = old(3);
            temp_out(i) = old_energy;
        end
    end
    %return array of calculated energies.
    energy_array = temp_out;
end

