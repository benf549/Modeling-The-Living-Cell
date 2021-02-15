function [energy_array] = MCMC(num_cycles, delta, beta, L, epsilon, sigma, data_matrix)
    temp_out = zeros(1, num_cycles);
    
    for i = 1:num_cycles
    
        %select a random particle
        selection = randi([1,30]);
        old = [data_matrix(selection, 1) data_matrix(selection, 2) data_matrix(selection, 3)];
    
        %calculate energy of old state
        old_energy = 0;
        for j = 1:30
            if j ~= selection
                old_energy = old_energy + calc_LJ_Potential(L, epsilon, sigma, old, [data_matrix(j,1) data_matrix(j,2) data_matrix(j,3)]);
            end
        end
        disp(old_energy)
        %Randomly independently adjust x, y, z coordinates of old
        new = zeros(1,3);
        new(1) = old(1) + delta*(rand()-0.5);
        new(2) = old(2) + delta*(rand()-0.5);
        new(3) = old(3) + delta*(rand()-0.5);
    
        %Calculate energy of new state
        new_energy = 0;
        for j = 1:30
            if j ~= selection
                new_energy = new_energy + calc_LJ_Potential(L, epsilon, sigma, new, [data_matrix(j,1) data_matrix(j,2) data_matrix(j,3)]);
            end
        end

        %Define acceptance criterion with Boltzmann factor.
        acc = min([1, exp(-beta * (new_energy - old_energy))]);
    
        %If a randomly generated number is less than the acceptance criterion, we accept the old to new transition
        if rand() < acc
            %Update old in data matrix with new value
            data_matrix(selection, 1) = new(1); 
            data_matrix(selection, 2) = new(2); 
            data_matrix(selection, 3) = new(3); 

            %Append the calculated energy of the new state to energy array for calculation of average
            temp_out(i) = new_energy;
        else
            %If we do not accept the transition, we keep the old energy and append to energy array for average calculation
            temp_out(i) = old_energy;
        end
    end
    energy_array = temp_out;
end

