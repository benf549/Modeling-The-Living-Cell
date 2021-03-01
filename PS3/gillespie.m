%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/28/21, coded on MATLAB _R2020b_ 
% 
% Implements the Gillespie algorithm as generally as possible for future
% applications. This algorithm stochastically solves ordinary differential
% equations such as reaction kinetics equations. 
% To reapply this solution to another system just need to
% adjust the reaction propensities below.
%
% Input is the initial conditions of each speices in a vector and the
% stochastic rate constants in a vector (for both of which, index of
% vector corresponds to number of species), an update matrix that is based
% on the net result of each reaction its rows are the mu^th reaction and
% the columns are the species number in ascending order. Finally, takes in
% the number of loops to apply the algorithm. 
%
% Output is a matrix of doubles. The first row is the time steps taken and
% the subsequent rows are the copy numbers for the i^th species.% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function out = gillespie(y0, stoch_rxn_consts, update_matrix, num_loops)
    out = zeros(length(y0) + 1, num_loops + 1);
    for i = 1:length(y0)
        out(i+1, 1) = y0(i);
    end
    
    for i = 1:num_loops
        %reaction propensities
        a1 = out(2, i) * out(3, i) * stoch_rxn_consts(1);
        a2 = out(2, i) * stoch_rxn_consts(2);
        a3 = out(3, i) * stoch_rxn_consts(3);
        as = [a1, a2, a3];

        %total propensity and random time step.
        atot = sum(as);
        tau = (-1/atot)*log(rand());
        
        %Calulate and sample from cdf
        cdf = cumsum(as);
        breakval = rand()*atot;
        for j = 1:length(as)
           if breakval < cdf(j)
                break
           end
        end
        
        %Update the previous values based on the selected reaction and
        %increment time by time step.
        out(1, i+1) = out(1, i) + tau;
        out(2, i+1) = out(2, i) + update_matrix(j, 1);
        out(3, i+1) = out(3, i) + update_matrix(j, 2);
    end
end
