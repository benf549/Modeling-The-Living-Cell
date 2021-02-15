%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/14/21, coded on MATLAB _R2020b_ 
% 
% Calculates the Lennard Jones Potential for a box of particles using the
% minimum image convention for periodic boundary conditions. 

% Output is a double corresponding to the potential energy of the unique pairwise
% interactions of the particles
% 
% Input is L - the length of the box, epsilon -
% the reduced energy unit, sigma - the reduced unit of length and the
% matrix of values for which to sum over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function V = calc_LJ_Potential(L, epsilon, sigma, data_matrix)
    v = 0;
    for i = 1:length(data_matrix)
        for j = i+1:length(data_matrix)
            
            %Calculate the difference between each vector component
            Deltx = data_matrix(i,1) - data_matrix(j,1);
            Delty = data_matrix(i,2) - data_matrix(j,2);
            Deltz = data_matrix(i,3) - data_matrix(j,3);

            %Since the delta cannot exceed L, dividing the delta by L and
            %casting to int gives -1/0/1 as possible values. Then multiply by L to
            %get the component's nearest image distance on subtraction from delta.
            Deltx = Deltx - L*double(int8(Deltx/L));
            Delty = Delty - L*double(int8(Delty/L));
            Deltz = Deltz - L*double(int8(Deltz/L));

            %Calculate magnitude of the distance between vectors
            dist = sqrt(Deltx^2 + Delty^2 + Deltz^2);

            %Apply Lennard-Jones Potential formula and return calculated potential.
            v = v + 4*epsilon*( (sigma/dist)^12 - (sigma/dist)^6 );
        end
    end
    V = v;
end
