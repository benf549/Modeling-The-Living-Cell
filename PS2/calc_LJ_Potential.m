%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/14/21, coded on MATLAB _R2020b_ 
% 
% Calculates the Lennard Jones Potential for a box of particles using the
% minimum image convention for periodic boundary conditions. 

% Output is a double corresponding to the potential energy from the interaction between
% the two vectors ri and rj. 
% 
% Input is L - the length of the box, epsilon -
% the reduced energy unit, sigma - the reduced unit of length and the two
% R^3 vectors ri and rj you want to find the potential between. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function v = calc_LJ_Potential(L, epsilon, sigma, ri, rj)

    %Calculate the difference between each vector component
    Deltx = ri(1) - rj(1);
    Delty = ri(2) - rj(2);
    Deltz = ri(3) - rj(3);
    
    %Since the delta cannot exceed L, dividing the delta by L and
    %casting to int gives -1/0/1 as possible values. Then multiply by L to
    %get the component's nearest image distance on subtraction from delta.
    Deltx = Deltx - L*double(int8(Deltx/L));
    Delty = Delty - L*double(int8(Delty/L));
    Deltz = Deltz - L*double(int8(Deltz/L));
    
    %Calculate magnitude of the distance between vectors
    dist = sqrt(Deltx^2 + Delty^2 + Deltz^2);

    %Apply Lennard-Jones Potential formula and return calculated potential.
    v = 4*epsilon*( (sigma/dist)^12 - (sigma/dist)^6 );
end
