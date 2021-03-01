%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/28/21, coded on MATLAB _R2020b_ 
% 
% A function that expresses the Lotka-Volterra Predator prey system for
% solution by MATLAB's built in ODE45 solver. Function is formatted to
% work as an input for the 'ODEFUN' parameter required in diffeq solver.
%
% Input is a time, a vector containing the values of y1 and y2, and the
% rate constants k1, k2, and k4.
%
% Output is a vector of rates of change for y1 and y2 respectively. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = predatorprey(t, y, k1, k2, k4)
    dydt = zeros(2,1);
    %rate of change for species 1 given in problem
    dydt(1) = k1*y(1) - k2*y(1)*y(2);
    %rate of change for species 2
    dydt(2) = k2*y(1)*y(2) - k4*y(2);
end
