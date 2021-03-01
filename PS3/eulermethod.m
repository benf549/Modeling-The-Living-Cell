%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/28/21, coded on MATLAB _R2020b_ 
% 
% Implements Euler's numerical method for ODE approximation on the function
% y' = y which has an analytical solution of y = e^x. Works by calculating
% value of the function a small step away using a first order taylor
% series. 
%
% Input is the previous time value, the previous y value, and the step size h. 
%
% Output is a vector containing the newly generated y value and the newly generated t value.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ynew, tnew] = eulermethod(told, yold, h)
    %first order taylor series of the form ynew = yold + dt*dy/dt
    ynew = yold + h*yold;
    tnew = told + h;
end