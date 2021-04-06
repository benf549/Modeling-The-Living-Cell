%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 3/25/21, coded on MATLAB _R2020b_ 
% 
% Performs a Brownian Dynamics simulation for a random walk in a two-well
% potential for a single particle. 
%
% Output is a matrix of coordinates of the random walk. The first column is
% the x values and the second column is the y values. The second output
% term is a boolean that indicates whether the walk came within 0.2 of the
% second well and can only be true if breakonwithin02 is set to true.
% 
% Input is the temperature in kT, the diffusion constant, the time step,
% the energy function constant lambda, the number of loops to perform, a
% boolean to decide whether to plot the results or not, a description that
% describes the conditions for a plot, the figure number of the plot, and a
% boolean for whether to break the simulation upon coming within 0.2 of the
% second energy well. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coords, within02] = BrownianDynamics(kT, D, dt, lambda, num_loops, doPlot, description, fignum, breakonwithin02)
    coords = zeros(num_loops, 2);
    coords(1, 1) = -1;
    coords(1, 2) = 0;
    

    for i = 2:num_loops
        %update x with euler method
        coords(i, 1) = coords(i-1, 1) - (( (dt*D)/(kT) ) * ( 20*coords(i-1,1)* ...
            (coords(i-1,1)^2 - 1))) + (sqrt(2*D*dt) * randn(1));
        %update y with euler method
        coords(i, 2) = coords(i-1, 2) - (( (dt*D)/(kT) ) * ( 20*coords(i-1,2)*lambda)) + ...
            (sqrt(2*D*dt) * randn(1));
        
        %detect if path comes within 0.2 of second basin minima at (1, 0)
        within02 = false;
        if (breakonwithin02)
            dist_from_basin10 = sqrt( (coords(num_loops, 1) - 1)^2 + (coords(num_loops, 2) - 0)^2 );
            if (dist_from_basin10 < 0.2) 
                within02 = true;
                break %break out of loop if comes within 0.2 if desired
            end
        end
    end
    
    %Plot if desired
    if (doPlot)
        PlotWalk(coords, lambda, description, fignum);
    end
    
end



function PlotWalk(coords, lambda, description, fignum) 
    %create contour plot of potential
    xs = linspace(-1.5, 1.5);
    ys = linspace(-2, 2);
    [X, Y] = meshgrid(xs, ys);
    Z = (5*( X.^2 - 1).^2) + (10 .* (Y.^2) .* lambda);
    
    figure(8+fignum)
    plot(coords(:, 1), coords(:, 2))
    hold on
    contour(X, Y, Z)
    axis([-1.5, 1.5, -2, 2])
    text(coords(1,1), coords(1, 2), "\leftarrow START", "Color", "red")
    text(coords(50001,1), coords(50001, 2), "\leftarrow END", "Color", "red")
    title(strcat("Plot of Brownian Motion In Bistable Potential: ", description))
    xlabel("X")
    ylabel("Y")
    hold off
end