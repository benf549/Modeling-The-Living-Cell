function coords = BrownianDynamics(kT, D, dt, lambda, num_loops, doPlot, description, fignum)
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
    end
    
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