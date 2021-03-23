function PlotOutput(output, xstep, xmax, plotnum, titlemod) 
    figure(plotnum);
    plot(0:xstep:xmax, output(1, :));
    hold on;
    plot(0:xstep:xmax, output(2, :));
    plot(0:xstep:xmax, output(1, :) + output(2, :));
    legend("Potential", "Kinetic",  "Total", 'Location', 'best');
    title(strcat("MD Energy as a Function of Time: ", titlemod))
    ylabel("Energy (Reduced Energy Units)")
    xlabel("Time (Reduced Time Units)")
    hold off;

    figure(plotnum+1);
    histogram(output(3, :))
    ylabel("Count")
    xlabel("Temperature (Reduced Temperature Units)")
    title(strcat("MD Temperature Measurement Distribution: ", titlemod))
end