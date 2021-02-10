%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjmain Fry (bfry2)
% 01/28/21, coded on MATLAB R2020b
% 
% INPUTS: mu, sigma, dx
% OUTPUTS: none
% plots the gaussian function given a mu, sigma, and dx. Saves the plot to 'mygauss.eps'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[]=plotgauss(mu, sigma, dx)
    xmin = mu - sigma*4;
    xmax = mu + sigma*4;

    x = xmin:dx:xmax;
    f = 1/sqrt(2*pi*sigma^2)*exp(-1/(2*sigma^2)*(x-mu).^2);

    a = figure();
    plot(x, f, 'LineWidth', 3)
    xlabel('x')
    ylabel('p(x)')
    set(gca,'FontSize',40)
    saveas(a, 'mygauss.eps', "psc2")
end