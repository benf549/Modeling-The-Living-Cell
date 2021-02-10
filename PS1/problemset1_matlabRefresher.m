function[]=prolemset1()

%%%%%%%%%% problem 1: help/functions in matlab%%%%%%%%%
%%%%%COMMANDS: help

%find out what randn is in matlab. what inputs does it require?
%what outputs does it produce?



%%%%%%%%%%%%%%%problem 2: load and plot %%%%%%%%%%%%%%
%%%%%COMMANDS: load, size, plot

  %load in a file, data_prob2.dat. It is xy data, column 1 is x, column 2 is y


%what is the size of this?

%how many rows?

%how many columns?

%rename column1 as xvar

%rename column2 as yvar

%multiply column 2 by 2, name as y2.

%plot xvar vs the new variable, y2.




%%%%%%%%%%%%%%%problem 3: matrices%%%%%%%%%
%%%%COMMANDS: magic, sum, size

%create a magic matrix magic(5)
% a = magic(5);
% size(a)

%what are the dimensions of the matrix?

%what does each row sum to?





%%%%% problem 4 , debug a function %%%%%%%% 


%open the script fix_me_script.m

%what inputs does it take?

%try and run the script on the command line, what happens?

%debug the script so that it works properly



%%%%%%%%problem 5, function definition %%%%%%%%%
%%%%COMMANDS: plot, xlabel, ylabel

%create a new script with the first line: 
%function[]=plotgauss(mu, sigma, dx)

%define xmin as mu-sigma*4

%define xmax=mu+sigma*4

%create an array of xvalues from xmin to xmax, at increments of dx

%define the gaussian f=1/sqrt(2*pi*sigma^2)*exp(-1/(2*sigma^2)*(x-mu)^2)

%plot x vs f.

%label xaxis x

%label yaxis p(x)

%change the fontsize to 40

%change the linewidth of the curve to 3

%save the figure using the saveas command to a file named 'mygauss.eps'




%%%%%%%%%%%%problem 6, surface & contour plots %%%%%%%%%%%%
%%%%COMMANDS: meshgrid, surf, contour

%create a meshgrid of x values from -5 to +5, and y values from -10 to +10

%define a function on this grid Z=X^2*cos(X)-Y^2;

%create a surface plot of Z vs X and Y

% hold this plot

%add a contour plot of Z vs X and Y



%%%%%%%%problem 7, calculate means and errors, plot%%%%%%%
%%%%%COMMANDS: load, size, mean, std, errorplot 

%load datapts3.dat file

%what is the size of the array?

%calculate the mean of the data pts

%calculate the standard deviation of the data pts

%calculate the standard error of the mean (SEM) for the first 10 data points.
%SEM is the sample standard deviation divided by the square root of the
%sample size.

%calculate the mean of the first 10 data points

%create an array that stores the number of data points you just sampled. (i.e.
%npts(1)=10;

%create an array to store the SEM you just calculated (i.e. semvec(1)=..)
%and another for the mean mymeanvec(1)=..)


%calculate the standard error of the mean for the first 100 data points.

%store these as the next entries in your datapoints array and SEM array

%repeat for the first 200 and 1000 datapoints

%make an errorplot of the mean as a function of number of datapoints, with
%the SEM used as the size of the bars



   %%%%%%%%%%%%problem 8. 
 %  Download Visual Molecular Dynamics (VMD) for free from:  www.ks.uiuc.edu)
   %%%load in a file (new molecule). 
 %%load in the clath7_col.psf. VMD should recognize it as a PSF file.
   %%Then load in the trajectory: clathrin713frames.xyz
%%   Under graphics->Representations, change the coloring method to a different ColorID (not gray). Change the drawing method to licorice, with a bond radius of 1.7. Watch the movie, and Render (under File) an image of one of the frames.



   %%%%%%%%%%problem 9, physical biology of the cell, problem 2.6%%%
%%%%%%estimating organelle sizes%%%%%
