%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjmain Fry (bfry2)
% 01/28/21, coded on MATLAB R2020b
% 
% Solutions to the various problems in PS1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1 - help functions

% Code: 
help randn

% Written answer: randn is a function in matlab that generates pseudorandom numbers sampled from the standard normal distribution. it can take in 1 input in which case it will output an NxN matrix of these random numbers.
% if passed two numbers separated by a comma, the randn(m, n) will output an MxN matrix of these random numbers.


%% Problem 2 - load and plot
load('data_prob2.dat');
size(data_prob2)
% The size of this data is 50 x 2
% This means 50 rows and 2 columns.

xvar = data_prob2(:, 1); %rename column1 as xvar
yvar = data_prob2(:, 2); %rename column2 as yvar
y2 = yvar * 2; %multiply column2 by 2, name as y2

plot(xvar, y2) %plot xvar vs the new variable, y2


%% Problem 3 - matrices
a = magic(5);
size(a) %the matrix is of size 5x5 (5 rows, 5 columns)

disp("sum of rows: ")
sum(a, 2)
%from the above code, it appears that each row sums to 65.


%% Problem 4 - debug a function
fix_me_script(pi, 0.2, 1)
% the function takes in two general numeric inputs corresponding to the frequency and damping constant for a damped oscillator. The third input is a positive integer corresponding to the figure number for the generated plot.
% when run before fixing, an error is thrown because the f1 and f2 vectors are the wrong shape for matrix multiplication. In this scenario, we probably want to perform element-wise multiplication so we change the operator to .* to do this.


%% Problem 5 - function definition
plotgauss(0, 1, 0.1)


%% Problem 6 - surface & contour plots
x  = -5:5;
y = -10:10;

[X, Y] = meshgrid(x, y);
Z = X.^2 .* cos(X) - Y.^2;

surf(X, Y, Z)
hold ON
contour(X, Y, Z)
hold OFF

%% Problem 7 - means, errors, and plot
load datapts3.dat

size(datapts3)
% The size of the array is 1000x1 (1000 rows, 1 column)

totmean = mean(datapts3)
% The mean of the data pts is 0.1911

datastd = std(datapts3)
% The standard deviation of the data pts is 5.7164

f10 = datapts3(1:10);
f10_sem = std(f10)/sqrt(10)
% The SEM of the first 10 datapoints is 1.7001

f10_mean = mean(f10)
% The mean of the first 10 datapoints is -1.5777

npts(1) = 10;
semvec(1) = f10_sem;
meanvec(1) = f10_mean;
% Created vectors to store the npts, sem, and mean respectively.

count = 1;
for i = [100 200 1000]
    count = count + 1;
    fi = datapts3(1:i);
    fi_sem = std(fi)/sqrt(i);
    fi_mean = mean(fi);
    npts(count) = i;
    semvec(count) = fi_sem;
    meanvec(count) = fi_mean;
end
% Calculate the statistics and add for the first 100, 200, and 1000 datapoints to respective vectors after calculation.
errorbar(npts, meanvec, semvec, "LineStyle", "none", "Marker", ".", "MarkerSize", 15, "MarkerEdgeColor", "red")
grid on

%% Problem 8 - in VMD

%% Problem 9 - see written work.
