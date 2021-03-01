%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/27/21, coded on MATLAB _R2020b_ 
% 
%Problem Set 3 - Modeling The Living Cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1
maxt = 10;

%
%Plot from 1 to 10 when h = 1
h1 = 1;
numsteps = maxt/h1 + 1;

%Initialize arrays for indexing
ys = zeros(1, numsteps);
ts = zeros(1, numsteps);

%Initial Condition y(0) = 1
ys(1) = 1;

for i = 2:(numsteps)
    [ynew, tnew] = eulermethod(ts(i-1), ys(i-1), h1);
    ts(i) = tnew;
    ys(i) = ynew;
end

%
%Plot from 1 to 10 when h = 0.5
h2 = 0.5;
numsteps2 = maxt/h2 + 1;
ys2 = zeros(1, numsteps2);
ts2 = zeros(1, numsteps2);
ys2(1) = 1;
for i = 2:(numsteps2)
    [ynew, tnew] = eulermethod(ts2(i-1), ys2(i-1), h2);
    ts2(i) = tnew;
    ys2(i) = ynew;
end

%
%Plot from 1 to 10 when h = 0.1
h3 = 0.1;
numsteps3 = maxt/h3 + 1;
ys3 = zeros(1, numsteps3);
ts3 = zeros(1, numsteps3);
ys3(1) = 1;
for i = 2:(numsteps3)
    [ynew, tnew] = eulermethod(ts3(i-1), ys3(i-1), h3);
    ts3(i) = tnew;
    ys3(i) = ynew;
end

figure(1)
plot(ts, ys,'Color', 'red')
hold on
plot(ts2, ys2,'Color', 'blue')
plot(ts3, ys3,'Color', 'green')
plot(ts3, exp(ts3),'Color', 'black')
legend("h=1.0","h=0.5","h=0.1","e^x")
title("Plot of Euler's Method Approximations for y' = y ODE")
xlabel("t")
ylabel("y")
hold off

figure(2)
err1 = exp(ts) - ys;
err2 = exp(ts2) - ys2;
err3 = exp(ts3) - ys3;
plot(ts, err1,'Color', 'red')
hold on
plot(ts2, err2,'Color', 'blue')
plot(ts3, err3,'Color', 'green')
legend("h=1.0","h=0.5","h=0.1")
title("Plot of Error as Known Value - Approximation")
xlabel("t")
ylabel("Error (e^t - appoximation(t))")
hold off

%% Problem 2
%Set up initial conditions matrix
y1_init = 55;
y2_init = 110;
y0 = [y1_init, y2_init];

%Initialize rate constants
k1 = 12; %tu-1
k2 = 0.04; %mi^2/tu-1
k4 = 8; %tu-1

%Arbitrarily generate data from t=0 to t=5
tspan = [0 5];

%Solve differential equations implemented in predatorprey.m
[t, y] = ode45(@(t,y) predatorprey(t,y, k1, k2, k4), tspan, y0);

%Plot solutions for predator and prey as functions of time. 
figure(3)
plot(t, y(:, 1))
hold on
plot(t, y(:, 2))
title("Lotka-Volterra System Solved With ODE45")
legend("Prey (y_1)", "Predator (y_2)")
ylabel("y_i (species i per mi^2)")
xlabel("time (tu)")
hold off

%% Problem 3
rng('shuffle')
y1_init = 55;
y2_init = 110;

k1 = 12; %tu-1
k2 = 0.04; %mi^2/tu-1
k4 = 8; %tu-1
A = 5; %mi^2

%copy numbers
num_y1 = A * y1_init; 
num_y2 = A * y2_init;
y0 = [num_y1, num_y2];

%stochastic rate constants
c1 = k2/A; %tu-1
c2 = k1; %tu-1
c3 = k4; %tu-1
srcs = [c1, c2, c3];

%Update matrix based on reaction equations. Row is reaction, col is y1 and y2
update_matrix = [-1 1; 1 0; 0 -1];

num_loops = 100000;

output_mtx1 = gillespie(y0, srcs, update_matrix, num_loops);

figure(4)
plot(output_mtx1(1, :), output_mtx1(2, :), 'r')
hold on
plot(output_mtx1(1, :), output_mtx1(3, :), 'b')


output_mtx2 = gillespie(y0, srcs, update_matrix, num_loops);
plot(output_mtx2(1, :), output_mtx2(2, :), 'r--')
plot(output_mtx2(1, :), output_mtx2(3, :), 'b--')
legend("prey (y1) trial 1", "predator (y2) trial 1", "prey (y1) trial 2", "predator (y2) trial 2", "location", "best")
ylabel("Copy Number")
xlabel("time (time units)")
title("Lotka-Volterra System Copy Numbers Solved With Gillespie (two trials)")
hold off

figure(5)
plot(t, y(:, 1).*A, 'r')
hold on
plot(t, y(:, 2).*A, 'b')
plot(output_mtx1(1, :), output_mtx1(2, :), 'r--')
plot(output_mtx1(1, :), output_mtx1(3, :), 'b--')

ymax = max([y(:, 1).*A'; y(:, 2).*A'; output_mtx1(2, :)'; output_mtx1(3, :)']);
xmax = max(output_mtx1(1, :));
axis([0, xmax, 0,ymax]) 
title("Lotka-Volterra System Solutions: Comparing Gillespie and ODE45")
ylabel("Copy Number")
xlabel("time (time units)")
legend("prey (y1) ODE45", "predator (y2) ODE45", "prey (y1) Gillespie", "predator (y2) Gillespie", "location", "best")
hold off

