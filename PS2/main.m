%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjamin Fry (bfry2)
% 2/14/21, coded on MATLAB _R2020b_ 
% 
%Problem Set 2 - Modeling The Living Cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Random Number Generation and Seeding

rng('shuffle')

%a. Adjust random normal distribution to mean of 5 and std.dev 2
r = 5 + 2.*randn(1000, 1);

%b. Plot histogram of Counts
figure(1)
histogram(r,35)
title('Frequency Plot For 1000 Points Sampled From N(5,4) Distribution');
ylabel('Number of Values Observed');
xlabel('Value Observed');

%c. Plot histogram of counts as normalized pdf.
figure(2)
histogram(r,35, 'Normalization', 'pdf')
title('PDF of 1000 Points Sampled From N(5,4) Distribution');
ylabel('Probability of Value Observed');
xlabel('Value Observed');
[N, edges] =  histcounts(r,35, 'Normalization', 'pdf');

%check that output of binning sums to 1 when multiplied by edge size.
binsize = edges(2) - edges(1);
disp("Normalized Sum :")
disp(sum(N) * binsize)

%% 2. Sampling From A Discrete Distribution

%Key
aa_key = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', ...
    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];

%Generate a 300 by 1 matrix of random integers on interval [1, 20]
num_seq = randi(20, 300, 1);
figure(1)
histogram(num_seq, 'Normalization', 'pdf')
set(gca,'xtick',[1:20],'xticklabel',aa_key(1:20)');
equal_freq_char_seq = aa_key(num_seq)
title("Residue Frequency in Randomly Generated 300 Residue String Assuming All Have Equal Liklihood")
xlabel("Amino Acid (One Letter Code)")
ylabel("Residue Frequncy (count/string length)")

%True Frequencies
true_freqs = [8.25, 5.53, 4.06, 5.45, 1.37, 3.93, 6.75, 7.07, 2.27, 5.96, ...
    9.66, 5.84, 2.42, 3.86, 4.70, 6.56, 5.34, 1.08, 2.92, 6.87] ./ 100;
sum(true_freqs);

%Calculate Cumulative Sum and plot CDF
cumul = cumsum(true_freqs);
figure(2)
h = stairs(1:20, cumul);
set(gca,'xtick',[1:20],'xticklabel',aa_key(1:20)');
title("Cumulative Distribution Function For Amino Acid Probabilities")
xlabel("Amino Acid (One Letter Code)")
ylabel("Probability")

%Generate 300 random numbers and use CDF to convert to residues
rand300 = rand(300,1);
result = zeros(1,length(rand300));
for i = 1:length(rand300)
    result(i) = getRes(rand300(i), true_freqs);
end
emperical_freq_char_seq = aa_key(result)

figure(3)
histogram(result, 'Normalization', 'pdf');
[model_freqs, edges] = histcounts(result, 'Normalization', 'pdf');
set(gca,'xtick',1:20,'xticklabel',aa_key(1:20)');
title("Residue Frequency in Randomly Generated 300 Residue String Using Expasy Residue Frequencies")
ylabel("Residue Frequncy (count/string length)")
xlabel("Amino Acid (One Letter Code)")


%% 3. Monte Carlo Sampling of Lennard-Jones Particles.

%Load in Data
load('init_crds_boxl_3.5.dat')

%Initialize Variables
L = 3.5; %Length of box
epsilon = 0.25;
sigma = 1;
delta = 0.5;

%Loop through points avoiding self and repeat interactions and sum calculated potential
V = calc_LJ_Potential(L, epsilon, sigma, init_crds_boxl_3_5);
fprintf("Lennard-Jones Potential With Periodic Boundary Conditions: %.4f\n", V)

%Markov Chain Monte Carlo
num_cycles = 100000;

%Calculate with kT = 1
kT1 = 1;
beta1 = 1/kT1;
energies1 = MCMC(num_cycles, delta, beta1, L, epsilon, sigma, init_crds_boxl_3_5);
average_V_1 = sum(energies1)/num_cycles;
fprintf("Calculated Average Potential Energy w/ kT=%.2f as: %.4f\n", kT1, average_V_1)

%Calculate with kT = 0.5 and 0.1
kT2 = 0.5;
beta2 = 1/kT2;
energies2 = MCMC(num_cycles, delta, beta2, L, epsilon, sigma, init_crds_boxl_3_5);
average_V_2 = sum(energies2)/num_cycles;
fprintf("Calculated Average Potential Energy w/ kT=%.2f as: %.4f\n", kT2, average_V_2)

kT3 = 0.1;
beta3 = 1/kT3;
energies3 = MCMC(num_cycles, delta, beta3, L, epsilon, sigma, init_crds_boxl_3_5);
average_V_3 = sum(energies3)/num_cycles;
fprintf("Calculated Average Potential Energy w/ kT=%.2f as: %.4f\n", kT3, average_V_3)
