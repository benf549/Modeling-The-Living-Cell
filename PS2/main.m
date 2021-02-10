%Benjamin Fry
%bfry2
%Problem Set 2 - Modeling The Living Cell

%% Random Number Generation and Seeding
rng('shuffle')                          

%Adjust random normal distribution to mean of 5 and std.dev 2
r = 5 + 2.*randn(1000, 1);

%Plot histogram of Counts
figure(1)
histogram(r,35)
title('Frequency Plot For 1000 Points Sampled From N(5,4) Distribution');
ylabel('Number of Values Observed');
xlabel('Value Observed');

%Plot histogram of counts as normalized pdf. 
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

%% Sampling From A Discrete Distribution

%Key
aa_key = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', ...
    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];

%Generate a 300 by 1 matrix of random integers between [1, 20]
num_seq = randi(20, 300, 1);
%histogram(sequence)
char_seq = aa_key(num_seq);

%True Frequencies
true_freqs = [8.25, 5.53, 4.06, 5.45, 1.37, 3.93, 6.75, 7.07, 2.27, 5.96, ...
    9.66, 5.84, 2.42, 3.86, 4.70, 6.56, 5.34, 1.08, 2.92, 6.87] ./ 100;
sum(true_freqs)

cumul = cumsum(true_freqs);
cdfplot(true_freqs)