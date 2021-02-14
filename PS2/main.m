%Benjamin Fry
%bfry2
%Problem Set 2 - Modeling The Living Cell

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

%% Sampling From A Discrete Distribution

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


%% Part 3 - Monte Carlo Sampling of Lennard-Jones Particles.
load('init_crds_boxl_3.5.dat')
L = 3.5; %Length of box
rho_star = 0.7; %Reduced density


%Initialize the Periodic Boundary
neg_x_image = zeros(30,3);
pos_x_image = zeros(30,3);
neg_y_image = zeros(30,3);
pos_y_image = zeros(30,3);
neg_z_image = zeros(30,3);
pos_z_image = zeros(30,3);
for i = 1:30
    r = [init_crds_boxl_3_5(i, 1) init_crds_boxl_3_5(i, 2) init_crds_boxl_3_5(i, 3)];
    
    neg_x_image(i,1) = -L + r(1);
    neg_x_image(i,2) = r(2);
    neg_x_image(i,3) = r(3);
    pos_x_image(i,1) = L + r(1);
    pos_x_image(i,2) = r(2);
    pos_x_image(i,3) = r(3);
    
    neg_y_image(i,1) = r(1);
    neg_y_image(i,2) = -L + r(2);
    neg_y_image(i,3) = r(3);
    pos_y_image(i,1) = r(1);
    pos_y_image(i,2) = L + r(2);
    pos_y_image(i,3) = r(3);
    
    neg_z_image(i,1) = r(1);
    neg_z_image(i,2) = r(2);
    neg_z_image(i,3) = -L + r(3);
    pos_z_image(i,1) = r(1);
    pos_z_image(i,2) = r(2);
    pos_z_image(i,3) = L + r(3);
    
end



%Visualize the input data and initialized periodic boundaries.
%Plot the given data in black
%Plot images with x for negative image and o for positive image and color by reflection axis
plot3(init_crds_boxl_3_5(:, 1), init_crds_boxl_3_5(:, 2), init_crds_boxl_3_5(:, 3), '.', 'Color', 'black')
hold on
plot3(neg_x_image(:, 1), neg_x_image(:, 2), neg_x_image(:, 3), 'x', 'Color', 'blue')
plot3(pos_x_image(:, 1), pos_x_image(:, 2), pos_x_image(:, 3), 'o', 'Color', 'blue')
plot3(neg_y_image(:, 1), neg_y_image(:, 2), neg_y_image(:, 3), 'x', 'Color', 'green')
plot3(pos_y_image(:, 1), pos_y_image(:, 2), pos_y_image(:, 3), 'o', 'Color', 'green')
plot3(neg_z_image(:, 1), neg_z_image(:, 2), neg_z_image(:, 3), 'x', 'Color', 'red')
plot3(pos_z_image(:, 1), pos_z_image(:, 2), pos_z_image(:, 3), 'o', 'Color', 'red')
plot3(0,0,0, 'x', 'Color', 'black')
hold off
xlabel('x')
ylabel('y')
zlabel('z')
axis([-1.75-L 1.75+L -1.75-L 1.75+L -1.75-L 1.75+L])
% axis([-1.75 1.75 -1.75 1.75 -1.75 1.75])
title("Plot of data and nearest images as periodic boundaries")



%Identify the nearest periodic boundary for each point
% V = 0;
% for i = 1:30
%     disp(V)
%     r = [init_crds_boxl_3_5(i, 1), init_crds_boxl_3_5(i, 2), init_crds_boxl_3_5(i, 3)];
%     for j = i:30 %setting to 30 gives correct answer for no p. bound so something is up with the way im working with the boundaries.
%         %Okay so see link in onenote, the problem is a misunderstanding of
%         %minimum image. Read the explanation then rework.
%         if (abs(r(1)) > abs(r(2)) && abs(r(1)) > abs(r(3)))
%             if r(1) > 0
%                 %                 disp('pos x')
%                 v = calc_LJ_Potential(r, j, init_crds_boxl_3_5, pos_x_image);
%             else
%                 %                 disp('neg x')
%                 v = calc_LJ_Potential(r, j, init_crds_boxl_3_5, neg_x_image);
%                 
%             end
%         elseif (abs(r(2)) > abs(r(1)) && abs(r(2)) > abs(r(3)))
%             if r(2) > 0
%                 %                 disp('pos y')
%                 v = calc_LJ_Potential(r, j, init_crds_boxl_3_5, pos_y_image);
%             else
%                 %                 disp('neg y')
%                 v = calc_LJ_Potential(r, j, init_crds_boxl_3_5, neg_y_image);
%             end
%         else
%             if r(3) > 0
%                 %                 disp('pos z')
%                 v = calc_LJ_Potential(r, j, init_crds_boxl_3_5, pos_z_image);
%             else
%                 %                 disp('neg z')
%                 v = calc_LJ_Potential(r, j, init_crds_boxl_3_5, neg_z_image);
%             end
%         end
%         V = V + v;
%     end
% end
% disp(V)


% function potential = calcPotential(ri, j, mat1, mat2)
%     epsilon = 0.25;
%     sigma = 1;
%     L = 3.5; %Length of box
%
%     rj = [mat1(j,1), mat1(j,2), mat1(j,3)];
%     dist = sqrt(sum((ri - rj) .^ 2));
%
%     if dist > (L/2)
%
%     end
%
% %     if j <= 30
% %         rj = [mat1(j,1), mat1(j,2), mat1(j,3)];
% %     else
% %         rj = [mat2(j-30,1), mat2(j-30,2), mat2(j-30,2)];
% %     end
%     if dist
%         potential = 4*epsilon*( (sigma/dist)^12 - (sigma/dist)^6 );
%     else
%         potential = 0;
%     end
% end

%%
load('init_crds_boxl_3.5.dat')

%Initialize Variables
L = 3.5; %Length of box
rho_star = 0.7; %Reduced density
epsilon = 0.25;
sigma = 1;


%Identify the nearest periodic boundary for each point
V = 0;
for i = 1:30
    disp(V)
    ri = [init_crds_boxl_3_5(i, 1), init_crds_boxl_3_5(i, 2), init_crds_boxl_3_5(i, 3)];
    for j = i+1:30
        rj = [init_crds_boxl_3_5(j, 1), init_crds_boxl_3_5(j, 2), init_crds_boxl_3_5(j, 3)];
        
        v = calc_LJ_Potential(L, epsilon, sigma, ri, rj);
        V = V + v;
    end
end
disp(V)


function v = calc_LJ_Potential(L, epsilon, sigma, ri, rj)

Deltx = ri(1) - rj(1);
    Delty = ri(2) - rj(2);
    Deltz = ri(3) - rj(3);
    
    Deltx = Deltx - L*double(int8(Deltx/L));
    Delty = Delty - L*double(int8(Delty/L));
    Deltz = Deltz - L*double(int8(Deltz/L));
    
    dist = sqrt(Deltx^2 + Delty^2 + Deltz^2);

    v = 4*epsilon*( (sigma/dist)^12 - (sigma/dist)^6 );
end

