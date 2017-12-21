%% The code implement decomposition by using non-negative matrix factorization
clear 
clc
%% Simulated data
% declaration of variable
q = 9; 
M = 100;
k = 8;
N_az = 100; N_ra = 100; % az==row, ran == col
N = N_az*N_ra; % Image with size (100x100) 
%C = random(a,M); % Coherent target space (qxM)
%A = random(M,k); % Assembly matrix (Mxk)
X = gen_map(k, N_az, N_ra,); % Spatial distribution matrix (kxN)


