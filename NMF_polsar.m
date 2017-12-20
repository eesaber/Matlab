%% The code implement decomposition by using non-negative matrix factorization
clear 
clc
%% Simulated data
% declaration of variable
q = 9; 
M = 100;
k = 8;
N_az = 100; N_ra = 100; N = N_az*N_ra; % Image with size (100x100)
%C = random(a,M); % Coherent target space (qxM)
%A = random(M,k); % Assembly matrix (Mxk)
%X = random(k,N); % Spatial distribution matrix (kxN)

% Plot the sptial distribution 
map = zeros(N_az, N_ra, k);
rng(1024)
t = randi([0 1],100);
for rr = 1 : k
    figure(k)
        imagesc(map(:,:,rr))
        %xlabel('range')
        %ylabel('azimuth')
        title(['# of ' num2str(k)])
end
pause
close all