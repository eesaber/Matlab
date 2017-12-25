%% The code implement decomposition by using non-negative matrix factorization
clear 
clc
test = 1;
q = 9; 
M = 100;
k = 8;
%% Read Data
if test 
	% Simulated data
	% declaration of variable
	N_az = 100; N_ra = 100; % az==row, ran == col
	N = N_az*N_ra; % Image with size (100x100) 
	%C = random(a,M); % Coherent target space (qxM)
	%A = random(M,k); % Assembly matrix (Mxk)
	X = gen_map(N_az, N_ra, k); % Spatial distribution matrix (kxN)
	D = 1;
	Y = D*X;
else
	% Real data import
	[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = data_io();
	n = numel(hh_hh);
	Y = [reshape(hh_hh, [1, n]); reshape(hv_hv, [1, n]); reshape(vv_vv, [1, n]); ......
		real(reshape(hh_hv, [1, n]))/sqrt(2); imag(reshape(hh_hv, [1, n]))/sqrt(2);
		real(reshape(hh_vv, [1, n]))/sqrt(2); imag(reshape(hh_vv, [1, n]))/sqrt(2);
		real(reshape(hv_vv, [1, n]))/sqrt(2); imag(reshape(hv_hv, [1, n]))/sqrt(2)];
	[N_az, N_ra] = Y.size;
	clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
end

%% Non-negative matrix factorization
[A_sol, X_sol] = nnmf(Y, k, 'Display', 'iter', 'UseParallel', true);
X_sol = reshape(X_sol, [N_az, N_ra, k]);
figure(1)
	for rr = 1 : k
		subplot(2, k/2, rr)
		imagesc(~map(:, :, rr))
        title(['# of ' num2str(rr)])
		set(gca, 'YDir','normal')
		axis off
        colormap gray
	end
	axis on
