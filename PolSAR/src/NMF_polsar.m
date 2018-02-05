%% The code implement decomposition by using non-negative matrix factorization
clear 
clc
test = 1;
q = 9; 
M = 100;
Q = 8;
%% Read Data
if test 
	% Simulated data
	% declaration of variable
	N_az = 100; N_ra = 100; % az==row, ran == col
	N = N_az*N_ra; % Image with size (100x100) 
	%A = random(M,Q); % Assembly matrix (MxQ)
	X = gen_map(N_az, N_ra, Q); % Spatial distribution matrix (QxN)
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
%% Build coherency target space 
%{
rho = linspace(0,1,101);
theta = linspace(0,2*pi-0.001,360);
[kp_2, kp_3] = pol2cart(meshgrid(rho, theta));
kp_2 = reshape(1, numel(kp_2));
kp_3 = reshape(1, numel(kp_3));
kp_1 = sqrt(1 - kp_2.^2 - kp_3.^2);
%} 
mu = 0;
sigma = 1;
size_K = 100;
k_p = normrnd(mu,sigma,[3, size_K]) + 1j*normrnd(mu,sigma,[3, size_K]);
C = zeros(9, size_K);
for r = 1 : size_K
	%k_p(:,r) = exp(-1j*angle(k_p(1,r)))*k_p(:,r);
	k_p(:,r) = k_p(:,r)/norm(k_p(:,r),2);
	t = k_p(:,r)*k_p(:,r)';
	C(:,r) = [real(t(1,1)); real(t(2,2)); real(t(3,3)); sqrt(2)*real(t(1,2)); sqrt(2)*imag(t(1,2)); ...
			sqrt(2)*real(t(1,3)); sqrt(2)*imag(t(1,3)); sqrt(2)*real(t(2,3)); sqrt(2)*imag(t(2,3))];
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
