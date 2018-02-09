%% The code implement decomposition by using non-negative matrix factorization
clear 
clc
simulation = 1;
%% import or simulate data
if simulation
	% Simulated data
    disp('Using simulation data')
	% declaration of variable
    rng(181); % number of random generator
	N_az = 100; N_ra = 100; % az==row, ran == col
	N = N_az*N_ra; % Image with size (100x100)
    size_Q = 8; % Numbers of atom in D
    size_M = 1000; % Numbers of vector in C
    %% generate coherency target space (C)
    mu = 0;
    sigma = 1;
    k_p = normrnd(mu,sigma,[3, size_M]) + 1j*normrnd(mu,sigma,[3, size_M]);
    C = zeros(9, size_M);
    for r = 1 : size_M
        %k_p(:,r) = exp(-1j*angle(k_p(1,r)))*k_p(:,r);
        k_p(:,r) = k_p(:,r)/norm(k_p(:,r),2);
        t = k_p(:,r)*k_p(:,r)';
        C(:,r) = [real(t(1,1)); real(t(2,2)); real(t(3,3)); sqrt(2)*real(t(1,2)); sqrt(2)*imag(t(1,2)); ...
                sqrt(2)*real(t(1,3)); sqrt(2)*imag(t(1,3)); sqrt(2)*real(t(2,3)); sqrt(2)*imag(t(2,3))];
    end
    Vis_Co(k_p,'ThreeD',false) % Visulization of the coherent target space
    %% generate spatial distribution (X)
    X = Gen_map(N_az, N_ra, size_Q); % Spatial distribution matrix (QxN)
    %% generate assembly matrix (A) MxQ
    sigma = Gen_Scatterer();
    A = zeros(size_M, size_Q);
    for q = 1 : size_Q
        A(:,q) = 1/(pi^2.5*det(sigma(:,:,q)))*exp(-real(sum((k_p'*sigma(:,:,q)).*k_p.',2)));
    end
    %% Synthesis the data (Y)
	Y = C*A*X;
else
	% Real data import
    disp('Using real data...')
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
[A_sol, X_sol] = nnmf(Y, size_Q, 'Display', 'iter', 'UseParallel', true);
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
