%% The code implements PolSAR decomposition by using non-negative matrix factorization
clear 
clc
simulation = 0;

%% import or simulate data
if simulation
	% Simulated data
    disp('Using simulation data...')
	% declaration of variable
    rng(99); % number of random generator
	N_az = 100; N_ra = 100; % az==row, ran == col
	size_N = N_az*N_ra; % Image with size (100x100)
    size_Q = 8; % Numbers of atom in D
    size_M = 100; % Numbers of vector in C
    %% generate coherency target space (C)
    [k_p, C, phi] = Gen_Cspace(size_M);

    %% generate assembly matrix (A) MxQ
    sigma = Gen_assem('Isplot',true);
    A = zeros(size_M, size_Q);
    for q = 1 : size_Q 
		temp = inv(sigma(:,:,q));
        for m = 1 : size_M
            A(m,q) = 1/(pi^2.5*det(sigma(:,:,q)))*exp(-k_p(:,m)'*temp*k_p(:,m));
        end
        A(:,q) = A(:,q)/sum(A(:,q));
    end
  
    %% generate spatial distribution (X)
    X = zeros(size_Q, size_N);
    temp = Gen_map(N_az, N_ra, size_Q); % Spatial distribution matrix (QxN)
    for k = 1 : size_Q
        X(k,:) = reshape(temp(:,:,k),[1, size_N]);
    end
    clear temp

    %% Synthesis the data (Y)
	Y = C*A*X;
	if (sum(sum(isnan(Y))) > 0) || (sum(sum(isfinite(Y))) > 0)
        disp('There are NaN or Inf in A')
		A( isnan(A)) = 0;
		Y = C*A*X;
	end
else
	% Real data import
    disp('Using real data...')
	[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('Test',true);
	[N_az, N_ra] = size(hh_hh);
	size_N = numel(hh_hh);
    size_Q = 4; % Numbers of atom in D
	T_11 = (hh_hh+vv_vv+hh_vv+conj(hh_vv))/2;
	T_22 = (hh_hh+vv_vv-hh_vv-conj(hh_vv))/2;
	T_33 = 2*hv_hv;
	T_12 = (hh_hh-vv_vv-hh_vv+conj(hh_vv))/2;
	T_13 = hh_hv + conj(hv_vv);
    T_23 = hh_hv - conj(hv_vv);
	Y = [reshape(T_11, [1, size_N]); reshape(T_22, [1, size_N]); reshape(T_33, [1, size_N]); ......
		reshape(real(T_12), [1, size_N])*sqrt(2); reshape(imag(T_12), [1, size_N])*sqrt(2);
		reshape(real(T_13), [1, size_N])*sqrt(2); reshape(imag(T_13), [1, size_N])*sqrt(2);
		reshape(real(T_23), [1, size_N])*sqrt(2); reshape(imag(T_23), [1, size_N])*sqrt(2)];
	clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
end
clear simulation
%Pauli decomposition
Pauli_decomp(reshape(Y(2,:),[N_az, N_ra]), reshape(Y(3,:),[N_az, N_ra]),......
    reshape(Y(1,:),[N_az, N_ra]), 'Pauli_decomp') 
%% Generate the redundant coding matrix 
disp('Generating R...')
size_M = 100;
[k_p, C, phi] = Gen_Cspace(size_M);%generate another coherency target space (C)
%%
R = zeros(size_M, size_N);
for n = 1 : size_N
    T_n = [Y(1,n), (Y(4,n)+1j*Y(5,n))/sqrt(2), (Y(6,n)+1j*Y(7,n))/sqrt(2);.....
            conj((Y(4,n)+1j*Y(5,n))/sqrt(2)), Y(2,n), (Y(8,n)+1j*Y(9,n))/sqrt(2);......
            conj((Y(6,n)+1j*Y(7,n))/sqrt(2)), conj((Y(8,n)+1j*Y(9,n))/sqrt(2)), Y(3,n)];
    kerl = inv(T_n/trace(T_n));
	for m = 1 : size_M
        %R(m,n) = real(-k_p(:,m)'*kerl*k_p(:,m))^(-4);
        R(m,n) = real(k_p(:,m)'*kerl*k_p(:,m))^(-7/2);
	end
	R(:,n) = R(:,n)/(Y(1,n)+Y(2,n)+Y(3,n));
end
% Compare to the original 
Y_sol = C*R;
Pauli_decomp(reshape(Y_sol(2,:),[N_az, N_ra]), reshape(Y_sol(3,:),[N_az, N_ra]),......
    reshape(Y_sol(1,:),[N_az, N_ra]), 'redun_re')

%% Non-negative matrix factorization
disp('NMF...')
opt = statset('MaxIter', 100, 'Display', 'iter', 'UseParallel', true);
size_Q = 8;
[A_sol, X_sol] = nnmf(R, size_Q,'algorithm', 'als', 'options', opt);

% Compare to the original 
Y_sol = C*A_sol*X_sol;
Pauli_decomp(reshape(Y_sol(2,:),[N_az, N_ra]), reshape(Y_sol(3,:),[N_az, N_ra]),......
    reshape(Y_sol(1,:),[N_az, N_ra]), 'result')

% distribution map
X_sol = reshape(X_sol, [N_az, N_ra, size_Q]);
figure
for rr = 1 : size_Q
	subplot(2, size_Q/2, rr)
    %subplot(1, size_Q, rr)
	imagesc(X_sol(:, :, rr))
	title(['# of ' num2str(rr)])
	set(gca, 'YDir','normal')
	%axis off
	colormap gray
end
plot_para('Fontsize', 12,'Filename','X_map','Maximize',true);
movefile X_map.jpg output/
