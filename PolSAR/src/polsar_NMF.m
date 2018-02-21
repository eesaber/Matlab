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
		A(isnan(A)) = 0;
		Y = C*A*X;
	end
else
	% Real data import
    disp('Using real data...')
	[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('Test',true);
	[N_az, N_ra] = size(hh_hh);
	size_N = numel(hh_hh);
	Y = [reshape(hh_hh, [1, size_N]); reshape(hv_hv, [1, size_N]); reshape(vv_vv, [1, size_N]); ......
		real(reshape(hh_hv, [1, size_N]))/sqrt(2); imag(reshape(hh_hv, [1, size_N]))/sqrt(2);
		real(reshape(hh_vv, [1, size_N]))/sqrt(2); imag(reshape(hh_vv, [1, size_N]))/sqrt(2);
		real(reshape(hv_vv, [1, size_N]))/sqrt(2); imag(reshape(hv_hv, [1, size_N]))/sqrt(2)];
	clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
end
clear simulation
%% Pauli decomposition
if(1)	
	up_ = 10; low_ = -20;
	Pauli = ones(N_az, N_ra, 3);	
	% |S_vv - S_hh|^2 -> double bounce scattering 
	t_p = 10*log10(sqrt(reshape(Y(2,:),[N_az, N_ra])));	
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,1) = (t_p-low_)/(up_-low_);	
	% |S_hv|^2 -> volume scattering
	t_p= 10*log10(sqrt(reshape(Y(3,:),[N_az, N_ra])));
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,2) = (t_p-low_)/(up_-low_);
	% |S_vv + S_hh|^2 -> single scattering
	t_p = 10*log10(sqrt(reshape(Y(1,:),[N_az, N_ra])));	
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,3) = (t_p-low_)/(up_-low_);
	figure(7)
		image(Pauli)
		set(gca,'Ydir','normal')
		xlabel('Azimuth (m)', 'Fontsize', 40)
		ylabel('Range (m)', 'Fontsize', 40)
		plot_para('Filename','My_pauli_c','Maximize',true, 'Ratio', [4 3 1]);
end
%% Generate the redundant coding matrix 
size_M = 100;
[k_p, C, phi] = Gen_Cspace(size_M);%generate another coherency target space (C)
R = zeros(size_M, size_N);
for m = 1 : size_M
    for n = 1 : size_N
        T_n = [Y(1,n), sqrt(2)*(Y(4,n)+1j*Y(5,n)), sqrt(2)*(Y(6,n)+1j*Y(7,n));.....
            conj(sqrt(2)*(Y(4,n)+1j*Y(5,n))), Y(2,n), sqrt(2)*(Y(8,n)+1j*Y(9,n));......
            conj(sqrt(2)*(Y(6,n)+1j*Y(7,n))), conj(sqrt(2)*(Y(8,n)+1j*Y(9,n))), Y(3,n)];
        R(m,n) = exp(k_p(:,m)'*(T_n\k_p(:,m)/sum(diag(T_n))))^(-7/2);
    end
end

%% Non-negative matrix factorization
opt = statset('MaxIter', 5, 'Display', 'iter', 'UseParallel', true);
[A_sol, X_sol] = nnmf(R, size_Q,'algorithm', 'als', 'options', opt);
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

