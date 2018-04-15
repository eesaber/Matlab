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
    area = 'area_3';
	[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('Test', area);
	[N_az, N_ra] = size(hh_hh);
	size_N = numel(hh_hh);
    size_Q = 6; % Numbers of atom in D
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
    % 4-component decomposition
    span = hh_hh+vv_vv+2*hv_hv;
    FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp')
	clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
end
clear simulation
temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
%% Pauli decomposition
figure
Pauli_decomp(reshape(Y(2,:),[N_az, N_ra]), reshape(Y(3,:),[N_az, N_ra]),......
    reshape(Y(1,:),[N_az, N_ra]),'Filename','Pauli_decomp','saibu',false)
%Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv))/2,'Filename','Pauli_decomp')
close all
%% H_alpha decomposition 
ind_ = randperm(size_N);
figure
H_Alpha(temp_T(:,:,ind_(1:2000)))
clear ind_

% Plot each entropy of each pixel
t = cputime;
[~,~,num]= size(temp_T);
H = zeros(1, num);
for r = 1 : num
    [U, L] = eig(temp_T(:,:,r));
    P = L/sum(L);
    H(r) = -sum(P.*log(P)/log(3));
    %alpha = sum(P'.*acosd(abs(U(1,:))));
end
e = cputime -t;
disp(e)
H  = reshape(H, size(T_11));
close all
line = zeros(1,size(H,1));
global dir;
load([dir area '_hgt.mat']);
for a=1 : size(H,1)
    line(a) = find(hgt(a,:)+10000,1);
end
figure
    imagesc(H)
    set(gca,'Ydir','normal','Clim',Clim, 'View',[90 90], 'XTick',1:300:1401)
    xt=arrayfun(@num2str,get(gca,'xtick')+600-1,'un',0);
    yt=arrayfun(@num2str,get(gca,'ytick')+26000,'un',0);
    set(gca,'xticklabel',xt,'yticklabel',yt)
    colormap jet; colorbar
    hold on 
    plot(line,1:size(H,1),'magenta', 'Linewidth', 2 )
    hold off
    xlabel('east')
    ylabel('north')
    plot_para('Filename','output/Entropy', 'Maximize',true)
%% Obtain terrain slope in azimuth and induced angle by terrain slope. 

close all
thm_angle = Find_angle(temp_T);
clear T_11 T_22 T_33 T_12 T_13 T_23 
%clear temp_T 
thm_angle = resize(thm_angle, [1,size_N]);

%% Generate the redundant coding matrix 
disp('Generating R...')
size_M = 100;
[k_p, C, phi] = Gen_Cspace(size_M);%generate another coherency target space (C)
%%
R = single(zeros(size_M, size_N));
for n = 1 : size_N
    T_n = [Y(1,n), (Y(4,n)+1j*Y(5,n))/sqrt(2), (Y(6,n)+1j*Y(7,n))/sqrt(2);.....
            conj((Y(4,n)+1j*Y(5,n))/sqrt(2)), Y(2,n), (Y(8,n)+1j*Y(9,n))/sqrt(2);......
            conj((Y(6,n)+1j*Y(7,n))/sqrt(2)), conj((Y(8,n)+1j*Y(9,n))/sqrt(2)), Y(3,n)];
    
    kerl = inv(T_n/trace(T_n));
	for m = 1 : size_M
        %R(m,n) = real(-k_p(:,m)'*kerl*k_p(:,m))^(-4);
        R(m,n) = real(k_p(:,m)'*kerl*k_p(:,m))^(-7/2);
	end
	R(:,n) = R(:,n)/sum(R(:,n))*trace(T_n);
    if mod(n,1e5) == 0
        disp(n)
    end
end
%Compare to the original 
Y_cod = C*R;
        
Pauli_decomp(reshape(Y_cod(2,:),[N_az, N_ra]), reshape(Y_cod(3,:),[N_az, N_ra]),......
    reshape(Y_cod(1,:),[N_az, N_ra]),'Filename', 'Result_Redund');
mse_msg = {'|S_hh + S_vv|^2', '|S_hh - S_vv|^2', '|S_hv|^2'};
disp('Comparison between the encoded image and the received image.')
for n = 1 : 3
    err_Y = 1/numel(Y_cod(n,:))*sum((Y_cod(n,:) - Y(n,:)).^2);
    fprintf(['mean of            ', char(mse_msg(n)), ' : %f, '], mean(Y(n,:)))
    fprintf(['std of             ', char(mse_msg(n)), ' : %f \n'], std(Y(n,:)))
    fprintf(['after NMF, mean of ', char(mse_msg(n)), ' : %f, '], mean(Y_cod(n,:)))
    fprintf(['after NMF, std of  ', char(mse_msg(n)), ' : %f \n'], std(Y_cod(n,:)))
    fprintf(['MSE of             ', char(mse_msg(n)), ' : %f\n'], err_Y)
    disp('---------------------------------------------------------')
end

%% Non-negative matrix factorization
opt = statset('MaxIter', 100, 'Display', 'final', 'TolX', 1e-6, 'TolFun', 1e-6,'UseParallel', false);
fprintf('Q = %i, NMF...\n', size_Q)

rng('shuffle') % random seed
%w_0 = single(rand(size_M, size_Q)>0.5);
%x_0 = single(rand(size_Q, size_N)>0.5);
%[A_sol, X_sol] = nnmf(R, size_Q,'w0',w_0, 'h0',x_0, 'algorithm', 'mult', 'options', opt);
[A_sol, X_sol] = nnmf(R, size_Q, 'algorithm', 'mult', 'options', opt);

%% Compare to the original 
Y_sol = C*A_sol*X_sol;
%{
T_temp = cat(1,cat(2, reshape(Y_sol(1,:),[1,1,size_N]), reshape((Y_sol(4,:)+1j*Y_sol(5,:))/sqrt(2),[1,1,size_N]), reshape((Y_sol(6,:)+1j*Y_sol(7,:))/sqrt(2),[1,1,size_N])), ......
             cat(2,reshape((Y_sol(4,:)-1j*Y_sol(5,:))/sqrt(2),[1,1,size_N]), reshape(Y_sol(2,:),[1,1,size_N]), reshape((Y_sol(8,:)+1j*Y_sol(9,:))/sqrt(2),[1,1,size_N])),.......
             cat(2,reshape((Y_sol(6,:)-1j*Y_sol(7,:))/sqrt(2),[1,1,size_N]), reshape((Y_sol(8,:)-1j*Y_sol(9,:))/sqrt(2),[1,1,size_N]), reshape(Y_sol(3,:),[1,1,size_N])));
Find_angle(T_temp)
%}
figure
Pauli_decomp(reshape(Y_sol(2,:),[N_az, N_ra]), reshape(Y_sol(3,:),[N_az, N_ra]),......
    reshape(Y_sol(1,:),[N_az, N_ra]),'Filename', 'Result_aftNMF')
D_sol = C*A_sol;


disp('Comparison between the reconstructed image and the received image.')
mse_msg = {'|S_hh + S_vv|^2', '|S_hh - S_vv|^2', '|S_hv|^2'};
for n = 1 : 3
    err_Y = 1/numel(Y_sol(n,:))*sum((Y_sol(n,:) - Y(n,:)).^2);
    fprintf(['mean of            ', char(mse_msg(n)), ' : %f, '], mean(Y(n,:)))
    fprintf(['std of             ', char(mse_msg(n)), ' : %f \n'], std(Y(n,:)))
    fprintf(['after NMF, mean of ', char(mse_msg(n)), ' : %f, '], mean(Y_sol(n,:)))
    fprintf(['after NMF, std of  ', char(mse_msg(n)), ' : %f \n'], std(Y_sol(n,:)))
    fprintf(['MSE of             ', char(mse_msg(n)), ' : %f\n'], err_Y)
    disp('---------------------------------------------------------')
end

%% Visualize dictionary
D_T = cat(1,cat(2, reshape(D_sol(1,:),[1,1,size_Q]), reshape((D_sol(4,:)+1j*D_sol(5,:))/sqrt(2),[1,1,size_Q]), reshape((D_sol(6,:)+1j*D_sol(7,:))/sqrt(2),[1,1,size_Q])), ......
             cat(2,reshape((D_sol(4,:)-1j*D_sol(5,:))/sqrt(2),[1,1,size_Q]), reshape(D_sol(2,:),[1,1,size_Q]), reshape((D_sol(8,:)+1j*D_sol(9,:))/sqrt(2),[1,1,size_Q])),.......
             cat(2,reshape((D_sol(6,:)-1j*D_sol(7,:))/sqrt(2),[1,1,size_Q]), reshape((D_sol(8,:)-1j*D_sol(9,:))/sqrt(2),[1,1,size_Q]), reshape(D_sol(3,:),[1,1,size_Q])));
Vis([],'A','T',D_T)

% Distribution map of atoms 
    % Aboslute value 
X_map = reshape(X_sol', [N_az, N_ra, size_Q]);
alphabets = char(97:96+size_Q).'; % Label for subplot, start from a to ...
subplot_label = strcat({'('}, alphabets, {')'});
X_avg = mean(X_map, 3);
for rr = 1 : size_Q
	figure
	%subplot(2, size_Q/2, rr)
    %subplot(1, size_Q, rr)
	%imagesc(X_map(:, :, rr), [0, 0.0001])
    imagesc(10*log10(span.*X_map(:, :, rr)./X_avg))
	ylabel(subplot_label(rr,:),'Interpreter', 'latex')
    colormap jet; colorbar
    set(gca,'Ydir','normal','Clim',[-30 10], 'View',[90 90], 'XTick',1:300:1401)
        xt=arrayfun(@num2str,get(gca,'xtick')+600-1,'un',0);
        yt=arrayfun(@num2str,get(gca,'ytick')+26000,'un',0);
    set(gca,'xticklabel',xt,'yticklabel',yt)
    %plot_para('Fontsize',22,'Ratio',[1,1,1]);
    %axis off
	plot_para('Fontsize',40,'Filename',['output/X_map_Q_', int2str(size_Q),'_',char(96+rr)],'Maximize',true);
end
%plot_para('Fontsize',22,'Filename',['X_map_Q_', int2str(size_Q)],'Maximize',true);
%movefile(['X_map_Q_', int2str(size_Q),'.jpg'], 'output')

    % Pauli decomposition of \bar{\bar{D}} \cdot \bar{\bar{X}}
subplot_label = char(97:96+size_Q).';
for rr = 1 : size_Q
	map_q = D_sol(:,rr)*X_sol(rr,:);
    figure
	Pauli_decomp(reshape(map_q(2,:),[N_az, N_ra]), 2*reshape(map_q(3,:),[N_az, N_ra]),......
    reshape(map_q(1,:),[N_az, N_ra]),'Filename', ['Map_Pauli_Q_' int2str(size_Q) '_Atom_', subplot_label(rr)])
end