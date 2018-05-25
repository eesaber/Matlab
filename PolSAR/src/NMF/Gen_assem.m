function [sigma] = Gen_assem(varargin)
% GEN_ASSEM is used to generate dictionary matrix for simulation.
% The return of GEN_ASSEM is an 3D matrix of size [row, column, # of
% atom]. 
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'logical'},{});
	validationFcn_2_ = @(x) validateattributes(x,{'logical'},{});
    validationFcn_3_ = @(x) validateattributes(x,{'logical'},{});
	addParameter(parse_,'Isplot',true,validationFcn_1_);
	addParameter(parse_,'SubPlot',false ,validationFcn_2_);    
	addParameter(parse_,'TwoD',true, validationFcn_3_);
	parse(parse_,varargin{:})

	rng(99); % seed of random number generator
	size_Q = 8; % number of atoms
	k_p_mat = zeros(3,size_Q);
    sigma = zeros(3,3,size_Q);
	delta = 0.00001; % make the \bar{\bar{T}} become full rank.
	% surface 3 atoms
	for k = 1 : 3
		S = rand(2).*eye(2);%.*exp(j*2*pi*rand(2));
		k_p = 1/sqrt(2)*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
		k_p = k_p/norm(k_p);
		psi = 120*rand(1)-60;
		R_k = [1, 0, 0; 0, cosd(2*psi), -sind(2*psi); 0, sind(2*psi), cosd(2*psi)];
		k_p_mat(:,k) = R_k*k_p;
		sigma(:,:,k) = k_p_mat(:,k)*k_p_mat(:,k)';
		for q = 1:2
			temp = diag(sigma(:,:,k));
			[~, ind]= sort(temp);
			temp(ind<3) = temp(ind<3) + delta;
			sigma(:,:,k) = sigma(:,:,k) - diag(diag(sigma(:,:,k))) + temp;
		end
		%{
		[V, D]= eig(k_p*k_p');
		D(diag((diag(D)<0.0001))) = 0.0001;
		k_p(:,k) = 1/sqrt(2)*R_k*k_p;
        sigma(:,:,k) = R_k*V*D*inv(V)*R_k.';
		%}
		%{
		disp('-----------------------------')
		disp(sigma(:,:,k))
		disp(R_k*(k_p*k_p')*R_k.')
		disp(inv(sigma(:,:,k)))
		disp(eig(sigma(:,:,k)))
		disp(rank(sigma(:,:,k)))
		%}
	end
	% 1st and 2nd are narrow dihedral (z = -1/2)
	% 3rd and 4th are dihedral (z = -1)
	S_temp = [2,0,0,-1;2,0,0,-1;1,0,0,-1;1,0,0,-1];
	for k = 4 : 7
		S = reshape(S_temp(k-3,:),[2,2]);
		k_p = 1/sqrt(2)*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
		k_p = k_p/norm(k_p);
		theta = 180*rand(1)-90;
		R_k = [1, 0, 0; 0, cosd(2*theta), -sind(2*theta); 0, sind(2*theta), cosd(2*theta)];
		k_p_mat(:,k) = R_k*k_p;
        sigma(:,:,k) = k_p_mat(:,k)*k_p_mat(:,k)';
		for q = 1:2
			temp = diag(sigma(:,:,k));
			[~, ind]= sort(temp);
			temp(ind<3) = temp(ind<3) + delta;
			sigma(:,:,k) = sigma(:,:,k) - diag(diag(sigma(:,:,k))) + temp;
		end
	end

	% volume
    % horizontal dipole S = [1, 0; 0, 0]
	phi = linspace(0,180,362); % Orientation angle
    for  k = 1 : numel(phi)
		temp =  1/sqrt(2)*[1; -cosd(2*phi(k)); -sind(2*phi(k))];
		k_p_mat(:,8) = k_p_mat(:,8) + temp*abs(cosd(phi(k)))/2*(phi(2)-phi(1))/180*pi;
		sigma(:,:,8) = sigma(:,:,8) + temp*temp'*abs(cosd(phi(k)))/2*(phi(2)-phi(1))/180*pi;
        %sigma(:,:,8) = sigma(:,:,8);
    end
	
	% Check if there exist invers matrix of A? i.e. det(A) ≠ 0
	for k = 1 : size_Q
        %disp(det(sigma(:,:,k)))
		%sigma(:,:,k) = sigma(:,:,k) + 10^-2*(abs(sigma(:,:,k)) < 0.0001).*eye(3);
		%sigma(:,:,k) = diag(diag(sigma(:,:,k)));
		%sigma(:,:,k) = sigma(:,:,k)/10;
	end
	
	%% plot
	if parse_.Results.Isplot 
        close all
		Vis(k_p, 'A', 'A_Subplot', parse_.Results.SubPlot, ......
				'A_2D', parse_.Results.TwoD)
	end
end