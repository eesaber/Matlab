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
	addParameter(parse_,'ThreeD',false, validationFcn_3_);
	parse(parse_,varargin{:})

	rng(99); % seed of random number generator
	size_Q = 8; % number of atoms
	k_p = zeros(3,size_Q);
    sigma = zeros(3,3,size_Q);
	% surface 3 atoms
	for k = 1 : 3
		S = abs(rand(2).*eye(2));
		psi = 120*rand(1)-60;
		R_k = [1, 0, 0; 0, cosd(2*psi), -sind(2*psi); 0, sind(2*psi), cosd(2*psi)];
        %temp = 1/sqrt(2)*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
		%temp = temp/norm(temp);
        %sigma(:,:,k) = temp*temp';
		k_p(:,k) = 1/sqrt(2)*R_k*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
        sigma(:,:,k) = k_p(:,k)*k_p(:,k)';
	end
	% 1st and 2nd are narrow dihedral (z = -1/2)
	% 3rd and 4th are dihedral (z = -1)
	S = [2,0,0,-1;2,0,0,-1;1,0,0,-1;1,0,0,-1];
	for k = 4 : 7
		S_temp = reshape(S(k-3,:),[2,2])*exp(j*2*pi*rand(1));
		theta = 180*rand(1)-90;
		R_k = [1, 0, 0; 0, cosd(2*theta), -sind(2*theta); 0, sind(2*theta), cosd(2*theta)];
		%temp = 1/sqrt(2)*[S_temp(1,1)+S_temp(2,2); S_temp(1,1)-S_temp(2,2);......
		%				2*S_temp(1,2)];
		%temp = temp/norm(temp);
		%sigma(:,:,k) = temp*temp';
		k_p(:,k) = 1/sqrt(2)*R_k*[S_temp(1,1)+S_temp(2,2);......
						S_temp(1,1)-S_temp(2,2); 2*S_temp(1,2)];
        sigma(:,:,k) = k_p(:,k)*k_p(:,k)';
	end

	% volume
    % horizontal dipole S = [1, 0; 0, 0]
	phi = linspace(0,180,362); % Orientation angle
    for  k = 1 : numel(phi)
		temp =  1/sqrt(2)*[1; -cosd(2*phi(k)); -sind(2*phi(k))];
		k_p(:,8) = k_p(:,8) + temp*abs(cosd(phi(k)))/2*(phi(2)-phi(1))/180*pi;
		sigma(:,:,8) = sigma(:,:,8) + temp*temp'*abs(cosd(phi(k)))/2*(phi(2)-phi(1))/180*pi;
        %sigma(:,:,8) = sigma(:,:,8);
    end
	
    for k = 1 : size_Q
        %sigma(:,:,k) = sigma(:,:,k) + 10^-2*(abs(sigma(:,:,k)) < 0.0001).*eye(3);
		%sigma(:,:,k) = diag(diag(sigma(:,:,k)));
		%sigma(:,:,k) = sigma(:,:,k)/10;
	end
	
	%% plot
	if parse_.Results.Isplot 
        close all
		Vis_Assem(k_p, sigma,'SubPlot', parse_.Results.SubPlot, ......
				'ThreeD', parse_.Results.ThreeD)
	end
end