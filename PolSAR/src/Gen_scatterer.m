function [k_p, sigma] = Gen_scatterer()
%function Gen_scatterer()
% GEN_SCATTERER is used to generate dictionary matrix for simulation.
% The return of GEN_SCATTERER is an 3D matrix of size [row, column, # of atom].
	rng(182); % seed of random number generator
	size_Q = 8; % number of atoms
	range = 10; 
	phi = linspace(0,180,362);
	k_p = zeros(3,size_Q);
    sigma = zeros(3,3,size_Q);
	% surface 3 atoms
	for k = 1 : 3
		S = abs(rand(2).*eye(2));
		theta = 120*rand(1)-60;
		R_k = [1, 0, 0; 0, cosd(2*theta), -sind(2*theta); 0, sind(2*theta), cosd(2*theta)];
        temp = 1/sqrt(2)*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
		temp = temp/norm(temp);
        sigma(:,:,k) = temp*temp';
		k_p(:,k) = 1/sqrt(2)*R_k*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
	end
	% dihedral 2 big sigma 2 small sigma
	for k = 4 : 7
		S = (2*randn(2)-1 + j*(2*randn(2)-1)).*eye(2);
		theta = 180*rand(1)-90;
		R_k = [1, 0, 0; 0, cosd(2*theta), -sind(2*theta); 0, sind(2*theta), cosd(2*theta)];
		temp = 1/sqrt(2)*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
		temp = temp/norm(temp);
		sigma(:,:,k) = temp*temp';
		k_p(:,k) = 1/sqrt(2)*R_k*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
        %sigma(:,:,k) = abs(k_p(:,k)*k_p(:,k)');
	end
	% volume
    Base = [1, 0; 0, 0]; % horizontal dipole
    S = zeros(2,2);
    for  k = 1 : numel(phi)
        R_i = [cosd(phi(k)), sind(phi(k)); -sind(phi(k)), cosd(phi(k))];
        %S = S + sind(phi(k))*R_i*Base*R_i.'*(phi(2)-phi(1))/180*pi;
		%S = S + R_i*Base*R_i.';
		S = R_i*Base*R_i.';
		temp =  1/sqrt(2)*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
		k_p(:,8) = k_p(:,8) + temp*sind(phi(k))/2*(phi(2)-phi(1))/180*pi;
		sigma(:,:,8) = sigma(:,:,8) + temp*temp'*sind(phi(k))/2*(phi(2)-phi(1))/180*pi;
        %sigma(:,:,8) = sigma(:,:,8)/100;
    end
	
    for k = 1 : size_Q
        sigma(:,:,k) = sigma(:,:,k) + 10^-2*(abs(sigma(:,:,k)) < 0.0001).*eye(3);
		%sigma(:,:,k) = diag(diag(sigma(:,:,k)));
    end

	%% plot
	if 1 
        close all
		Vis_Assem(k_p, sigma,'SubPlot',true)
	end
end