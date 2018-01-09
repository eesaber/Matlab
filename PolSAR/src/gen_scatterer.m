function gen_scatterer()
%function [S, Sigma] = gen_scatterer()
% GEN_SCATTERER is used to generate dictionary matrix for simulation.
% The return of GEN_SCATTERER is an 3D matrix of size [row, column, # of atom].
	rng(77);
	n_atom = 8;
	range = 10;
	phi = linspace(0,pi,180);
	R = reshape([cos(phi); sin(phi); -sin(phi); cos(phi)], [2,2,numel(phi)]);
	R_i = reshape([cos(phi); -sin(phi); sin(phi); cos(phi)], [2,2,numel(phi)]);
	A = zeros(3,8);
	I = [1, 0; 0, 1];
	% surface 3 atoms
	for k = 1 : 3
		S = range*rand(2).*I;
		theta = rand(1)*pi;
		R_k = [1, 0, 0; 0, cos(2*theta), -sin(2*theta); 0, sin(2*theta), cos(2*theta)];
		A(:,k) = 1/sqrt(2)*R_k*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
	end
	% dihedral 2 big sigma 2 small sigma
	for k = 4 : 7
		S = range*(1-randn(2) + 1j-1j*randn(2)).*I;
		theta = rand(1)*pi;
		R_k = [1, 0, 0; 0, cos(2*theta), -sin(2*theta); 0, sin(2*theta), cos(2*theta)];
		A(:,k) = 1/sqrt(2)*R_k*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
	end
	% volume
	S = reshape([cos(phi); zeros(size(phi)); zeros(size(phi)); zeros(size(phi))], [2,2,numel(phi)]);
	S = sum(R_i.*S.*R,3)*(phi(2)-phi(1));
	A(:,8) = 1/sqrt(2)*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
	
	for qq = 1 : n_atom
		subplot(2, n_atom/2, qq)
		Vis_Co(A(:,qq),'xlabel','','ylabel')
		set(gca, 'YDir','normal')
        colormap gray
        plot_para('Fontsize', 24,'Ratio', [1 1 1])
	end
	axis on
    plot_para('Fontsize', 24,'Ratio', [1 1 1], 'Maximize',true, 'Filename','SimMap')
	%{
	% All atom random
	rng(97); % control random seed
	range = 10; % -range/2 to range/2
	S = (-range+range*rand([2,2,n_atom])) + 1j*(-range+range*rand([2,2,n_atom]));
	S(2,1,:) = S(1,2,:);
	k_p = 1/sqrt(2)*reshape([S(1,1,:)+S(2,2,:); S(1,1,:)-S(2,2,:); 2*S(1,2,:)], [n_atom,3]).';
	for qq = 1 : n_atom
		Vis_Co(k_p(:,qq))
	end
	Sigma = zeros(2,2,n_atom);
	%}
end