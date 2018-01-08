function gen_scatterer()
%function [S, Sigma] = gen_scatterer()
% GEN_SCATTERER is used to generate dictionary matrix for simulation.
% The return of GEN_SCATTERER is an 3D matrix of size [row, column, # of atom].
	rng(10);
	n_atom = 8;
	psi = rand(1)*pi;
	R = [1, 0, 0; 0, cos(2*psi), -sin(2*psi); 0, sin(2*psi), cos(2*psi)];
	
	for qq = 1 : n_atom
		subplot(2, n_atom/2, qq)
		Vis_Co(k_p(:,qq))
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