function gen_scatterer(varargin)
%function [S, Sigma] = gen_scatterer(varargin)
% GEN_SCATTERER is used to generate dictionary matrix for simulation.
% The return of GEN_SCATTERER is an 3D matrix of size [row, column, # of atom].
% How to call? Variable is GEN_SCATTERER(number of atom).
	switch nargin 
    	case 0
	        n_atom = 8; 
	    case 1
	    	optargs = {eps 17 @magic};
        	optargs(1) = varargin;
        	n_atom = optargs{:};
        otherwise
        	error('INPUT VARIABLE ERROR')
	end
	rng(97); % control random seed
	range = 10; % -range/2 to range/2
	S = (-range+range*rand([2,2,n_atom])) + 1j*(-range+range*rand([2,2,n_atom]));
	S(2,1,:) = S(1,2,:);
	k_p = 1/sqrt(2)*reshape([S(1,1,:)+S(2,2,:); S(1,1,:)-S(2,2,:); 2*S(1,2,:)], [n_atom,3]).';
	for qq = 1 : n_atom
		Vis_Co(k_p(:,qq))
	end
	Sigma = zeros(2,2,n_atom);
	
end