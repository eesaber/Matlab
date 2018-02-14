function [k_p, C] = Gen_Cspace(size_M)
% GEN_CSPACE is used to generate coherent target space
% and return matrix, C, and corresponding k_p
	
	rng(99); % number of random generator
    size_M = 100; % Numbers of vector in C
	% generate coherency target space (C)
	mu = 0;
    sigma = 1;
	k_p = normrnd(mu,sigma,[3, size_M]) + 1j*normrnd(mu,sigma,[3, size_M]);
    phi = zeros(3,size_M);
    C = zeros(9, size_M);
    for r = 1 : size_M
        k_p(:,r) = k_p(:,r)/norm(k_p(:,r),2)*exp(-j*atan2(imag(k_p(1,r), real(k_p(1,r)))));
        phi(1,r) = acosd(abs(k_p(:,1)));
        phi(2,r) = 
        t = k_p(:,r)*k_p(:,r)';
        C(:,r) = [real(t(1,1)); real(t(2,2)); real(t(3,3)); sqrt(2)*real(t(1,2)); sqrt(2)*imag(t(1,2)); ...
                sqrt(2)*real(t(1,3)); sqrt(2)*imag(t(1,3)); sqrt(2)*real(t(2,3)); sqrt(2)*imag(t(2,3))];
    end
	% Visulization of the coherent target space
    Vis_Co(k_p,'ThreeD',false) 
end