function [k_p, C, phi] = Gen_Cspace(size_M)
% GEN_CSPACE is used to generate coherent target space
% and return matrix, C, and corresponding k_p
	
	rng(99); % number of random generator
    %size_M = 100; % Numbers of vector in C
	% generate coherency target space (C)
	mu = 0;
    sigma = 1;
    phi = zeros(3,size_M); % Unit: degree, (NOT RADIAN!)
    C = zeros(9, size_M);
	k_p = normrnd(mu,sigma,[3, size_M]) + 1j*normrnd(mu,sigma,[3, size_M]);
    k_p = k_p./repmat(sqrt(sum(abs(k_p).^2, 1)).*exp(-1j*angle(k_p(1,:))), [3,1]);
    phi(1,:) = acosd(real(k_p(1,:)));
    phi(2,:) = acosd(real(k_p(2,:))./sind(phi(1,:)));
    phi(3,:) = acosd(imag(k_p(2,:))./(sind(phi(1,:)).*sind(phi(2,:))));
    phi(4,:) = acosd(real(k_p(3,:))./(sind(phi(1,:)).*sind(phi(2,:)).*sind(phi(3,:)) ));
    for r = 1 : size_M
        t = k_p(:,r)*k_p(:,r)';
        C(:,r) = [real(t(1,1)); real(t(2,2)); real(t(3,3)); sqrt(2)*real(t(1,2)); sqrt(2)*imag(t(1,2)); ...
                sqrt(2)*real(t(1,3)); sqrt(2)*imag(t(1,3)); sqrt(2)*real(t(2,3)); sqrt(2)*imag(t(2,3))];
    end
	% Visulization of the coherent target space
    Vis(k_p,'C', 'C_2D',false) 
end