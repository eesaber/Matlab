function Vis_Co(k_p)
    % Generate hemisphere.
	SP_NUM = 300;
	thetavec = linspace(0,pi/2,SP_NUM/2);
	phivec = linspace(0,2*pi,2*SP_NUM);
	[th, ph] = meshgrid(thetavec,phivec);
	R = ones(size(th)); % should be your R(theta,phi) surface in general
	x = R.*sin(th).*cos(ph);
	y = R.*sin(th).*sin(ph);
	z = R.*cos(th);
   
    % Testing
    test = 1;
    if test
        psi = 90*(pi/180);
        R = [cos(psi), sin(psi); -sin(psi), cos(psi)];
        S = R.'*[1, 0; 0, 0]*R
        k_p = 1/sqrt(2)*[S(1)+S(2), S(1)-S(2), 2*S(3)].';
        k_p = k_p/sum(abs(k_p));
    end
    % Approach 1. Brute Force 
    n = 180;
    psi = linspace(0,pi/2,n);
    c = reshape(cos(2*psi),[1,1,n]);
    s = reshape(sin(2*psi),[1,1,n]);
    R = reshape([ones(1,1,n), zeros(1,1,n), zeros(1,1,n); zeros(1,1,n), c, -s; zeros(1,1,n), s, c],[3,3*n]).';
    temp = R*k_p;
    [~, ind] = min(abs(temp(3:3:end)));
    psi_br = psi(ind);
    % Approach 2. Analytical Solution
    alpha = acos(k_p(1));
    beta = acos(k_p(2)/sin(alpha));
    phi_2 = angle(k_p(2));
    phi_3 = angle(k_p(3));
    psi = mod((2*(cos(phi_2-phi_3)>=0)-1)*(2*beta - mod(2*beta,pi) ......
    	+ mod(atan(tan(2*beta)*abs(cos(phi_2-phi_3))), pi))/4, pi/2);
    fprintf('Brute force: %f, Analytical: %f \n', psi_br/pi*180, psi/pi*180)
    
    % Generate the pattern 
    R = [1, 0, 0; 0, cos(2*psi), -sin(2*psi); 0, sin(2*psi), cos(2*psi)];
    k_p = (R*k_p).*[1, cos(2*psi), sin(2*psi)].';  
    mu = [k_p(2) k_p(3)];
	sig = [0.001, 0; 0, 0.001];
	F = mvnpdf([x(:) y(:)],mu, sig);
	F = reshape(F,size(th));
	F = ind2rgb(uint8((F/max(max(F)))*255),jet);
    
	figure
	surf(x,y,z,F,'EdgeColor','none')
	%surf(x,y,z,'FaceColor','w')
	view(45,30);
	xlabel('$|k''_{p3}| \sin (2 \psi_m)$','Interpreter', 'latex')
	ylabel('$|k''_{p2}| \cos (2 \psi_m)$','Interpreter', 'latex')
	zlabel('$|k''_{p1}|$','Interpreter', 'latex')
	plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph_Dipo')
end


