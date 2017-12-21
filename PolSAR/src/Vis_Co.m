function Vis_Co(k_p)
%VIS_CO   Visualization of coherent target sphere.
%	Vis_Co(k_p) is a function to visualize the scattering mechanism.

    % Generate hemisphere.
	SP_NUM = 30;
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
        fprintf('Vis_Co IS IN TESTING MODE!\n')
        psi = 0*(pi/180);
        R = [cos(psi), sin(psi); -sin(psi), cos(psi)]; % Rotate counter-clockwisely with psi.
        S = R.'*[1, 0; 0, 0]*R;
		if isreal(S)
			fprintf('Scattering matrix is:\n %f %f\n %f %f \n',S)
		else
			fprintf('Scattering matrix is:\n %f%+fi, %f%+fi\n %f%+fi, %f%+fi \n',[real(S), imag(S)].')
		end
        k_p = 1/sqrt(2)*[S(1,1)+S(2,2), S(1,1)-S(2,2), 2*S(1,2)].';    
	end
	
	k_p = k_p/sqrt(sum(abs(k_p).^2));
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
    alpha = acos(abs(k_p(1)));
    beta = acos(abs(k_p(2))/sin(alpha));
	if ~isreal(beta)
		fprintf('beta is not real number, alpha = %f, beta = %f+%fi \n', alpha/pi*180, real(beta), imag(beta) )
		beta = real(beta);
	end
    phi_2 = angle(k_p(2));
    phi_3 = angle(k_p(3));
    psi = mod((2*(cos(phi_2-phi_3)>=0)-1)*(2*beta - mod(2*beta,pi) ......
    	+ mod(atan(tan(2*beta)*abs(cos(phi_2-phi_3))), pi))/4, pi/2);
    fprintf('Brute force: %f, Analytical: %f \n', psi_br/pi*180, psi/pi*180)
    %psi = psi_br;
    % Generate the pattern 
    R = [1, 0, 0; 0, cos(2*psi), -sin(2*psi); 0, sin(2*psi), cos(2*psi)];
    k_p = abs((R*k_p)).*[1, cos(2*psi), sin(2*psi)].';  
    mu = [k_p(2) k_p(3)]; % Let k_p(2) and k_p(3) be x-axis and y-axis respectively.
	sig = [0.001, 0; 0, 0.001];
	F = mvnpdf([x(:) y(:)],mu, sig);
	F = reshape(F,size(th));
	F = ind2rgb(uint8((F/max(max(F)))*255),jet);
    
    [x_plain, y_plain] = meshgrid(linspace(-1,1,SP_NUM),linspace(-1,1,SP_NUM));
    F_plain = mvnpdf([x_plain(:) y_plain(:)],mu, sig);
    F_plain = reshape(F_plain,size(x_plain));
    
    close all
    cd /home/akb/Code/Matlab
    
	figure
	%surf(x,y,z,F,'EdgeColor','none')
	surf(x,y,z,'FaceColor','w')
	view(135,30);
	xlabel('$|k''_{p2}| \cos (2 \psi_m)$','Interpreter', 'latex')
	ylabel('$|k''_{p3}| \sin (2 \psi_m)$','Interpreter', 'latex')
	zlabel('$|k''_{p1}|$','Interpreter', 'latex')
	plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph')
    figure
    imagesc(x_plain(1,:), y_plain(:,1),-F_plain/max(max(F_plain)))
    xlabel('$|k''_{p2}| \cos (2 \psi_m)$','Interpreter', 'latex')
	ylabel('$|k''_{p3}| \sin (2 \psi_m)$','Interpreter', 'latex')
    set(gca,'Ydir','normal','XGrid','on','YGrid','on')
    colormap gray
    plot_para('Ratio',[4 3 1],'Maximize', true,'Filename', 'VisSph_Dipo')
end