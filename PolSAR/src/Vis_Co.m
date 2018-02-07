function Vis_Co(k_p, varargin)
%VIS_CO   Visualization of coherent target sphere.
%	Vis_Co(k_p) is a function to visualize the scattering mechanism.
    if isunix
		cd /home/akb/Code/Matlab
    else 
		cd D:\Code\Simu\PolSAR
    end
	%% Parse input parameter
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'logical'},{});
	addParameter(parse_,'ThreeD',false, validationFcn_1_);
	parse(parse_,varargin{:})
	
    %% Generate hemisphere.
	SP_NUM = 500;
	thetavec = linspace(0,pi/2,SP_NUM/2);
	phivec = linspace(0,2*pi,2*SP_NUM);
	[th, ph] = meshgrid(thetavec,phivec);
	R = ones(size(th)); % should be your R(theta,phi) surface in general
	x = R.*sin(th).*cos(ph);
	y = R.*sin(th).*sin(ph);
	z = R.*cos(th);
    
    %%
	[~, num_k_p] = size(k_p);

	n = 722;
	psi = linspace(0,pi/2,n);
	c = reshape(cos(2*psi),[1,1,n]);
	s = reshape(sin(2*psi),[1,1,n]);
	R = reshape([ones(1,1,n), zeros(1,1,n), zeros(1,1,n); zeros(1,1,n), c, -s; zeros(1,1,n), s, c],[3,3*n]).';
	sig = 1e-4*eye(2);
	F = zeros(size(th));
	[x_plain, y_plain] = meshgrid(linspace(-1,1,SP_NUM),linspace(-1,1,SP_NUM));
	F_plain = zeros(size(x_plain));

    for qq = 1 : num_k_p
        % Normalize k_p
        k_p(:,qq) = k_p(:,qq)/norm(k_p(:,qq),2);
        
		% Approach 1. Brute Force 	
		temp = R*k_p(:,qq);
		[sup_, ind] = min(abs(temp(3:3:end)));
		psi_br = psi(ind);
		% Approach 2. Analytical Solution
        
		alpha = acos(abs(k_p(1,qq)));
		beta = acos(abs(k_p(2,qq))/sin(alpha));
		if ~isreal(beta)
			fprintf('beta is not real number, alpha = %f, beta = %f+%fi \n', alpha/pi*180, real(beta), imag(beta) )
			beta = real(beta);
		end
		phi_2 = angle(k_p(2,qq));
		phi_3 = angle(k_p(3,qq));
		psi_ana = mod((2*(cos(phi_2-phi_3)>=0)-1)/4*(2*beta - mod(2*beta,pi) ......
			+ mod(atan(tan(2*beta)*abs(cos(phi_2-phi_3))), pi)), pi/2);
		%fprintf('Brute force: %f, Analytical: %f \n', psi_br/pi*180, psi_ana/pi*180)    
        
        %fprintf('Brute force: %f \n', psi_br/pi*180)
        
		%% Generate the pattern
        
		R_t = [1, 0, 0; 0, cos(2*psi_ana), -sin(2*psi_ana); 0, sin(2*psi_ana), cos(2*psi_ana)];
		k_pp = abs(k_p(:,qq)).*[1, cos(2*psi_ana), sin(2*psi_ana)].' 
		mu = [k_pp(2) k_pp(3)]; % Let k_p(2) and k_p(3) be x-axis and y-axis respectively.  		
		F = F + reshape(mvnpdf([x(:) y(:)],mu, sig), size(th));
		F_plain = F_plain + reshape(mvnpdf([x_plain(:) y_plain(:)],mu, sig), size(x_plain));
        
        %{
		allkp_3 = temp(3:3:end);
		psi_p = psi(abs(abs(allkp_3) - sup_) < 2.2251e-20);
		for n = 1 : numel(psi_p)
			R_t = [1, 0, 0; 0, cos(2*psi_p(n)), -sin(2*psi_p(n)); 0, sin(2*psi_p(n)), cos(2*psi_p(n))];
            %R_t*k_p(:,qq)
			k_pp = abs(k_p(:,qq)).*[1, cos(2*psi_p(n)), sin(2*psi_p(n))].';  
			mu = [k_pp(2) k_pp(3)]; % Let k_p(2) and k_p(3) be x-axis and y-axis respectively.
			F = F + reshape(mvnpdf([x(:) y(:)],mu, sig), size(th));
			F_plain = F_plain + reshape(mvnpdf([x_plain(:) y_plain(:)],mu, sig), size(x_plain));
		end
        %}
    end 
	F = ind2rgb(uint8((F/max(max(F)))*255),jet);
	F_plain = F_plain/max(max(F_plain));

    if parse_.Results.ThreeD
		figure
			surf(x,y,z,F,'EdgeColor','none')
			view(135,30);
			xlabel('$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$','Interpreter', 'latex')
			ylabel('$|k''_{p3}| \sin (2 \psi_m)/ \|\bar{k}_p\|_2 $','Interpreter', 'latex')
			zlabel('$|k''_{p1}| / \|\bar{k}_p\|_2 $','Interpreter', 'latex')
			%plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph','Fontsize',32)
            colormap gray
	else 
        figure
            imagesc(x_plain(1,:), y_plain(:,1),-F_plain/max(max(F_plain)))
            xlabel('$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$','Interpreter', 'latex')
            ylabel('$|k''_{p3}| \sin (2 \psi_m)/ \|\bar{k}_p\|_2 $','Interpreter', 'latex')
            grid on
            colormap gray
            %plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph','Fontsize',32)
    end
	%% Display the sphere coordinate
	if 0
	figure
		surf(x,y,z,'FaceColor','w')
		view(135,30);
		xlabel('$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$','Interpreter', 'latex')
		ylabel('$|k''_{p3}| \sin (2 \psi_m)/ \|\bar{k}_p\|_2 $','Interpreter', 'latex')
		zlabel('$|k''_{p1}| / \|\bar{k}_p\|_2 $','Interpreter', 'latex')
		plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph','Fontsize',32)
	end
end