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
	SP_NUM = 700;
	thetavec = linspace(0,pi/2,SP_NUM/2);
	phivec = linspace(0,2*pi,2*SP_NUM);
	[th, ph] = meshgrid(thetavec,phivec);
	R = ones(size(th)); % should be your R(theta,phi) surface in general
	x = R.*sin(th).*cos(ph);
	y = R.*sin(th).*sin(ph);
	z = R.*cos(th);
    
    %% Parameters
	[~, num_k_p] = size(k_p);
	sig = 1e-4*eye(2);
	F = zeros(size(th));
	[x_plain, y_plain] = meshgrid(linspace(-1,1,SP_NUM),linspace(-1,1,SP_NUM));
	F_plain = zeros(size(x_plain));

	%% Plot all the k_p
    for qq = 1 : num_k_p
        % Normalize k_p
        k_p(:,qq) = k_p(:,qq)/norm(k_p(:,qq),2);
        % Judge if it is symmetry
		S = 1/sqrt(2)*[k_p(1,qq)+k_p(2,qq); sqrt(2)*k_p(3,qq); k_p(1,qq)-k_p(2,qq)];
		chi = 0.5*atand((conj(k_p(2,qq))*k_p(3,qq) + k_p(2,qq)*conj(k_p(3,qq)))/(abs(k_p(2,qq))^2 + abs(k_p(3,qq))^2));
		DS = dot(S, S_a)*S_a + dot(S, cosd(chi)*S_b+sind(chi)*S_c)*(cosd(chi)*S_b+sind(chi)*S_c);
		tau = acos(dot(S,DS)/norm(S)/norm(DS));
		symmetry = (tau <= 8/pi);
		if tau > 8/pi
			symmetry = 0;
			disp('asymmetry target')
		end
		% find the deorientation angle
		alpha = acosd(abs(k_p(1,qq)));
		beta = acosd(abs(k_p(2,qq))/sind(alpha));
		if ~isreal(beta)
			fprintf('beta is not real number, alpha = %f, beta = %f+%fi \n', alpha, real(beta), imag(beta) )
			beta = real(beta);
		end
		phi_2 = atan2d(imag(k_p(2,qq)), real(k_p(2,qq)));
		phi_3 = atan2d(imag(k_p(3,qq)), real(k_p(3,qq)));
		psi_ana = mod((2*(cosd(phi_2-phi_3)>=0)-1)/4*(2*beta - mod(2*beta, 180) ......
			+ mod(atand(tand(2*beta)*abs(cosd(phi_2-phi_3))), 180)), 90);
		S_hh = (k_p(1,qq)+k_p(2,qq))/sqrt(2);
        S_vv = (k_p(1,qq)-k_p(2,qq))/sqrt(2);
        S_hv = k_p(3,qq)/sqrt(2);
        a = atan2d(abs(S_vv),abs(S_hh));
        %b = 0.5*(atan2d(imag(S_vv),real(S_vv)-atan2d(imag(S_hh),real(S_hh));
        c = acosd(sqrt(2)*abs(S_hv)/norm([S_hh, sqrt(2)*S_hv, S_vv]));
        u = sind(c)*cosd(2*a);
        psi_ana = psi_ana+90*(u<0);
        %fprintf(Analytical: %f \n', psi_ana)
	
		mu = [k_pp(2) k_pp(3)]; % Let k_p(2) and k_p(3) be x-axis and y-axis respectively.
		% "mvnpdf()" IS THE BOTTLENECK when num_k_p is large.
		F = F + reshape(mvnpdf([x(:) y(:)],mu, sig), size(th));
		F_plain = F_plain + reshape(mvnpdf([x_plain(:) y_plain(:)],mu, sig), size(x_plain));
    end 
	F_color = ind2rgb(uint8((F/max(max(F)))*255),jet);
    F_black = -(F/max(max(F)))*255;
	F_plain = F_plain/max(max(F_plain));

    if parse_.Results.ThreeD
		figure
			surf(x,y,z,F_black,'EdgeColor','none')
            colormap gray
			view(135,30);
			xlabel('$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$','Interpreter', 'latex')
			ylabel('$|k''_{p3}| \sin (2 \psi_m)/ \|\bar{k}_p\|_2 $','Interpreter', 'latex')
			zlabel('$|k''_{p1}| / \|\bar{k}_p\|_2 $','Interpreter', 'latex')
            hold on 
            %% Generate mesh grid
            SP_NUM = 20;
            thetavec = linspace(0,pi/2,SP_NUM/2);
            phivec = linspace(0,2*pi,2*SP_NUM);
            [th, ph] = meshgrid(thetavec,phivec);
            R = ones(size(th)); % should be your R(theta,phi) surface in general
            x = R.*sin(th).*cos(ph);
            y = R.*sin(th).*sin(ph);
            z = R.*cos(th);			
            surf(x,y,z,'FaceColor','none','FaceAlpha',1);
            hold off
            %plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph','Fontsize',32)
	else 
        figure
            imagesc(x_plain(1,:), y_plain(:,1),-F_plain)
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