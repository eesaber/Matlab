function Vis_Assem(k_p, sigma, varargin)
%VIS_ASSEM   Visualization of assembly matrix.
% k_p is a two-dimensional matrix with size 3xN. sigma is a three-dimensional matrix 
% with size 3x3xN

	if isunix
		cd /home/akb/Code/Matlab
	else 
		cd D:\Code\Simu\PolSAR
	end

	% Parse input parameter
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'logical'},{});
	validationFcn_2_ = @(x) validateattributes(x,{'logical'},{});
    validationFcn_3_ = @(x) validateattributes(x,{'logical'},{});
	addParameter(parse_,'Label_On',false, validationFcn_1_);
    addParameter(parse_,'SubPlot',false ,validationFcn_2_);    
	addParameter(parse_,'ThreeD',false, validationFcn_3_);
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
	n = 91;
	psi = linspace(0,90,n);
	c = reshape(cosd(2*psi),[1,1,n]);
	s = reshape(sind(2*psi),[1,1,n]);
	R = reshape([ones(1,1,n), zeros(1,1,n), zeros(1,1,n); zeros(1,1,n), c, -s; zeros(1,1,n), s, c],[3,3*n]).';
	
	[x_plain, y_plain] = meshgrid(linspace(-1,1,SP_NUM),linspace(-1,1,SP_NUM));

	alphabets = char(97:96+num_k_p).'; % Label for subplot, start from a to ...
	subplot_label = strcat({'('}, alphabets, {')'});
	
	for qq = 1 : num_k_p
        k_p(:,qq) = k_p(:,qq)/norm(k_p(:,qq),2); % Normalize k_p
		% Caulculate the deorientation angle
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
        %fprintf('Orientation angle: %f \n', psi_ana)
	
		%% Generate the pattern
		k_pp = abs(k_p(:,qq)).*[1, cosd(2*psi_ana), sind(2*psi_ana)].';
		mu = [k_pp(2) k_pp(3)]; % Let k_p(2) and k_p(3) be x-axis and y-axis respectively.  		
		F = reshape(mvnpdf([x(:) y(:)],mu, sigma(2:3,2:3,qq)), size(th));
		F_plain = reshape(mvnpdf([x_plain(:) y_plain(:)],mu, sigma(2:3,2:3,qq)), size(x_plain));
		F = ind2rgb(uint8((F/max(max(F)))*255),jet);
		F_plain = F_plain/max(max(F_plain));
		%% plot 
		if parse_.Results.SubPlot
			% We can plot a circular plot by using SURF
			% and view topdown
			subplot(2, 4, qq)
			surf(x,y,z,F,'EdgeColor','none')
			set(gca,'View',[0,90],'TickLength',[0,0],'YColor','none',......
                'XTickLabelMode','manual','GridColor','none')
            xlabel(subplot_label(qq),'Interpreter', 'latex','Fontsize',28)
            pbaspect([1 1 1])
		else
			% Plot each \bar{A} in different figure
			figure
			if parse_.Results.ThreeD
				surf(x,y,z,F,'EdgeColor','none')
				view(135,30);
				xlabel('$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$','Interpreter', 'latex')
				ylabel('$|k''_{p3}| \sin (2 \psi_m)/ \|\bar{k}_p\|_2 $','Interpreter', 'latex')
				zlabel('$|k''_{p1}| / \|\bar{k}_p\|_2 $','Interpreter', 'latex')
			else
				imagesc(x_plain(1,:), y_plain(:,1),-F_plain/max(max(F_plain)))
				xlabel('$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$','Interpreter', 'latex')
				ylabel('$|k''_{p3}| \sin (2 \psi_m)/ \|\bar{k}_p\|_2 $','Interpreter', 'latex')
				grid on
				colormap gray
            end
        end
	plot_para('Maximize',true,'Filename','SimAtom')
	movefile SimAtom.jpg PolSAR/output
end