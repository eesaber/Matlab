function Vis_Co(k_p, varargin)
%VIS_CO   Visualization of coherent target sphere.
%	Vis_Co(k_p) is a function to visualize the scattering mechanism.
	
	% Parse input parameter
	
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'char'},{}); 
	validationFcn_2_ = @(x) validateattributes(x,{'char'},{});
    validationFcn_3_ = @(x) validateattributes(x,{'logical'},{});
	validationFcn_4_ = @(x) validateattributes(x,{'char'},{});
    validationFcn_5_ = @(x) validateattributes(x,{'numeric'},{'size',[3, 3]});
	addParameter(parse_,'xlabel','$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$', validationFcn_1_);
	addParameter(parse_,'ylabel','$|k''_{p3}| \sin (2 \psi_m) / \|\bar{k}_p\|_2$', validationFcn_2_);
    addParameter(parse_,'Subplot',0 ,validationFcn_3_);
	addParameter(parse_,'Filename','1', validationFcn_4_);
    addParameter(parse_,'Sigma',0.1*eye(3), validationFcn_5_);
	parse(parse_,varargin{:})
	
    % Generate hemisphere.
	SP_NUM = 500;
	thetavec = linspace(0,pi/2,SP_NUM/2);
	phivec = linspace(0,2*pi,2*SP_NUM);
	[th, ph] = meshgrid(thetavec,phivec);
	R = ones(size(th)); % should be your R(theta,phi) surface in general
	x = R.*sin(th).*cos(ph);
	y = R.*sin(th).*sin(ph);
	z = R.*cos(th);
    % Testing
    test = 0;
    if test
        fprintf('Vis_Co IS IN TESTING MODE!\n')
        psi = 0*(pi/180);
        R = [cos(psi), sin(psi); -sin(psi), cos(psi)]; % Rotate counter-clockwisely with psi.
        %S = R.'*[1, -1j; -1j, -1]*R;
		S = R.'*[1, 0; 0, 0]*R;
		if isreal(S)
			fprintf('Scattering matrix is:\n %f %f\n %f %f \n',S)
		else
			fprintf('Scattering matrix is:\n %f%+fi, %f%+fi\n %f%+fi, %f%+fi \n',[real(S), imag(S)].')
		end
        k_p = 1/sqrt(2)*[S(1,1)+S(2,2), S(1,1)-S(2,2), 2*S(1,2)].';    
    end
    % Normalize k_p
    %k_p = k_p/norm(k_p,2);
    
    % Approach 1. Brute Force 
    n = 90;
    psi = linspace(0,pi/2-0.001,n);
    c = reshape(cos(2*psi),[1,1,n]);
    s = reshape(sin(2*psi),[1,1,n]);
    R = reshape([ones(1,1,n), zeros(1,1,n), zeros(1,1,n); zeros(1,1,n), c, -s; zeros(1,1,n), s, c],[3,3*n]).';
    temp = R*k_p;
    [sup_, ind] = min(abs(temp(3:3:end)));
    psi_br = psi(ind);
    % Approach 2. Analytical Solution
    alpha = acos(abs(k_p(1))/sqrt(2));
    beta = acos(abs(k_p(2)) / (sqrt(2)*norm(k_p,2)*sin(alpha)));
	if ~isreal(beta)
		fprintf('beta is not real number, alpha = %f, beta = %f+%fi \n', alpha/pi*180, real(beta), imag(beta) )
		beta = real(beta);
	end
    phi_2 = angle(k_p(2));
    phi_3 = angle(k_p(3));
    psi_ana = mod((2*(cos(phi_2-phi_3)>=0)-1)/4*(2*beta - mod(2*beta,pi) ......
    	+ mod(atan(tan(2*beta)*abs(cos(phi_2-phi_3))), pi)), pi/2);
    fprintf('Brute force: %f, Analytical: %f \n', psi_br/pi*180, psi_ana/pi*180)    
	
	% Generate the pattern 
	% Check if there exist more than one min value
	F = zeros(size(th));
	[x_plain, y_plain] = meshgrid(linspace(-1,1,SP_NUM),linspace(-1,1,SP_NUM));
	F_plain = zeros(size(x_plain));
    
    
    R = [1, 0, 0; 0, cos(2*psi_ana), -sin(2*psi_ana); 0, sin(2*psi_ana), cos(2*psi_ana)];
    %k_pp = abs((R*k_p)).*[1, cos(2*psi_ana), sin(2*psi_ana)].';  
    k_pp = abs(k_p).*[1, cos(2*psi_ana), sin(2*psi_ana)].';  
    
    k_pp = k_pp/norm(k_pp, 2);
    mu = [k_pp(2) k_pp(3)]; % Let k_p(2) and k_p(3) be x-axis and y-axis respectively.  		
    sig = parse_.Results.Sigma(2:3, 2:3);
    temp_F = mvnpdf([x(:) y(:)],mu, sig);
    F = F + reshape(temp_F,size(th));
    temp_P = mvnpdf([x_plain(:) y_plain(:)],mu, sig);
    F_plain = F_plain + reshape(temp_P,size(x_plain));

    
    %{
	allkp_3 = temp(3:3:end);
    min_low = 2.2251e-10;
	psi_p = psi(abs(abs(allkp_3) - sup_) < min_low);
	for n = 1 : numel(psi_p)
		R = [1, 0, 0; 0, cos(2*psi_p(n)), -sin(2*psi_p(n)); 0, sin(2*psi_p(n)), cos(2*psi_p(n))];
		k_pp = abs((R*k_p)).*[1, cos(2*psi_p(n)), sin(2*psi_p(n))].';  
        k_pp = k_pp/norm(k_pp, 2);
		mu = [k_pp(2) k_pp(3)]; % Let k_p(2) and k_p(3) be x-axis and y-axis respectively.  		
        sig = parse_.Results.Sigma(2:3, 2:3);
		temp_F = mvnpdf([x(:) y(:)],mu, sig);
		F = F + reshape(temp_F,size(th));
		temp_P = mvnpdf([x_plain(:) y_plain(:)],mu, sig);
		F_plain = F_plain + reshape(temp_P,size(x_plain));
	end
    %}
	F = ind2rgb(uint8((F/max(max(F)))*255),jet);
	F_plain = F_plain/max(max(F_plain));
    %close all
	if isunix
		cd /home/akb/Code/Matlab
	end
    %{
	figure
	surf(x,y,z,F,'EdgeColor','none')
	%surf(x,y,z,'FaceColor','w')
	view(135,30);
	xlabel('$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$','Interpreter', 'latex')
	ylabel('$|k''_{p3}| \sin (2 \psi_m)/ \|\bar{k}_p\|_2 $','Interpreter', 'latex')
	zlabel('$|k''_{p1}| / \|\bar{k}_p\|_2 $','Interpreter', 'latex')
	%plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph','Fontsize',32)
    %}
    if parse_.Results.Subplot
        figure(1)
        ax_1 = axes;
        imagesc(ax_1, x_plain(1,:), y_plain(:,1),F_plain/max(max(F_plain)));
        ax_1.Visible = 'off';
        ax_1.XTick = [];
        ax_1.YTick = [];
        plot_para('Ratio',[1 1 1])

        ax_2 = axes;
        imagesc(ax_2, x_plain(1,:), y_plain(:,1),sqrt(x_plain.^2 + y_plain.^2)>1,.....
            'AlphaData', sqrt(x_plain.^2 + y_plain.^2)>1);
        linkaxes([ax_1,ax_2])
        ax_2.Visible = 'off';
        ax_2.XTick = [];
        ax_2.YTick = [];
        colormap(ax_1, 'jet')
        colormap(ax_2, 'gray')
        xlabel(ax_2, parse_.Results.xlabel,'Interpreter', 'latex')
        ylabel(ax_2, parse_.Results.ylabel,'Interpreter', 'latex')
        set(ax_1,'Ydir','normal','XGrid','on','YGrid','on','GridAlpha', .5, 'GridColor', 'w')
        set(ax_2,'Ydir','normal')
        axis off
        plot_para('Ratio',[1 1 1],'Maximize',true,'Filename',parse_.Results.Filename)
        close all
    else
        figure
            imagesc(x_plain(1,:), y_plain(:,1),-F_plain/max(max(F_plain)))
            xlabel('$|k''_{p2}| \cos (2 \psi_m) / \|\bar{k}_p\|_2$','Interpreter', 'latex')
            ylabel('$|k''_{p3}| \sin (2 \psi_m)/ \|\bar{k}_p\|_2 $','Interpreter', 'latex')
            grid on
            colormap gray
            %plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph','Fontsize',32)
    end
end