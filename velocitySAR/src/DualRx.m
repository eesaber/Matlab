function [v_xt, v_yt, a_xt, y_0t] = DualRx(vx,vy,ax,ay)
% This function generate the SAR signal regarding to the input parameter.
% Usage: Gen_signal(vx, vy, ax, ay), 'vx' is range velocity, 'vy' is azimuth
% velocity, 'ax' is range accelaration and 'ay' is azimuth accelration. If no input parameters, the velocity in both direction are set
% to zero
% Noth that the parameters of SAR need to be modified in the function.

    %% Parameters - C Band airborne SAR parameters
    % Target

    x_0 = [2199.6, 2196, 2185, 2166]; % Squint angle: 0, 2.3, 4.5, 7.0
    y_0 = [0, 116, 249, 382];
	
    x_0 = x_0(4);
    y_0 = y_0(4);
    d = 200; % Length of the target area 
	v_x = vx; a_x = ax; % rangecl -16
	v_y = vy; a_y = ay; % azimuth 

%	v_x = -15; a_x = 3; % rangecl -16
%	v_y = -15; a_y = 0; % azimuth 

    % Platform
    h = 2200;
	R_0 = sqrt(x_0^2 + y_0^2 + h^2);
    d_a = 1;     
    v_p = 90; 
	% f_0 need to be large than f_1
    f_0 = 9.6e9; c = 299792458 ; lambda = c/f_0 ;
	f_1 = 8.9e9; lambda1 = c/f_1;
    dur = 2 ; %
    PRF = 4000; %%2.35e3
    K_r = 1e13 ;
    T_p = 5e-6; % Pulse width
    B =  K_r*T_p;
	fsamp = 1/16/B;
	%save('para.mat')
	%%
    % Signal 
	aa = 0;
    eta = (aa - dur/2): 1/PRF: (aa + dur/2);  % Slow time -6
	%dd = 20;
	%eta = (aa - dd): 1/10: (aa + dd);  % Slow time -6
    upran = sqrt( (x_0 + d)^2 + h^2) ; downran =  sqrt( (x_0 - d)^2 + h^2); %Set upper limit and down limit 
    tau =  2*(downran)/c: fsamp : 2*(upran)/c + T_p  ;  % fast time space
    R1 = sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation 
    R2 = sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (d_a + y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
	%R3 =  sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (2*d_a + y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ;
	R2 = R2 + R1;
	%R3 = R3 + R1;
	R1 = R1 * 2;
	
	temp = repmat(tau,length(eta),1);
	w_r1 = ( temp >= repmat(R1.', 1, length(tau)) / c  & temp <= (repmat(R1.', 1, length(tau))/c + T_p)) ;
	w_r2 = ( temp >= repmat(R2.', 1, length(tau)) / c  & temp <= (repmat(R2.', 1, length(tau))/c + T_p)) ;
	%w_r3 = ( temp >= repmat(R3.', 1, length(tau)) / c  & temp <= (repmat(R3.', 1, length(tau))/c + T_p)) ;
	s1 = exp(-1i*2*pi* repmat(R1.', 1, length(tau)) / lambda + 1i*pi*K_r*((temp - repmat(R1.'/c, 1, length(tau)) ).^2)) .* w_r1;
	s2 = exp(-1i*2*pi* repmat(R2.', 1, length(tau)) / lambda + 1i*pi*K_r*((temp - repmat(R2.'/c, 1, length(tau)) ).^2)) .* w_r2;
	s3 = exp(-1i*2*pi* repmat(R1.', 1, length(tau)) / lambda1 + 1i*pi*K_r*((temp - repmat(R1.'/c, 1, length(tau)) ).^2)) .* w_r1;
	s4 = exp(-1i*2*pi* repmat(R2.', 1, length(tau)) / lambda1 + 1i*pi*K_r*((temp - repmat(R2.'/c, 1, length(tau)) ).^2)) .* w_r2;
	add_noise = 0;
	if add_noise 
		w_noise = wgn(length(eta), length(tau),0,'complex');
		s1 = s1 + w_noise;
		s2 = s2 + w_noise;
		s3 = s3 + w_noise;
		s4 = s4 + w_noise;
		clear w_noise
	end
	clear w_r1 w_r2 R_2 R_1
	
	%{
	purinto(s1) 
	xlabel('$\tau / \Delta f_\tau $', 'Interpreter', 'latex')
	ylabel('$\eta / \Delta \eta$', 'Interpreter', 'latex') 
	%caxis([0 1])
	%export_fig 1.jpg
	
	figure
	plot(R2 - R1,'k','Linewidth',1.5)
	hold
	plot(d_a*(y_0/R_0 - (v_p-v_y)/R_0 *eta ),'k-.','Linewidth',1.5)
	xlim([0 8000])
	xlabel('$\eta / \Delta \eta $', 'Interpreter', 'latex')
	ylabel('$\Delta R(\eta)$', 'Interpreter', 'latex') 
	%plot_para('Maximize',true,'Filename','1')
	%}
	%% Range compression 
	ref_time = -T_p : fsamp : 0 ;
	h_m = repmat(exp(1j * 4 * pi * f_0 * R_0 / c) *  exp(- 1j * pi * K_r * ref_time.^2),length(eta),1 ); 
	%f_tau = linspace(0, 2*B, length(tau));
	
	%Fh_m = repmat(exp(-1j * pi * f_tau.^2 / K_r) .* (f_tau <= B & f_tau >= 0 ), length(eta), 1);
	%tau_nt2 = nextpow2(length(tau)) ;  % Make total number to power of 2
	tau_nt2 = length(tau) ;
    Fh_m = fft( h_m, tau_nt2, 2); 
    Fs1_rc = fft(s1, [], 2) .* Fh_m; % Range compression
    %s1_rc = ifft( Fs1_rc,[], 2);
    Fs1_rc = fftshift(Fs1_rc,2);
    
    Fs2_rc = fftshift(fft(s2, [], 2) .* Fh_m, 2); % Range compression
	Fs3_rc = fftshift(fft(s3, [], 2) .* Fh_m, 2); % Range compression
	Fs4_rc = fftshift(fft(s4, [], 2) .* Fh_m, 2); % Range compression
	clear s1 s2 s3 s4 h_m Fh_m
	
	%{
	fprintf('%f\n',2*(v_x*x_0 - y_0*(v_p - v_y))/R_0/lambda)
	purinto(fftshift(fft(fftshift(Fs1_rc,2)),1))
	xlabel('$\tau / \Delta \tau $', 'Interpreter', 'latex')
	ylabel('$f_\eta / \Delta f_\eta$', 'Interpreter', 'latex')
	
	figure
	plot(unwrap(angle(Fs1_rc(:,48))),'r','linewidth',1)
	figure
	plot(abs(Fs1_rc(:,48)))
	close all
	
	Fh_m = fft(h_m,length(tau), 2); 
    Fs1_rc = fft(s1,[], 2) .* Fh_m; % Range compression
    s1_rc = ifft( Fs1_rc,[], 2);
    Fs1_rc = fftshift(Fs1_rc,2);
    
    Fs2_rc = fftshift(fft(s2, [], 2) .* Fh_m, 2); % Range compression
	Fs3_rc = fftshift(fft(s3, [], 2) .* Fh_m, 2); % Range compression
	Fs4_rc = fftshift(fft(s4, [], 2) .* Fh_m, 2); % Range compression
	%}
	%{
	purinto(s1_rc)
	xlabel('$\tau / \Delta \tau $', 'Interpreter', 'latex')
	ylabel('$\eta / \Delta \eta$', 'Interpreter', 'latex')
	export_fig src.jpg
	
	purinto(fftshift(fft(fftshift(Fs1_rc,2)),1))
	xlabel('$\tau / \Delta \tau $', 'Interpreter', 'latex')
	ylabel('$f_\eta / \Delta f_\eta$', 'Interpreter', 'latex')
	export_fig src1.jpg
	
	%clear Fh_m s1 s2
	%h_de = exp(j*4*pi(f_tau + f_0)/c)
	
	figure, plot((atan2(imag(Fs1_rc(:,48)), real(Fs1_rc(:,48)))),'k')
	xlim([0 PRF*dur])
	xlabel('$\tau / \Delta \tau $', 'Interpreter', 'latex')
	ylabel('Phase', 'Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','1')
		
	purinto(Fs1_rc)
	xlabel('$f_\tau / \Delta f_\tau $', 'Interpreter', 'latex')
	ylabel('$\eta / \Delta \eta$', 'Interpreter', 'latex')
	%export_fig 2.jpg
	%}

	%% RCMC for range curvture 
	%-----------------------%
	%  Key Stone transform  %
	%-----------------------%
	
	%Interpolation 
	do_key = 1;
	if do_key
		
		f_tau = linspace(-1/fsamp/2, 1/fsamp/2 , tau_nt2) ;
		%f_tau = linspace(0, 1/fsamp/2 , 2^tau_nt2 ) ;
		t = eta.' * sqrt((f_0 + f_tau)/f_0);
		t = [t(1,:); t(end,:)];
		up_bound = max(max(t)) ;
		low_bound = min(min(t));
		Fs1_1 = zeros(floor( PRF * up_bound) - floor(PRF * low_bound) + 1, tau_nt2);
		Fs2_1 = Fs1_1;
		Fs3_1 = Fs1_1;
		Fs4_1 = Fs1_1;
		for m = 1 : tau_nt2
			
			up_ = floor(PRF * t(1,m)) - floor(PRF * low_bound) + 1;
			lo_ = floor(PRF * t(2,m)) - floor(PRF * low_bound) + 1;
			
			Fs1_1(up_:lo_, m) = ...
				interp1(eta, abs(Fs1_rc(:,m)).', linspace(eta(1),eta(end), lo_ - up_ + 1), 'spline', 0).';
			Fs1_1(up_:lo_, m) = Fs1_1(up_:lo_, m) .* exp(1j* ...
				interp1(eta, unwrap(angle(Fs1_rc(:,m))).', linspace(eta(1),eta(end), lo_  - up_ + 1), 'spline', 0).');
			Fs2_1(up_:lo_, m) = ...
				interp1(eta, abs(Fs2_rc(:,m)).', linspace(eta(1),eta(end), lo_ - up_ + 1), 'spline', 0).';
			Fs2_1(up_:lo_, m) = Fs2_1(up_:lo_, m) .* exp(1j* ...
				interp1(eta, unwrap(angle(Fs2_rc(:,m))).', linspace(eta(1),eta(end), lo_ - up_ + 1), 'spline', 0).');
			Fs3_1(up_:lo_, m) = ...
				interp1(eta, abs(Fs3_rc(:,m)).', linspace(eta(1),eta(end), lo_ - up_ + 1), 'spline', 0).';
			Fs3_1(up_:lo_, m) = Fs3_1(up_:lo_, m) .* exp(1j* ...
				interp1(eta, unwrap(angle(Fs3_rc(:,m))).', linspace(eta(1),eta(end), lo_ - up_ + 1), 'spline', 0).');
			Fs4_1(up_:lo_, m) = ...
				interp1(eta, abs(Fs4_rc(:,m)).', linspace(eta(1),eta(end), lo_ - up_ + 1), 'spline', 0).';
			Fs4_1(up_:lo_, m) = Fs4_1(up_:lo_, m) .* exp(1j* ...
				interp1(eta, unwrap(angle(Fs4_rc(:,m))).', linspace(eta(1),eta(end), lo_ - up_ + 1), 'spline', 0).');
		end
		
		s1_1 = ifft( ifftshift(Fs1_1, 2), [],2);
%{
		s2_1 = ifft( ifftshift(Fs2_1, 2), [],2);
		s3_1 = ifft( ifftshift(Fs3_1, 2), [],2);
		s4_1 = ifft( ifftshift(Fs4_1, 2), [],2);
	
		purinto(Fs1_1)
		xlabel('$f_\tau / \Delta f_\tau$', 'Interpreter', 'latex')
		ylabel('$t/\Delta t$', 'Interpreter', 'latex')
		caxis([0 160])
		%export_fig 1.jpg	

		figure
		plot(unwrap(angle(Fs1_1(:,48))),'k','linewidth',1)
		hold on 
		plot(unwrap(angle(Fs1_rc(:,48))),'r','linewidth',1)
		%xlim([0 length(Fs1_1(:,48))])
		xlabel('$t/\Delta t$', 'Interpreter', 'latex')
		ylabel('phase', 'Interpreter', 'latex')
		%plot_para('Maximize',true,'Filename','1')

		figure
		plot(diff(unwrap(angle(Fs1_1(:,48))),1),'k','Linewidth',2)
		xlabel('$t/\Delta t$', 'Interpreter', 'latex')
		ylabel('$\Delta$phase', 'Interpreter', 'latex')
		%xlim([280 520])
		%plot_para(1,'2') 

		purinto(s1_1)
		%title(['$(v_x,v_y)= ($' int2str(v_x) ',' int2str(v_y) '$)$'],'Interpreter', 'latex')
		xlabel('$\tau/ \Delta \tau$','Interpreter', 'latex')
		ylabel('$ t / \Delta t $','Interpreter', 'latex')
		%export_fig(['(' int2str(v_x) '_' int2str(v_y) ')_rc.jpg'])
		%}
	end	
	
	%------------------------%
	% Range Curvature filter %
	%------------------------%
    
	if ~do_key
		f_tau = [linspace(0, 1/fsamp/2, 2^(tau_nt2-1) ) linspace(0, 1/fsamp/2, 2^(tau_nt2-1))]; 
			
		H_c =  exp(  1j* (2* pi /c )* ( (v_y-v_p)^2 / R_0) *eta'.^2 * f_tau); %range curveture filter
		%Multiply H_c to Fs_r
		Fs1_1 = Fs1_rc .* H_c ;
		Fs2_1 = Fs2_rc .* H_c ;
		Fs3_1 = Fs3_rc .* H_c ;
		Fs4_1 = Fs4_rc .* H_c ;
		s1_1 = ifft( ifftshift(Fs1_1, 2), [],2);
		%s2_1 = ifft( ifftshift(Fs2_1, 2), [],2);
		%s3_1 = ifft( ifftshift(Fs3_1, 2), [],2);
		%s4_1 = ifft( ifftshift(Fs4_1, 2), [],2);
		%purinto(s1_1) 
		%xlabel('$\tau /\Delta \tau$', 'Interpreter', 'latex')
		%ylabel('$t/\Delta t$', 'Interpreter', 'latex')
	end 
	
	clear Fs1_rc Fs2_rc Fs3_rc Fs4_rc

	%% Radon transform to estimate the slop
	% real v_r  and the Doppler frequency
	% v_r_real_ = v_x * (h/R_0/sqrt(1-y_0^2/R_0^2) );
	
	do_radon = 1;
	if do_radon 
		if do_key 
			theta = 175:0.005:181;
		else
			theta = 170:0.05:181;
		end
		parral =0;
		if parral 
			[R,xp] = radon(gpuArray(abs(s1_1)),theta); % Use GPU to compute
			R = gather(R); xp = gather(xp);
		else
			[R,~] = radon(abs(s1_1(:, 4700:5000)),theta);
		%{
		figure
		imagesc(theta,xp/PRF,R)
		colormap('jet'),colorbar
		xlabel('$\theta$ (degrees)','Interpreter', 'latex')
		ylabel('$t''$','Interpreter', 'latex')
		%xlim([178 182]), ylim([0.32 0.4])
		plot_para('Maximize',true,'Filename','radon')
		%}
		end
		[~,ind_c] = max(max(R));
		ell = -1/(tand(theta(ind_c)+90)*1/fsamp/PRF);
	end
	clear s1_1
	%% RCMC for range walk 
	if do_key 
		if do_radon 
			Fh_rw = exp(1j * 2 * pi  * ell * linspace(t(1,end), t(end,end), length(Fs1_1(:,1))).' * f_tau );    
		else
			Fh_rw = exp(1j * 4 * pi / c * (v_x * x_0 - y_0*(v_p - v_y))/R_0 * ...
			linspace(t(1,end), t(end,end), length(Fs1_1(:,1))).' * f_tau);    
		end
	else 
		if do_radon 
			Fh_rw = exp(j * 2 * pi  * ell * eta.' * f_tau );    
		else
			Fh_rw = exp(j * 4 * pi / c * (v_x * x_0 - y_0*(v_p - v_y))/R_0 * eta.' * f_tau); 
		end
	end
    Fs1_2 = Fs1_1 .* Fh_rw ;
	s1_2 = ifft(ifftshift(Fs1_2, 2).').';
	clear Fs1_2 Fs1_1
	Fs2_2 = Fs2_1 .* Fh_rw ; 
	s2_2 = ifft(ifftshift(Fs2_2, 2).').';
	clear Fs2_2 Fs2_1
	Fs3_2 = Fs3_1 .* Fh_rw ;
	s3_2 = ifft(ifftshift(Fs3_2,2 ).').';
	clear Fs3_2 Fs3_1
	Fs4_2 = Fs4_1 .* Fh_rw ; 
    s4_2 = ifft(ifftshift(Fs4_2,2 ).').';
	clear Fs4_2 Fs4_1 Fh_rw
	%{
	purinto(s1_2)
	title(['$(v_x,v_y)= ($' int2str(v_x) ',' int2str(v_y) '$)$'],'Interpreter', 'latex')
	xlabel('$\tau/ \Delta \tau$','Interpreter', 'latex')
	ylabel('$ t / \Delta t $','Interpreter', 'latex')
	%export_fig(['(' int2str(v_x) '_' int2str(v_y) ')_rw.jpg'])
	%}
	%purinto(s1_2(:,1000:end))
	%purinto(s2_2(:,1000:end))
	%% Azimuth velocity Estimation
	[~, ind] = max(abs(s1_2),[],2);
	ind = floor(sum(ind)/length(s1_2(:,1)));
	R_0t = tau(ind- int32(T_p / fsamp))*c/2;
	method = string({'Corr filter', 'WVD', 'Phase Slope', 'GAF'}); 

	method = method(3);
	switch method 
		case 'Corr filter'  % Matched Filter bank
			v_ySpace = -20: 0.1: 20 ;
			qq = v_ySpace;
			t_temp = linspace(t(1,end), t(end, end), length(s1_2(:,ind)));
			for i = 1 : length(v_ySpace)
				if do_key
					s_sample = exp(-1j * (2* pi * f_0 / c)* ((v_ySpace(i) - v_p)^2 / R_0) *	t_temp.^2);
					s1_2_Fdc = s1_2(:,ind).' .* exp(1j * 2 * pi * f_0 * ell * t_temp);
				else
					s_sample = exp(-1j* (2* pi * f_0 / c)* ((v_ySpace(i) - v_p)^2 / R_0) *eta.^2);
					s1_2_Fdc = s1_2(:,ind).' .* exp(1j * 2 * pi * f_0 * ell * eta);
				end
				rr = max(abs(ifft(fft(s1_2_Fdc).*fft(conj(s_sample)))));
				qq(i) = rr;
				if  rr > temp
					temp = rr;
					v_yt = v_ySpace(i);
					v_y_crr = v_yt;
				end 
			end
			%figure 
			%plot(real(s1_2_Fdc))
			%figure
			%spectrogram(s1_2_Fdc,128,120,8196,PRF,'centered','yaxis')
		case 'WVD'  % WVD 
			if do_key 
				t_temp = linspace(t(1,end), t(end, end), length(s1_2(:,ind)));
				temp= spectrogram(s1_2(:,ind).* exp(1j * 2 * pi * f_0 * ell * t_temp).',128,120,8196,PRF,'centered','yaxis');
				x = [t(1,end), t(end,end)];
				t_len_ = length(temp(1,:));
				delta_t =  (x(2) - x(1))/2 ;
				[~ , upf] = max(temp(:,int16(t_len_/4)));
				[~ , dwf] = max(temp(:,int16(3*t_len_/4)));
				t_K_a = (dwf - upf) * PRF / 8196 / delta_t ;
			else
				temp= spectrogram(s1_2(:,ind).* exp(1j * 2 * pi * f_0 * ell * eta).',128,120,8196,PRF,'centered','yaxis');
				x = [(aa - dur/2) ,(aa + dur/2)];
				t_len_ = length(temp(1,:));
				[~ , upf] = max(temp(:,int16(t_len_/4)) );
				[~ , dwf] = max(temp(:,int16(3*t_len_/4)) );
				t_K_a = (dwf - upf) * PRF / 8196 / dur * 2;
			end
			v_yt = v_p - sqrt( - t_K_a * c * R_0 / f_0 / 2);
				
		case 'Phase Slope' % Search the parameters by phase slope
			temp = s1_2(:,ind).*conj(s2_2(:,ind));
			b_1 = lambda*angle(temp(int32(length(temp)/2)));
			temp = s3_2(:,ind).*conj(s4_2(:,ind));
			b_2 = lambda1*angle(temp(int32(length(temp)/2)));
			temp = s1_2(:,ind).*conj(s2_2(:,ind)) ./ (s3_2(:,ind).*conj(s4_2(:,ind)));
			b_w = c/(f_0 -f_1) * angle(temp( int32( length(temp)/2)));
			temp = s1_2(:,ind).*conj(s2_2(:,ind)) .* (s3_2(:,ind).*conj(s4_2(:,ind)));
			b_n = c/(f_0 +f_1) * angle(temp( int32( length(temp)/2)));
			
			temp = ambi_re([b_n, b_1, b_2, b_w], c./[f_0+f_1, f_0, f_1, f_0 - f_1], [d_a, R_0, v_p, v_y, y_0]);
			y_0t =  temp / (2*pi*d_a/R_0t);

			if do_key  
				t2fit = linspace(t(1,end), t(end,end), length(s1_2(:,1))) ;
			else
				t2fit = eta ;
			end
			[fit1,~,~] = fit(t2fit.', unwrap(angle(s1_2(:,ind).*conj(s2_2(:,ind)))), 'poly1');
			v_yt = v_p + fit1.p1 * lambda / 2 / pi *R_0 / d_a;
			
			
			%{
			figure
			plot(unwrap(angle(s1_2(:,ind).*conj(s2_2(:,ind)))))
			hold on 
			plot(-2*pi/lambda*d_a*( y_0/R_0 + (v_p - v_y)/R_0 * t2fit)+ 2*pi)	
			2*pi/lambda*d_a*( (v_p - v_y)/R_0 + y_0*(v_x*x_0 - y_0*(v_p - v_y))/R_0^3);
		
			figure
			plot(fit1)
			hold on 
			plot(t2fit,unwrap(angle(s1_2(:,ind).*conj(s2_2(:,ind)))))
			xlabel('$t$', 'Interpreter', 'latex')
			ylabel('Phase', 'Interpreter', 'latex')
			plot_para('Maximize',true)
			%}
		case 'GAF'
			temp = abs(GAF(s1_2(:, ind), 2, 2, 2));
			[~,ind_f] = max(temp);
			f_t = linspace(-PRF/2, PRF/2, length(temp));
			if do_key
				c_2 = f_t(ind_f) / abs(t(1,end) - t(end,end));
			else 
				c_2 = f_t(ind_f) / dur  ;
			end
			v_yt = v_p - sqrt(-c_2 * R_0 * lambda - a_x * x_0);
	end
	
	%fprintf('Method: %s, Actual: %f, Estimate: %f\n', method, v_y, v_yt)
	%% Range velocity Estimation
	v_rt = ell * c ;
	f_dc = 2*(v_x*x_0 - y_0*(v_p - v_y))/R_0t/lambda;
	
	x_0t = sqrt(R_0t^2 - h^2 - y_0t^2);
	v_xt = (ell * c * R_0t + y_0t*(v_p - v_y)) / x_0t;
	%fprintf('Actual v_x: %f, Estimate v_x: %f\n',v_x, v_xt)
	
	
	coeff_ = zeros(2,2);
    b_ = zeros(2,1);
    coeff_(1,1) = -x_0t*y_0t/R_0t^3 ;
    coeff_(1,2) = 1/R_0t - y_0t^2/R_0t^3;
    b_(1) = v_p/R_0t - y_0t^2*v_p/R_0t^3 + (c/2/pi/f_0/d_a)*fit1.p1;
	coeff_(2,1) = 1;
    coeff_(2,2) = y_0t/x_0t;
    b_(2) = (c*R_0t*ell + y_0t*v_p)/x_0t;
	A = coeff_ \ b_;
    v_xt = A(1);
	v_yt = A(2);
	
	%% Range acceleration Estimation
	f_eta = [linspace(0,PRF/2,length(s1_2)/2) linspace(-PRF/2,0-1/PRF,length(s1_2)/2)];
	a_xSpace = -10: 0.025: 10 ;
		temp = -1;
		for i = 1 : length(a_xSpace)		
			if do_key
				K_a = -2 * f_0/c * ( (v_yt-v_p)^2 + a_xSpace(i)*x_0t)/ R_0t;
				H_ac = exp(j*pi*v_rt/lambda*f_eta + j*pi*f_eta.^2/K_a).';
			else
				fprintf('If not keystone it is impossible to correct range curvature\n')
			end 
			rr = max(abs(ifft(fft(s1_2(:,ind)) .* H_ac)));
			if  rr > temp
				temp = rr;
				a_xt = a_xSpace(i);
			end
		end
	%fprintf('Actual a_x: %f, Estimate a_x: %f\n',a_x, a_xt)
	
	%% Azimuth compression
	%{
	K_a = -2 * f_0/c * ( (v_yt-v_p)^2 + a_xt*x_0 )/ R_0t;
	H_ac = exp(1j*pi*v_rt/lambda*f_eta + 1j*pi*f_eta.^2/K_a).';
	Fs_a = fft(s1_2) .* repmat(H_ac,1,length(s1_2(1,:)));
	s_a = ifft(Fs_a,[],1);
	purinto(s_a)
	%}
end
