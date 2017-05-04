function [v_yt] = DualRx(vx,vy,ax,ay)
% This function generate the SAR signal regarding to the input parameter.
% Usage: Gen_signal(vx, vy, ax, ay), 'vx' is range velocity, 'vy' is azimuth
% velocity, 'ax' is range accelaration and 'ay' is azimuth accelration. If no input parameters, the velocity in both direction are set
% to zero
% Noth that the parameters of SAR need to be modified in the function.
%{
    delete parameter.mat
    if nargin == 0
        vx = 0; vy = 0; 
        ax = 0; ay = 0;
    end
%}
	close all 
    %% Parameters - C Band airborne SAR parameters
    % Target
    x_0 = 2160;
    y_0 = 400;
    d = 100; % Length of the target area 
    %v_x = vx; a_x = ax; % rangecl -16
    %v_y = vy; a_y = ay; % azimuth 
	v_x = -4; a_x = 0; % rangecl -16
    v_y = -4; a_y = 0; % azimuth 

    % Platform
    h = 2200;
	R_0 = sqrt(x_0^2 + y_0^2 + h^2);
	l_0 = sqrt( y_0^2 + h^2);
    d_a = 0.2;     
    v_p = 90; 
    f_0 = 9.6e9; c = 3e8 ; lambda = c/f_0 ;
    dur = 2 ; %
    PRF = 3000; %%2.35e3
    K_r = 5e14 ;
    T_p = 0.1e-6; % Pulse width
    B =  K_r*T_p;
	fsamp = 1/8/B;
	
    % Signal 
    aa = 0;
    eta = (aa - dur/2): 1/PRF: (aa + dur/2);  % Slow time -6
	
    upran = sqrt( (x_0 + d)^2 + h^2) ; downran =  sqrt( (x_0 - d)^2 + h^2); %Set upper limit and down limit 
    tau =  2*(downran)/c: fsamp : 2*(upran)/c + T_p  ;  % fast time space

    R1 = sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
    R2 = sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (d_a + y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
	
    s1 = zeros(length(eta), length(tau));
	s2 = s1;
	%zxc = zeros(1,length(tau));
	for k = 1 : length(eta)
        td1 = 2* R1(k)/ c ;
		td2 = (R1(k) + R2(k)) /c ;
		s1(k,:) = exp(-1i* 4* pi*R1(k)/ lambda + 1i* pi* K_r* (tau - td1).^2) .*(tau>=td1 & tau-td1<=T_p)  ;
		s2(k,:) = exp(-1i* 2* pi*(R1(k) + R2(k))/ lambda + 1i* pi* K_r* (tau - td2).^2) .*(tau>=td2 & tau-td2<=T_p)  ;
		%{
		zxc((tau>=td1 & tau-td1<=T_p)) = linspace(0, T_p, length(tau((tau>=td1 & tau-td1<=T_p))));
		s1(k,:) = exp(-1i* 4* pi*R1(k)/ lambda + 1i* pi* K_r* zxc.^2) .*(tau>=td1 & tau-td1<=T_p);
		zxc((tau>=td2 & tau-td2<=T_p)) = linspace(0, T_p, length(tau((tau>=td2 & tau-td2<=T_p))));
		s2(k,:) = exp(-1i* 2* pi*(R1(k) + R2(k))/ lambda + 1i* pi* K_r* zxc.^2) .*(tau>=td2 & tau-td2<=T_p) ;
		%}
	end
	%{
	purinto(s1)
	xlabel('$\tau / \Delta f_\tau $', 'Interpreter', 'latex')
	ylabel('$\eta / \Delta \eta$', 'Interpreter', 'latex')
	%caxis([0 1])
	%export_fig 1.jpg
	%}
	%purinto(s2)
	clear R1 R2
	
	%% Range compression 
	ref_time = -T_p : fsamp : 0 ;
	h_m = repmat(exp(j * 4 * pi * f_0 * R_0 / c) *  exp(- j * pi * K_r * ref_time.^2),length(eta),1 ) ; 

	tau_nt2 = nextpow2(length(tau)) ;  % Make total number to power of 2
    Fh_m = fft( h_m, 2^tau_nt2, 2); 
    Fs1_rc = fft(s1, 2^tau_nt2, 2) .* Fh_m; % Range compression
    s1_rc = ifft( Fs1_rc,[], 2);
    Fs2_rc = fft(s2, 2^tau_nt2, 2) .* Fh_m; % Range compression
	s2_rc = ifft( Fs2_rc,[], 2);
	%clear Fh_m s1 s2
	figure, plot((atan2(imag(Fs1_rc(:,48)), real(Fs1_rc(:,48)))),'k')
	xlim([0 PRF*dur])
	xlabel('$\tau / \Delta \tau $', 'Interpreter', 'latex')
	ylabel('Phase', 'Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','1')
	%{
	purinto(s1_rc)
	xlabel('$\tau / \Delta \tau $', 'Interpreter', 'latex')
	ylabel('$\eta / \Delta \eta$', 'Interpreter', 'latex')
	export_fig src.jpg
	
	purinto(Fs1_rc)
	xlabel('$f_\tau / \Delta f_\tau $', 'Interpreter', 'latex')
	ylabel('$\eta / \Delta \eta$', 'Interpreter', 'latex')
	export_fig 2.jpg
	
	purinto(s1_rc)
	xlabel('$\tau / \Delta \tau $', 'Interpreter', 'latex')
	ylabel('$\eta / \Delta \eta$', 'Interpreter', 'latex')
	export_fig src_f.jpg
	%}
	figure
	plot(-20:1:20,v_Est_Err_Vy_(17,:),'k','Linewidth',2)
	xlabel('$v_y$', 'Interpreter', 'latex')
	ylabel('$\tilde{v}_y - v_y $', 'Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','1')
	
	%% Key Stone transform - RCMC for range curveture 
	% Interpolation 
	
    f_tau = linspace(0, 1/fsamp , 2^tau_nt2 ); 
	t = eta.' * sqrt((f_0 + f_tau)/f_0);
	%Fs1_1 = zeros(2*floor(t(end,end)*PRF)+4*PRF+1,2^tau_nt2);
	Fs1_1 = zeros(2*floor(t(end,end)*PRF)+1,2^tau_nt2);
	Fs2_1 = zeros(2*floor(t(end,end)*PRF)+1,2^tau_nt2);
	add_ = floor(PRF*t(end,end))-PRF;
	
	for m = 1 : 2^tau_nt2
		[min_p,ind] = min(unwrap(angle(Fs1_rc(:,m))));
		if ~(ind == 1 || ind == length(Fs1_rc(:,m)))
			 %Fs1_rc(ind+1:end,m) = Fs1_rc(ind+1:end,m) .* exp(-1j*2*unwrap(angle(Fs1_rc(ind+1:end,m))));
			 Fs1_rc(1:ind,m) = Fs1_rc(1:ind,m) .* exp(-1j*2*unwrap(angle(Fs1_rc(1:ind,m))));
		end
		cou_ = floor(PRF*t(end,end))-floor(PRF*t(end,m));
		Fs1_1(:, m) = [zeros(cou_,1); ...
			interp1(eta, abs(Fs1_rc(:,m)).', linspace(-1,1, 2*add_+ PRF*dur - 2*cou_+1), 'spline', 0).'; zeros(cou_,1)];
		Fs1_1(:, m) = Fs1_1(:, m) .* exp(1j*[zeros(cou_,1); ...
			interp1(eta, unwrap(angle(Fs1_rc(:,m))).', linspace(-1,1, 2*add_+ PRF*dur - 2*cou_+1), 'spline', 0).'; zeros(cou_,1)]);
		[~,ind] = min(abs(unwrap(angle(Fs1_1(:,m)))-min_p));
		Fs1_1(ind+1:end,m) = Fs1_1(ind+1:end,m) .* exp(-j*2*unwrap(angle(Fs1_1(ind+1:end,m))));
		%Fs1_1(1:ind,m) = Fs1_1(1:ind,m) .* exp(-j*2*unwrap(angle(Fs1_1(1:ind,m))));
		%{%}
		Fs2_1(:, m) = [zeros(cou_,1); ...
			interp1(eta, abs(Fs2_rc(:,m)).', linspace(-1,1, 2*add_+ PRF*dur - 2*cou_+1), 'spline', 0).'; zeros(cou_,1)];
		Fs2_1(:, m) = Fs2_1(:,m) .* exp(1j * [zeros(cou_,1); ...
			interp1(eta, unwrap(angle(Fs2_rc(:,m))).', linspace(-1,1, 2*add_+ PRF*dur - 2*cou_+1), 'spline', 0).'; zeros(cou_,1)]);
	end

	s1_1 = ifft(Fs1_1, [],2);
	s2_1 = ifft(Fs2_1, [],2);
	
	figure
	plot(unwrap(angle(Fs1_1(:,48))),'k','linewidth',1)
	hold on 
	plot(unwrap(angle(Fs1_rc(:,48))),'k','linewidth',1)
	xlim([0 length(Fs1_1(:,48))])
	xlabel('$t/\Delta t$', 'Interpreter', 'latex')
	ylabel('phase', 'Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','1')
	%{
	purinto(Fs1_1)
	xlabel('$f_\tau / \Delta f_\tau$', 'Interpreter', 'latex')
	ylabel('$t/\Delta t$', 'Interpreter', 'latex')
	caxis([0 160])
	%export_fig 1.jpg

	
	
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
	
	%  Filter h_c - RCMC for range curveture 
    rcc_f = 0;
	if rcc_f  
		f_tau = linspace(0, 1/fsamp, 2^tau_nt2 ); 
			
		H_c =  exp(  1j* (2* pi /c )* ( (v_y-v_p)^2 / R_0) *eta'.^2 * f_tau); %range curveture filter
		%Multiply H_c to Fs_r
		Fs1_1 = Fs1_rc .* H_c ;
		Fs2_1 = Fs2_rc .* H_c ;
		s1_1 = ifft( Fs1_1,[],2);
		s2_1 = ifft( Fs2_1,[],2);
		
		%purinto(s1_1) 
		%xlabel('$\tau /\Delta \tau$', 'Interpreter', 'latex')
		%ylabel('$t/\Delta t$', 'Interpreter', 'latex')
	end 
	K_a = 2 * f_0/c * (v_y-v_p)^2 / R_0;
	(v_x * x_0 + y_0*(v_p - v_y))/R_0 ;
	plot((v_x * x_0 + y_0*(v_p - v_y))/R_0 * eta + 2 * f_0/c * (v_y-v_p)^2 / R_0 * eta.^2)

	clear Fs1_rc Fs2_rc

	%% Radon transform to estimate the slop
	% real v_r  and the Doppler frequency
	v_r_real_ = v_x * (h/R_0/sqrt(1-y_0^2/R_0^2) );
	%fprintf('Doppler frequency %f\n',f_dc)
	
	do_radon = 1;
	if do_radon 
		iptsetpref('ImshowAxesVisible','on')
		theta = [177:0.005:181];
		parral = 0;
		if parral 
			[R,xp] = radon(gpuArray(abs(s1_1)),theta); % Use GPU to compute
			R = gather(R); xp = gather(xp);
		else
			[R,xp] = radon(abs(s1_1),theta);
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
	else
		[~, down] = max(abs(s1_1(1,:)) );
		[~, up] = max(abs(s1_1(length(s1_1(:,1)),:)));
		x_sl = (up - down) / 4 / B / dur  ;    
		v_rt = c * x_sl ;
	end
    %% RCMC for range walk 
	%t = repmat(linspace(-t(end,end),t(end,end),2*floor(PRF*t(end,end))),length(f_tau),1).';
    if do_radon 
		%Fh_rw = exp(j * 2 * pi  * ell * t .* repmat(f_tau, length(t(:,1)), 1) );    
		Fh_rw = exp(j * 2 * pi  * ell * linspace(t(1,end), t(end,end), length(s1_1(:,1))).' * f_tau );    
	else 
		Fh_rw = exp(j * 4 * pi / c * t / 2 * v_rt .* repmat(f_tau, length(eta), 1) );    
	end
    Fs1_2 = Fs1_1 .* Fh_rw ; 
	Fs2_2 = Fs2_1 .* Fh_rw ; 
    s1_2 = ifft(Fs1_2.').';
	s2_2 = ifft(Fs2_2.').';
	
	purinto(s1_2)
	title(['$(v_x,v_y)= ($' int2str(v_x) ',' int2str(v_y) '$)$'],'Interpreter', 'latex')
	xlabel('$\tau/ \Delta \tau$','Interpreter', 'latex')
	ylabel('$ t / \Delta t $','Interpreter', 'latex')
	%export_fig(['(' int2str(v_x) '_' int2str(v_y) ')_rw.jpg'])
	%{%}
	clear Fh_rw
	%% Search the parameters by phase 
	%f = fittype('a*x+b');	
	%[fit1,~,~] = fit(eta.', unwrap(angle(s1_2(:,148).*conj(s2_2(:,148)))), f,'StartPoint',[1 1]);
	%v_yt = v_p + fit1.a * lambda / 2 / pi *R_0 / d_a;
	%fprintf('Method 0, Actual: %f, Estimate: %f\n', v_y, v_yt)
	%clear f fit1
		

		L =201;
		ind = 295; %588; 148; 
		filt_ = [hamming(L).'/sum(hamming(L)); rectwin(L).'/L];
		for gg = 1 : 1
			s_filter = ifft(fft(unwrap(angle(s1_2(:,ind).*conj(s2_2(:,ind))))) .*  fft(filt_(gg,:),length(s1_2(:,ind))).' );
			
			figure
			plot(s_filter,'k','Linewidth',2)
			%hold on
			%plot(unwrap(angle(s1_2(:,ind).*conj(s2_2(:,ind)))),'Linewidth',2)
			%xlabel('$\Delta t$', 'Interpreter', 'latex')
			%ylabel('Phase', 'Interpreter', 'latex')
			%plot_para('Maximize',1,'Filename',int2str(gg))
			
			t2fit = linspace(t(1,end), t(end,end), length(s1_1(:,1))) ;
			[fit1,~,~] = fit(t2fit(L:end).',s_filter(L:end) ,'poly1');
			v_yt = v_p + fit1.p1 * lambda / 2 / pi *R_0 / d_a;
			%fprintf('Method %d, Actual: %f, Estimate: %f\n', gg, v_y, v_yt)
			
			%figure
			%plot(fit1,t(L:end,148),s_filter(L:end))
			%xlabel('$t$', 'Interpreter', 'latex')
			%ylabel('Phase', 'Interpreter', 'latex')
			%plot_para('Maximize',1,'Filename','fit')
			clear f fit1
		end
end
%{	
	%%
	ind = 148;
	plot(ifft(fft(unwrap(angle(s1_2(:,ind).*conj(s2_2(:,ind))))) .*  fft(filt_(gg,:),length(eta)).'  ),'Linewidth',1)
	hold on 
	plot(ifft(fft(unwrap(angle(s1_2(:,ind+1).*conj(s2_2(:,ind+1))))) .*  fft(filt_(gg,:),length(eta)).'  ),'Linewidth',1)
	hold on
	plot(ifft(fft(unwrap(angle(s1_2(:,ind+2).*conj(s2_2(:,ind+2))))) .*  fft(filt_(gg,:),length(eta)).'  ),'Linewidth',1)
	plot_para(1,1,'fit')

%%
close all
rrr= exp(j*2*pi/lambda* (d_a * y_0 / R_0 - d_a*(v_p - v_y)/R_0 * t(:,148)));  
I2 = exp(-j * 4* pi /lambda *(R_0 + d_a * y_0 / 2 / R_0 - d_a*(v_p - v_y)/ 2 / R_0 *t(:,148) ...
        + (v_x * x_0 - y_0 *(v_p - v_y))/ R_0 * t(:,148) + ((v_y - v_p)^2 + a_x * x_0)/2/R_0 * t(:,148).^2 ));
I1 = exp(-j * 4* pi /lambda *(R_0 + ...
        + (v_x * x_0 - y_0 *(v_p - v_y))/ R_0 * t(:,148) + ((v_y - v_p)^2 + a_x * x_0)/2/R_0 * t(:,148).^2 ));

	
figure 
plot(angle(I1))
hold on 
plot(angle(I2))
plot_para(1,1,'1')

figure 
plot(unwrap(angle(I1)))
hold on 
plot(unwrap(angle(I2)))
plot_para(1,1,'2')

figure
yyaxis left
plot(angle(I1) - angle(I2),'Linewidth',2)
yyaxis right
plot(unwrap(angle(I1)) - unwrap(angle(I2)),'Linewidth',2)
plot_para(1,1,'3')


close all
temp = 1 : length(s1_2(:,148));

figure
plot(angle(s1_2(:,148)),'Linewidth',1)
hold on 
plot(angle(s2_2(:,148)),'r','Linewidth',1)
plot_para(1,1,'4')

figure
plot(unwrap(angle(s1_2(:,148))),'Linewidth',1)
hold on 
plot(unwrap(angle(s2_2(:,148))),'r','Linewidth',1)
plot_para(1,1,'5')

figure
yyaxis left
plot(angle(s1_2(:,148)) - angle(s2_2(:,148)),'Linewidth',1)
yyaxis right
plot(unwrap(angle(s1_2(:,148))) - unwrap(angle(s2_2(:,148))),'Linewidth',1)
plot_para(1,1,'6')

figure
plot(angle(s1_2(:,148).*conj(s2_2(:,148))),'Linewidth',1)
hold on 
plot(unwrap(angle(s1_2(:,148))) - unwrap(angle(s2_2(:,148))),'r','Linewidth',1)
plot_para(1,1,'7')

	
figure
plot(t(:,148),angle(s1_2(:,148).*conj(s2_2(:,148))),'k','Linewidth',2)
xlabel('t')
ylabel('Phase')
plot_para(1,1,'7')
%}

%{
	%% Spectrum of the filter
	L = 101;
	lon = 2^10;
	filter = [rectwin(L).'; triang(L).'; hamming(L).'];
	f_spect = linspace(-0.5,0.5,lon);
	figure 
	plot(f_spect,10*log10(fftshift(abs(fft(filter(1,:),lon)  .* conj(fft(filter(1,:),lon))))/max(abs(fft(filter(1,:),lon) .* conj(fft(filter(1,:),lon))))) ,'Linewidth',2)
	hold on 
	plot(f_spect,10*log10(fftshift(abs(fft(filter(2,:),lon)  .* conj(fft(filter(2,:),lon))))/max(abs(fft(filter(2,:),lon) .* conj(fft(filter(2,:),lon))))) ,'Linewidth',2)
	hold on 
	plot(f_spect,10*log10(fftshift(abs(fft(filter(3,:),lon)  .* conj(fft(filter(3,:),lon))))/max(abs(fft(filter(3,:),lon) .* conj(fft(filter(3,:),lon))))) ,'k','Linewidth',2)
	ylim([-50 0])
	xlabel('Normalized frequency','Interpreter', 'latex')
	ylabel('Spectrum (dB)','Interpreter', 'latex')
	plot_para(1,1,'f_respon')
	
	f_spect = linspace(-0.5,0.5,length(eta));
	plot(f_spect,10*log10(fftshift(abs(fft(unwrap(angle(s1_2(:,148).*conj(s2_2(:,148)))))  .* conj(fft(unwrap(angle(s1_2(:,148).*conj(s2_2(:,148))))))))/max(abs(fft(unwrap(angle(s1_2(:,148).*conj(s2_2(:,148))))) .* conj(fft(unwrap(angle(s1_2(:,148).*conj(s2_2(:,148))))))))),'k' ,'Linewidth',2)
	xlabel('Normalized frequency','Interpreter', 'latex')
	ylabel('Spectrum (dB)','Interpreter', 'latex')
	plot_para(1,1,'ATI_Pha_f')
	
%}

