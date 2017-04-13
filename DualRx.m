function DualRx()
% This function generate the SAR signal regarding to the input parameter.
% Usage: Gen_signal(vx, vy, ax, ay), 'vx' is range velocity, 'vy' is azimuth
% velocity, 'ax' is range accelaration and 'ay' is azimuth accelration. If no input parameters, the velocity in both direction are set
% to zero
% Noth that the parameters of SAR need to be modified in the function.

    delete parameter.mat
    if nargin == 0
        vx = 0; vy = 0; 
        ax = 0; ay = 0;
    end


    %% Parameters - C Band airborne SAR parameters
    % Target
    x_0 = 2160;
    y_0 = 400;
	
    d = 100; % Length of the target area 
    v_x = vx; a_x = ax; % rangecl
    v_y = vy; a_y = ay; % azimuth 
    % Platform
    h = 2200;
	R_0 = sqrt(x_0^2 + y_0^2 + h^2);
	l_0 = sqrt( y_0^2 + h^2);
    d_a = 0.2;     
    v_p = 90; 
    f_0 = 9.6e9; c = 3e8 ; lambda = c/f_0 ;
    dur = 2 ; %
    PRF = 1000; %%2.35e3
    K_r = 5e14 ;
    T_p = 0.1e-6; % Pulse width
    B =  K_r*T_p;

    %% Signal 
    aa = 0;
    eta = linspace( aa - dur/2, aa + dur/2,PRF*dur) ;  % Slow time -6

    upran = sqrt( (x_0 + d)^2 + h^2) ; downran =  sqrt( (x_0 - d)^2 + h^2); %Set upper limit and down limit 
    tau =  2*(downran)/c: 1/4/B : 2*(upran)/c + T_p  ;  % fast time space

    R1 = sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
    R2 = sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (d_a + y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
   
    s1 = zeros(length(eta), length(tau));
	s2 = zeros(length(eta), length(tau));
    for k = 1 : length(eta)
        td1 = 2* R1(k)/ c ;
		td2 = (R1(k) + R2(k)) /c ;
        s1(k,:) = exp(-i* 4* pi*R1(k)/ lambda + i* pi* K_r* (tau - 2* R1(k)/ c).^2) .*(tau>=td1 & tau-td1<=T_p)  ;
		s2(k,:) = exp(-i* 2* pi*(R1(k) + R2(k))/ lambda + i* pi* K_r* (tau - (R1(k) + R2(k) )/ c).^2) .*(tau>=td2 & tau-td2<=T_p)  ;
	end
	clear R1 R2
	%% range compression 
	ref_time = -T_p : 1/4/B : 0 ;
	h_m = repmat(exp(j * 4 * pi * f_0 * R_0 / c) *  exp(- j * pi * K_r * ref_time.^2),length(eta),1 ) ; 

	tau_nt2 = nextpow2(length(tau)) ;  % Make total number to power of 2
    Fh_m = fft( h_m.', 2^tau_nt2 ).'; 
    Fs1_rc = fft(s1.', 2^tau_nt2).' .* Fh_m; % Range compression
    s1_rc = ifft( Fs1_rc.' ).';
    Fs2_rc = fft(s2.', 2^tau_nt2).' .* Fh_m; % Range compression
	s2_rc = ifft( Fs2_rc.' ).';
	clear Fh_m s1 s2
	%%
	% Key Stone transform - RCMC for range curveture 
    % Interpolation 
    f_tau = linspace(0, 4*B -1, 2^tau_nt2 ); 
    shift = 0; % FFTshift false
    if shift == 1
    	f_tau = linspace(-2*B, 2*B -1, 2^tau_nt2 ); 
    	Fs_rc = fftshift(Fs_rc, 2);
    end
    for m = 1 : 2^tau_nt2
    	t(:, m) = eta * sqrt(f_0 / (f_0 + f_tau(m)) ).' ;    
	end
	
	
    for m = 1 : 2^tau_nt2
    	Fs1_1(:, m) = interp1(eta, Fs1_rc(:,m).', t(:,m).', 'spline', 0).';
		Fs2_1(:, m) = interp1(eta, Fs2_rc(:,m).', t(:,m).', 'spline', 0).';
    end
    if shift == 1
    	s1_1 = ifft(ifftshift(Fs1_1, 2).').'; % If the central of f_tau is 0, then ifftshift is needed. If not so, do not do ifftshift 
		s2_1 = ifft(ifftshift(Fs2_1, 2).').';
    else
    	s1_1 = ifft(Fs1_1.').';
		s2_1 = ifft(Fs2_1.').';
	end
	clear Fs1_rc Fs2_rc
	%% Estimate the slope 
	
    [~, down] = max(abs(s1_1(1,:)) );
    [~, up] = max(abs(s1_1(length(s1_1(:,1)),:)));
    x_sl = (up - down) / 4 / B / dur  ;    
    v_rt = c * x_sl ;
  		
    %% RCMC for range walk 
    v_r = v_x * x_0 / R_0;
    Fh_rw = exp(j * 4 * pi / c * t / 2 * v_rt .* repmat(f_tau, length(eta), 1) );    
    Fs1_2 = Fs1_1 .* Fh_rw ; 
	Fs2_2 = Fs2_1 .* Fh_rw ; 
    if shift == 1
        s1_2 = ifft(ifftshift(Fs1_2, 2).').';
		s2_2 = ifft(ifftshift(Fs2_2, 2).').';
    else
        s1_2 = ifft(Fs1_2.').';
		s2_2 = ifft(Fs2_2.').';
	end
	%% Search the parameters by phase 
	f = fittype('a*x+b');
	[fit1,~,~] = fit(t(:,148),unwrap(angle(s1_2(:,148).*conj(s2_2(:,148)))),f,'StartPoint',[1 1]);
	v_yt = v_p + fit1.a * lambda / 2 / pi *R_0 / d_a;
	%fprintf('Actual: %f, Estimate: %f\n', v_y, v_yt)
	fprintf('%f\n, ',v_yt)
	clear f fit1
	
    %% Search the parameters by Fourier transform
	N = 2^18;
	[~, ind] = max(abs(fftshift(fft(s1_2(:,1468).*conj(s2_2(:,1468)), N ))));
	temp =linspace(-PRF/2,PRF/2, N) ;
	f_tp = temp(ind);
    clear temp
    v_yt = v_p + f_tp * R_0 * c /(2 * f_0 * d_a);
	fprintf('The ideal f_tp: %f, The estimation f_tp: %f', -2 * d_a / lambda * d_a *  (v_p - v_y) / R_0, f_tp)
    
	plot(linspace(-20,20,41), v_yt -linspace(-20,20,41),'k','Linewidth',3)
    xlabel('$v_y$', 'Interpreter', 'latex')
    ylabel('$\tilde{v}_y - v_y$', 'Interpreter', 'latex')
    plot_para(1,1)
    

	%% Analyze the phase characteristic
	close all
	% Ideal signal 
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

end
