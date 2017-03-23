function [co_3, c_3, v_rt, v_yt, a_rt] = SAR_key(vx, vy, ax, ay) 
% Signal generation is coded in other function.
% This function implement SAR from range compression, 
% Keystone transform, which is for range curvature correction. Then estimate the slope of the line to 
% apply range walk correction. The last step is azimuth compression.
% 
% Usage: SAR_key(vx, vy), where "vx" is range velocity and "vy" is azimuth
% velocity.
	
	%clear 
    %clc
    close all
    %delete parameter.mat
    if ispc 
		cd D:\Code\Simu
    else
		cd ~/Code/Matlab 
	end
    %% Signal 
    %s = Gen_signal(20,-10,5,0);
	s = Gen_signal(vx, vy, ax, ay); 
    load('parameter.mat');
    ref_time = -T_p : 1/4/B : 0 ;
	R_m = sqrt(x_n^2 + (y_n - v_p * eta).^2 + h^2); % Static target range equation
    %plot(R_m,'k')
    %hold on 
    %plot(R)
    %% Range Compression
	% Matched filter.
    for i = 1 : length(eta)
        h_m(i,:) = exp(j * 4 * pi * f_0 * R_m(length(eta)/2) / c) *  exp(- j * pi * K_r * ref_time.^2) ; 
    end
    %sprintf('?�是?�固定�?R_0 不是?�R(eta)')
	%purinto(h_m)
    %purinto(s)
    tau_nt2 = nextpow2(length(tau)) ;  % Make total number to power of 2
    Fh_m = fft( h_m.', 2^tau_nt2 ).';
    Fs = fft(s.', 2^tau_nt2).'; 
    Fs_rc = Fs .* Fh_m; % Range compression
    s_rc = ifft( Fs_rc.' ).';
    %purinto(s_rc); % Print the s_r
    
    %% Key Stone transform - RCMC for range curveture 
    % Interpolation 
    f_tau = linspace(0, 4*B -1, 2^tau_nt2 ); 
    shift = 0; % FFTshift false
    if shift == 1
    	f_tau = linspace(-2*B, 2*B -1, 2^tau_nt2 ); 
    	Fs_rc = fftshift(Fs_rc, 2);
    end
    for i = 1 : 2^tau_nt2
    	t(:, i) = eta * sqrt(f_0 / (f_0 + f_tau(i)) ).' ;    
    end
    for m = 1 : 2^tau_nt2
    	Fs_1(:, m) = interp1(eta, Fs_rc(:,m).', t(:,m).', 'spline', 0).';
    end
    if shift == 1
    	s_1 = ifft(ifftshift(Fs_1, 2).').'; % If the central of f_tau is 0, then ifftshift is needed. If not so, do not do ifftshift 
    else
    	s_1 = ifft(Fs_1.').';
    end
    %purinto(s_1);
    [~, down] = max(abs(s_1(1,:)) );
    [~, up] = max(abs(s_1(length(s_1),:)));
    x_sl = (up - down) / 4 / B / dur  ;    
    v_rt = c * x_sl ;
    
    %% RCMC for range walk 
    v_r = v_x * x_n / R_0; 
    Fh_rw = exp(j * 4 * pi / c * t / 2 * v_r .* repmat(f_tau, length(eta), 1) );    
    Fs_2 = Fs_1 .* Fh_rw ; 
    if shift == 1
        s_2 = ifft(ifftshift(Fs_2, 2).').';
    else
        s_2 = ifft(Fs_2.').';
	end
	
    %purinto(s_2)
    %export_fig s_2.jpg
    fprintf(' v_r: %f , v_rt: %f Estimation error: %f \n ', v_r , v_rt, v_r - v_rt )
	
%% General Ambiguity function (m = 3)
	%Chose the maximum value when R = R_0
	[~, index] = max(abs(s_2(:)));
	[~, I_col] = ind2sub(size(s_2),index);
	s_3 = s_2(:,I_col).';
	close all
	%plot(real(s_2(:,663)))
	
	%{
	AF = GAF(s_3,4,4,3);
	[~, f_t] = max(AF);
	xi = dur / 8 ;
	f = linspace(-PRF/2, PRF/2, length(AF)); 
	f(f_t) / 2^(3) / xi^(3)
	%}
	
	%s_3 = exp( j * 2 * pi * ( co_0 + co_1 * eta + co_2 * eta.^2 / 2 + co_3 * eta.^3 / 6 ) );
	order = 3;
	Ec = zeros(1, order + 1);
	for n =  order : -1 : 1
		AF = GAF(s_3,3,n,3);
		[~, f_t] = max(AF);
		xi = dur / 6 ;
		f = linspace(-PRF/2, PRF/2, length(AF)); 
		Ec(n+1) = f(f_t) / 2^(n-1) / xi^(n-1);
		s_3 = s_3 .* exp(- j*2*pi * Ec(n+1) / factorial(n) * eta.^(n));
		%figure
		%plot(f,abs(AF))
	end
	v_yt = v_p - sqrt( Ec(4) * c / f_0 * R_0^2 / 6 / v_rt); % There are bug.
	c_3 = Ec(4);
	fprintf('v_y: %f, t_vy: %f , error: %f \n',v_y, v_yt , v_y - v_yt );
	a_rt = - Ec(3) -2 / lambda * (v_yt -v_p)^2 / R_0;
	a_r = a_x * x_n / R_0;
	fprintf('a_r: %f, t_ar: %f , error: %f \n', a_r, a_rt , a_r - a_rt );
	fprintf('------------------------------------------\n');
	
	% Analyze the coefficient
    co_0 = -2 / lambda * R_0;
    co_1 = -2 / lambda * v_r; 
    co_2 = -2 / lambda * ((v_y -v_p)^2 + a_x * x_n) / R_0;
    co_3 = 1 / lambda * 6* v_r * (v_y - v_p)^2 / R_0^2 ;
	%ahq = [co_0 co_1 co_2 co_3];
    %refSignal = exp(j * 2 * pi * ( co_0 + co_1 * eta + co_2 * eta.^2 / 2 + co_3 * eta.^3 / 6 ) );
    %{
	%%
    close all
    figure
	subplot(2,1,1) 
	plot(eta, real(refSignal))
	xlim([-.8 -0.5])
    subplot(2,1,2)
    plot(eta , -real(s_2(:,I_col))/max(abs(s_2(:,I_col))))
	xlim([-.8 -.5 ])
	%plot(eta,real(s_2(:,662)))
    %fprintf('The coefficient of the signal is: \n  co_0 = %d , co_1 = %d \n  co_2 = %d , co_3 = %d \n', co_0, co_1, co_2, co_3);
	
	%plot(real(k_3))
	% The simulation result
	figure
		plot(f,abs(AF_2),'k','Linewidth',3.5)
		set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
		set(gcf,'color','w');
		xlabel('Hz', 'Interpreter', 'latex')
		%xlim([-1000 0])
		%pause(0.00001);
		%frame_h = get(handle(gcf),'JavaFrame');
		%set(frame_h,'Maximized',1); 
		%export_fig AF_2.jpg
	figure
		plot(f,abs(AF_3),'k','Linewidth',3.5)
		set(gcf,'color','w');
		set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
		xlabel('Hz', 'Interpreter', 'latex')
		%xlim([-100 100])
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1); 
		%export_fig AF_3.jpg
    
   
    %% Hough trnaform
    I_rgb = imread('HAF.tif');
    
    I_cmyk = imread('HAF.tif');   
    inprof = iccread('USSheetfedCoated.icc');
    outprof = iccread('sRGB.icm');
    C = makecform('icc',inprof,outprof);
    I_rgb = applycform(I_cmyk,C);   
    
    I = rgb2gray(I_rgb);
    BW = edge(imrotate(I,0,'crop'),'canny');
    [H,T,R] = hough(BW);
    P  = houghpeaks(H,2);
    imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;
    plot(T(P(:,2)),R(P(:,1)),'s','color','white');
    
    lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);
    figure, imshow(I), hold on
    max_len = 0;
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
   
       end
    %}
end 