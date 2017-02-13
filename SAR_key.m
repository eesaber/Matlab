function [v_r, v_rt, v_yt] = SAR_key(vx, vy, ax, ay) 
% Signal generation is coded in other function.
% This function implement SAR from range compression, 
% Keystone transform, which is for range curvature correction. Then estimate the slope of the line to 
% apply range walk correction. The last step is azimuth compression.
% 
% Usage: SAR_key(vx, vy), where "vx" is range velocity and "vy" is azimuth
% velocity.
	
	%clear 
    %clc
    %close all
    %delete parameter.mat
    if ispc 
		cd D:\Code\Simu
    else
		cd ~/Matlab 
	end
    %% Signal 
    s = Gen_signal(-10,0,0,0);
    %s = Gen_signal(vx, vy, ax, ay); 
    load('parameter.mat');
    ref_time = -T_p : 1/4/B : 0 ;
	R_m = sqrt(x_n^2 + (y_n - v_p * eta).^2 + h^2); % Static target range equation
    plot(R_m,'k')
    hold on 
    plot(R)
    %% Range Compressio
	% Matched filter.
    for i = 1 : length(eta)
        h_m(i,:) = exp(j * 4 * pi * f_0 * R_m(length(eta)/2) / c) *  exp(- j * pi * K_r * ref_time.^2) ; 
    end
    sprintf('這是用固定的R_0 不是用R(eta)')
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
    Fh_rw = exp(j * 4 * pi / c * t / 2 * v_rt .* repmat(f_tau, length(eta), 1) );    
    Fs_2 = Fs_1 .* Fh_rw ; 
    if shift == 1
        s_2 = ifft(ifftshift(Fs_2, 2).').';
    else
        s_2 = ifft(Fs_2.').';
    end
    purinto(s_2)
    export_fig s_2.jpg
    fprintf(' v_r: %f , v_rt: %f \n Estimation error: %f \n ', v_r , v_rt, v_r - v_rt )
    
    %% Ambiguity function
    xi_4 = H_AF(s_2(:, 663)',4);
    purinto(xi_4, '$t$(s)', '$f_\xi$(Hz)', linspace(-1, 1, 8192), linspace(1000, -1000 ,4096) )
    %xi_44 = H_AF(s_2(:, 664)',4);
    %purinto(xi_44, '$t$(s)', '$f_\xi$(Hz)', linspace(-1, 1, 8192), linspace(1000, -1000 ,4096) )
    
    % Analyze the coefficient
    co_0 = 4 * pi / lambda * R_0;
    co_1 = 4 * pi / lambda * v_rt; 
    co_2 = 4 * pi / lambda * ((v_y -v_p)^2 + a_x * x_n) / 2 / R_0;
    co_3 = 4 * pi / lambda * v_rt^2 * (v_y - v_p)^2 / 2 / R_0^2 ;
    fprintf('The coefficient of the signal is: \n  co_0 = %d , co_1 = %d \n  co_2 = %d , cp_3 = %d \n', co_0, co_1, co_2, co_3);
    xi_0 = dur * 0.089 ;
    y_sl = -24 * xi_0 * f_0 / c * v_rt * (v_y - v_p)^2 / 2 / R_0^2;
    f_int = - 8 * xi_0 / lambda * co_2 / 4 / pi * lambda;
    fprintf('From the analytic solution\n  Slope: %f Intercept: %f, \n ', y_sl, f_int); 
    y_sl = (239.1 - 237.1) / (0.1193 + 0.05843)
    v_yt = v_p - sqrt(- y_sl * c * 2 * R_0^2 / (24 * xi_0 * f_0 * v_rt)) ;
    fprintf('Estimated v_y: %f m/s \n', v_yt)
    
    %{
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