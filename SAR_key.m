function [] = SAR_key(vx, vy) 
    %% Signal 
    if ispc 
		cd D:\Code\Simu
	else if islinux
		cd /Matlab
	end
    s = Gen_signal();
    %Gen_signal(vx, vy)
    load('parameter.mat');
	
    %% Range Compression
	td = 2* R_0 / c ;
    ref_time = -Tp : 1/4/B : 0 ;
    s_0 = exp(i* 4* pi*R_0/ lambda - j* pi* Kr* (ref_time).^2)  ;
	
	R_m = sqrt(x_n^2 + (y_n - v_p * eta).^2 + h^2); % Static target range equation
	h_m = exp(j * 4 * pi * f_0 * R_m / c) * exp(- j * pi * K_r * tau.^2) ; % Matched filter.
	
    M = nextpow2(length(tau)) ;  % Make total number to power of 2
    Fh_m = fft( h_m , 2^M );
    Fs = fft(s.',2^M).'; 
    Fs_rc = Fs .* Fh_m; % Range compression
	
    s_rc = ifft( Fs_rc.' ).';
    %purinto(s_r); % Print the s_r
    
   
    %% Key Stone transform - RCMC for range curveture 
    eta^2 = f_0 / f_0 + f_tau t^2
	


end