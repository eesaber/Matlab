function [] = SAR_key(vx, vy) 
    %% Signal 
    clear 
    close all
    delete parameter.mat
    if ispc 
		cd D:\Code\Simu
    else
		cd /Matlab
	end
    s = Gen_signal();
    %s = Gen_signal(vx, vy);
    load('parameter.mat');
    %% Range Compression
    ref_time = -T_p : 1/4/B : 0 ;
	R_m = sqrt(x_n^2 + (y_n - v_p * eta).^2 + h^2); % Static target range equation
	% Matched filter.
    for i = 1 : length(eta)
        h_m(i,:) = exp(j * 4 * pi * f_0 * R_m(i) / c).' * ( exp(- j * pi * K_r * tau.^2) .* ( tau >= 2* R_m(i) / c & tau <= 2* R_m(i) / c + T_p) ); 
    end
	h_m = flip(h_m, 2);
	purinto(abs(h_m))
    purinto(abs(s))
    %%
    M = nextpow2(length(tau)) ;  % Make total number to power of 2
    Fh_m = fft( h_m.', 2^M ).';
    Fs = fft(s.', 2^M).'; 
    Fs_rc = Fs .* Fh_m; % Range compression
    s_rc = ifft( Fs_rc.' ).';
    purinto(abs(s_rc)); % Print the s_r
    %% Key Stone transform - RCMC for range curveture 
    %eta * sqrt(f_0 / f_0 + f_tau) t^2 ;
    
    
    %% RCMC for range walk 


end