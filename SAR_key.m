function [] = SAR_key(vx, vy) 

% Signal generation is coded in other function.
% This function implement SAR from range compression, 
% Keystone transform, which is for range curvature correction. Then estimate the slope of the line to 
% apply range walk correction. The last step is azimuth compression.
% 
% Usage: SAR_key(vx, vy), where "vx" is range velocity and "vy" is azimuth
% velocity.
	
	clear 
    close all
    delete parameter.mat
    if ispc 
		cd D:\Code\Simu
    else
		cd /Matlab 
	end
    %% Signal 
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
    S_rc = fftshift(S_rc, 2) ; % Not sure if "fftshift" is needed ???????? ?çÊ??ª„ÅØÂøÖË?Ôº?
    % Interpolation 
    for i = 1 : length(f_tau)
    	t(:, i).' = eta * sqrt(f_0 / f_0 + f_tau(i) ) ;    
    end
    for m = 1 : length(f_tau)
    	for n = 1 : length(eta)
    		if t(1) < eta(n) && t(length(eta)) > eta(n)
    			S_1(m,n) = S_rc(n,) ;
    		elseif t(1) > eta(n) 
    			continue
    		else
    			break
    		end
    	end
    end

    %% RCMC for range walk 
    v_rt = 0 ; % ‰øÆÊ≠£Ë¶Å„?
    H_rw = exp(j * 4 * pi / c * f_tau / 2 *??_rt * t );
    s_2 = s_1 .* H_rw ; 
    %% Ambiguity function
    H_AF(s_2);

end