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
    s = Gen_signal(10,0);
    %s = Gen_signal(vx, vy); 
    load('parameter.mat');
    
    %% Range Compression
    ref_time = -T_p : 1/4/B : 0 ;
	R_m = sqrt(x_n^2 + (y_n - v_p * eta).^2 + h^2); % Static target range equation
	% Matched filter.
    for i = 1 : length(eta)
        h_m(i,:) = exp(j * 4 * pi * f_0 * R_m(length(eta)/2) / c) *  exp(- j * pi * K_r * ref_time.^2) ; 
    end
	%purinto(h_m)
    %purinto(s)
    tau_nt2 = nextpow2(length(tau)) ;  % Make total number to power of 2
    Fh_m = fft( h_m.', 2^tau_nt2 ).';
    Fs = fft(s.', 2^tau_nt2).'; 
    Fs_rc = Fs .* Fh_m; % Range compression
    s_rc = ifft( Fs_rc.' ).';
    %purinto(s_rc); % Print the s_r
    
    %% Key Stone transform - RCMC for range curveture 
    %Fs_rc = fftshift(Fs_rc, 2) ; % Not sure if "fftshift" is needed ????!!!!!!
    % Interpolation 
    f_tau = 0 : 2 * B / 2^tau_nt2 : 2*B -1 ;
    for i = 1 : 2^tau_nt2
    	t(:, i) = eta * sqrt(f_0 / (f_0 + f_tau(i)) ).' ;    
    end
    
    for m = 1 : 2^tau_nt2
    	Fs_1(:, m) = interp1(eta, Fs_rc(:,m).', t(:,m).', 'spline', 0).';
    	%{
    	for n = 1 : length(eta)
    		if t(1) < eta(n) && t(length(eta)) > eta(n)
    			[~, idx] = min(abs(eta(n) - t(:,m));
    			idx_u = idx + (eta(n) > t(:,idx)) ; 
    			idx_l = idx_u - 1;
    			S_1(n,m) = ((eta(n) - t(:,idx_l)) * S_rc(idx_l, m) + (t(:,idx_u) - eta(n)) * S_rc(idx_l, m))/ (t(:,idx_u) - t(:,idx_l)) ;
    		elseif t(1) > eta(n) 
    			continue
    		else
    			break
    		end
    	end
    	%}
    end
    s_1 = ifft(ifftshift(Fs_1, 2).').'; % If the central of f_tau is 0, then ifftshift is needed. If not so, do not do ifftshift 
    
 
    %% RCMC for range walk 
    v_rt = 10 ; % 修正?
    Fh_rw = exp(j * 4 * pi / c * t / 2 * v_rt .* repmat(f_tau, length(eta), 1) );
    Fs_2 = Fs_1 .* Fh_rw ; 
    s_2 = ifft(Fs_2.').';
    
    
   
    purinto(s_1);
    purinto(s_2);
    %% Ambiguity function
    H_AF(s_2);

end