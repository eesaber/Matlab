function [] = SAR_key(vx, vy) 
    %% Signal 
    Gen_signal(vx, vy)
    load('parameter.mat')
    %% Range Compression
    td = 2* R_0 / c ;
    ref_time = -Tp : 1/4/B : 0 ;
    s_0 = exp(i* 4* pi*R_0/ lambda - j* pi* Kr* (ref_time).^2)  ;
    % FFT to s_0 w.r.t tau
    M = nextpow2(length(tau)) ; 
    Fs_0 = fft( s_0, 2^M );
    % FFT to s w.r.t tau and "Do range compression"
    Fs = fft(s.',2^M).';
    for k = 1 : length(eta)
        Fs_r(k,:) = Fs(k,:) .* Fs_0 ;
    end
    s_r = ifft( Fs_r.' ).';
    %purinto(s_r); % Print the s_r
    %% RCMC for range curveture 
    % f_tau 
    f_tau = 0  : 1 : 2^M - 1;
    f_tau = f_tau/(2^M/4/B) ;
    H_c =  exp(  j* (2* pi /c )* ( (0-v_p)^2 / R_0) *eta'.^2 * f_tau); %range curveture filter
    %Multiply H_c to Fs_r
    Fs_1 = Fs_r .* H_c ;
    s_1 = ifft( Fs_1.' ).'; 
    %% Key Stone transform 
    

    %%

end