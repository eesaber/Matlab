function [r_PLSR, a_PLSR, rsandb, asandb, t_v_x, t_v_y] = SAR_rda (vx,vy) 
    %% Signal 
    %Gen_signal(vx, vy);
    s = Gen_signal(10, 0);
    load('parameter.mat');
    %% Range compression 
    % Reference signal 
    td = 2* R_0 / c ;
    ref_time = -T_p : 1/4/B : 0 ;
    s_0 = exp(i* 4* pi*R_0/ lambda - j* pi* K_r* (ref_time).^2)  ;
    % FFT to s_0 w.r.t tau
    M = nextpow2(length(tau)) ; 
    Fs_0 = fft( s_0, 2^M );
    % FFT to s w.r.t tau and "Do range compression"
    Fs = fft(s.',2^M).';
    for k = 1 : length(eta)
        Fs_r(k,:) = Fs(k,:) .* Fs_0 ;
    end
    s_r = ifft( Fs_r.' ).';
    purinto(s_r); % Print the s_r
    %% RCMC for range curveture 
    
    f_tau = 0  : 1 : 2^M - 1;
    f_tau = f_tau/(2^M/4/B) ;
    H_c =  exp(  j* (2* pi /c )* ( (0-v_p)^2 / R_0) *eta'.^2 * f_tau); %range curveture filter
    %Multiply H_c to Fs_r
    Fs_1 = Fs_r .* H_c ;
    s_1 = ifft( Fs_1.' ).';
    purinto(Fs_1); % Print the s_1
    purinto(s_1); % Print the s_1
    %% Compute v_r 
    control_rt = true;  % Control if the procedure go through hough transform 
    if control_rt == true
        [temp low_index] = max(abs(s_1(1,:)), [], 2);
    	[temp up_index] = max(abs(s_1( size(s_1(:,1)),: )), [], 2); 
        weiyi = up_index - low_index;
        kai = 1 / 4 / B / dur * weiyi(1) ;
        t_v_r = c / 2 * kai ;
        t_v_x = t_v_r * R_0 / x_n ;
        2 * t_v_r / lambda ;
         % Apply Hough tranform to fint the angle.
    end

    %% RCMC for range walk
    % ï¿½ï¿½ï¿½ï¿½ï¿?
    god = false;
    if god == true
        H_w = exp(1i* 4* pi /c * v_x / R_0 * x_n * eta' * f_tau);
    else 
        H_w = exp(1i* 4* pi /c * t_v_r * eta' * f_tau);
    end
    Fs_2 = Fs_1 .* H_w;  
    s_2 = ifft( Fs_2.' ).';
    %purinto(s_2); % Print the s_2

    %% WVD (spectrum)
    [dum Maxran ] = max(abs(s_2(3000,:)));
    figure
    temp= spectrogram(s_2(:,Maxran),128,120,8196,1000,'centered','yaxis');
    x = [-3 3];
    y = [-500 500];
    [dum , upf] = max(temp(:,1));
    [dum , dwf] = max(temp(:,length(temp(1,:))));
    t_K_a = (dwf - upf) * 1000 / 8196 / dur ;
    t_v_y = 100 - sqrt( - t_K_a * c * R_0 / f0 / 2);
    if imag(t_v_y) ~= 0
        [dum , upf] = max(temp(:,length(temp(1,:)) - 10 ));
        [dum , dwf] = max(temp(:,length(temp(1,:)) ));
        t_K_a = (dwf - upf) * 1000 / 8196 / dur;
        t_v_y = 100 - sqrt( - t_K_a * c * R_0 / f0 / 2);
    end
    imagesc(x, y, abs(temp))
    set(gca,'Ydir','normal')
    set(gca,'FontSize',32,'Fontname','CMU Serif')
    xlabel('$\eta$','Interpreter','latex')
    ylabel('$f_\eta$','Interpreter','latex')
    colormap('Jet')
    colorbar
    pbaspect([4 3 1])

    %% Azimuth compression 
    %[t_v_x t_v_y]
    control_ac = true ;
    if control_ac == true 
    	% Azimuth reference signal
    	eta0 = linspace(-dur/2, dur/2,length(s_2(:,1)) );
    	K_a = 2/lambda *(v_y - v_p)^2 / R_0 ;
    	s_a0 = exp( -j * 4 * pi * t_v_r / lambda * eta0 + j * pi* K_a* eta0.^2 ) ;
    	M = nextpow2(length(eta0));
    	Fs_a0 = fft(s_a0,2^M)'; 
    	Fs_a0 = conj(Fs_a0); 
        Fs_2a = fft(s_2 ,2^M); 
    	for k = 1 : length(s_2(1,:))
            Fs_a(:,k) = Fs_2a(:,k).* Fs_a0 ;
        end
        s_a = ifft(Fs_a);
        purinto(s_a)
        %clearvars eta0 M eta tau
    end
close all
    %% Pulse sidelobe ratio (PLSR) - Right sidelobe 
    %Azimuth 
    clear temp 
    [dum , max_aind] = max(abs(s_a(:,Maxran)));
    temp = diff(abs(s_a(:,Maxran)));
    coun = 0;
    for k = max_aind : length(s_a(:,1))
        if temp(k) * temp(k+1) < 0
            coun = coun + 1;
        end
        if coun == 2
            a_PLSR = 10* log10(abs(s_a(k + 1, Maxran)) / abs(s_a(max_aind, Maxran)));
            break
        end
    end
    % Range
    clear temp
    coun = 0;
    [dum max_rind] = max(abs(s_a(max_aind,:)));
    temp = diff(abs(s_a(max_aind,:)));
    for k = max_rind : length(s_a(1,:))
        if temp(k) * temp(k + 1) < 0
            coun = coun + 1;
        end
        if coun == 2
            r_PLSR = 10 * log10(abs(s_a(max_aind, k + 1) ) / abs(s_a(max_aind, Maxran)) );
            break
        end
    end

    %% Range 3dB  azimuth 3dB
    % azimuth
    for k = max_aind : length(s_a(:,1))
        if abs(s_a(max_aind, Maxran)) / 2  >= abs(s_a(k, Maxran))
            asandb = (k - max_aind) * v_p / PRF;
            break
        end
    end
    % range 
    for k = max_rind : length(s_a(1,:))
        if abs(s_a(max_aind, max_rind)) / 2 >= abs(s_a(max_aind, k))
            rsandb = (k - max_rind) * (upran - downran) / length(tau);
            break
        end
    end

    

end