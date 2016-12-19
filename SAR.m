function [asandb, PLSR, rsandb, t_v_x, t_v_y] = SAR (vx, vy) 
%% Parameters - C Band airborne SAR parameters
% Target
x_n = 3000;
y_n = 0;
d = 200; % Length of the target area 
v_x = vx; a_x = 0; % rangecl
v_y = vy; a_y = 0; % azimuth 
% Platform
h = 8500;
La = 2;     
v_p = 100; 
f0 = 8.85e9; c = 3e8 ; lambda = c/f0 ;
dur = 3 ; %
PRF = 1000; %%2.35e3
Kr = 20e12 ;
Tp = 1.5e-6; % Pulse width
B =  Kr*Tp;
Ta = 0.886 * (lambda * x_n)/(La *v_p *1);

%% Signal 
aa = 0;
eta = linspace( aa - dur/2, aa + dur/2,PRF*dur) ;  % Slow time -6

upran = sqrt( (x_n + d)^2 + h^2) ; downran =  sqrt( (x_n - d)^2 + h^2); %Set upper limit and down limit 
tau =  2*(downran)/c: 1/4/B : 2*(upran)/c + Tp  ;  % fast time space

R = sqrt(h^2 + (x_n + v_x* eta + 0.5* a_x* eta.^2).^2 + (y_n + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
R_0 = R(length(R)/2);

for k = 1 : length(eta)
    td = 2* R(k)/ c ;
    s(k,:) = exp(-i* 4* pi*R(k)/ lambda + i* pi* Kr* (tau - 2* R(k)/ c).^2) .*(tau>=td & tau-td<=Tp)  ;
end
%purinto(s); % Print the s
%% Range compression 
% Reference signal 
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

%purinto(s_1); % Print the s_1
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
% ������� 
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
t_K_a = (dwf - upf) * 1000 / 8196 / dur;
t_v_y = 100 - sqrt( - t_K_a * c * R_0 / f0 / 2);
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

%% Pulse sidelobe ratio (PLSR) - Right sidelobe 
clear temp 
[dum , max_aind] = max(s_a(:,Maxran));
temp = diff(abs(s_a(:,Maxran)));
coun = 0;
for k = max_aind : length(s_a(:,1))
    if temp(k) * temp(k+1) < 0
        coun = coun + 1;
    end
    if coun == 2
        PLSR = 10* log10(abs(s_a(k + 1, Maxran)) / abs(s_a(max_aind, Maxran)));
        break
    end
end
%% Range 3dB  azimuth 3dB
% azimuth
for k = max_aind : length(s_a(:,1))
    if abs(s_a(max_aind, Maxran)) / 2  >= abs(s_a(k, Maxran))
        asandb = k - max_aind ;
        break
    end
end
% range 
[dum max_rind] = max(abs(s_a(max_aind,:)));
for k = max_rind : length(s_a(1,:))
    if abs(s_a(max_aind, max_rind)) / 2 >= abs(s_a(max_aind, k))
        rsandb = k - max_rind;
        break
    end
end

close all

end