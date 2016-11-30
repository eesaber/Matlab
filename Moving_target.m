%% Parameters - C Band airborne SAR parameters
% Target
x_n = 5000;
y_n = 0;
d = 200; % Length of the target area 
v_x = 0; a_x = 0; % range
v_y = 20; a_y = 0; % azimuth 
% Platform
h = 5000;
La = 2;     
v_p = 100; 
f0 = 8.85e9; c = 3e8 ; lambda = c/f0 ;
dur = 6 ; %
PRF = 1000; %%2.35e3
Kr = 20e12 ;
Tp = 1.5e-6; % Pulse width
B =  Kr*Tp;
Ta = 0.886 * (lambda * x_n)/(La *v_p *1);

%% Signal 
eta = linspace( -0-dur/2,-0+dur/2,PRF*dur) ;  % Slow time -6

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
for k = 1 : length(eta)
    Fs(k,:) = fft(s(k,:),2^M);
    Fs_r(k,:) = Fs(k,:) .* Fs_0 ;
    s_r(k,:) = ifft( Fs_r(k,:) );
end
%purinto(s_r); % Print the s_r
%% RCMC for range curveture 
% f_tau 
f_tau = 0  : 1 : 2^M - 1;
f_tau = f_tau/(2^M/4/B) ;
H_c =  exp(  j* (2* pi /c )* ( (20-v_p)^2 / R_0) *eta'.^2 * f_tau); %range curveture filter
%Multiply H_c to Fs_r
Fs_1 = Fs_r .* H_c ;
for k = 1 : length(eta)
    s_1(k,:) = ifft( Fs_1(k,:) );
end
purinto(s_1); % Print the s_1
%% Compute v_r 
control_rt = false;  % Control if the procedure go through hough transform 
if control_rt == true
		[Index temp] = max(abs(s_1') ); 
     % Apply Hough tranform to fint the angle.
end

%% RCMC for range walk
H_w = exp(j* 4* pi /c * v_x /1.41* eta' * f_tau);
Fs_2 = Fs_1 .* H_w;  
for k = 1 : length(eta)
   s_2(k,:) = ifft( Fs_2(k,:) );
end
purinto(s_2); % Print the s_2
%% Likelihood computation
control_lc = false ;
if control_lc == true
    [Maxnum , Ind] = max(abs(s_2(:)));
    [I_row, I_col] = ind2sub(size(s_2),Ind) ;
    step = 0.1;
    Test_v_y = -20: step :20 ;
    likevector = linspace(0,0,length(Test_v_y) ) ;
    for  k = 1 : length(Test_v_y) 
        K_a = 2/lambda *(Test_v_y(k) - v_p)^2 / R_0 ;
        s_a0 = exp(-j * 2* pi / lambda * K_a* eta.^2 ).* (eta <= (max(eta)+ min(eta))/2) ;
        likevector(k) = max(abs( conv(s_a0,s_2(:,I_col) )));
    end
   [Maxnum , Ind] = max(likevector) ;
    Sim_v_y = (Ind - length(Test_v_y) /2 ) * step
end
%plot(real(s_2(:,I_col)))
%plot(Test_v_y,likevector)
%% Azimuth compression 
control_ac = false ;
if control_ac == true 
	% Azimuth reference signal
	eta0 = linspace(-dur/2, dur/2,length(s_2(:,1)) );
	K_a = 2/lambda *(v_y - v_p)^2 / R_0 ;
	s_a0 = exp(j * pi* K_a* eta0.^2 ) ;
	M = nextpow2(length(eta0));
	Fs_a0 = fft(s_a0,2^M)'; 
	Fs_a0 = conj(Fs_a0); 

	for k = 1 : length(s_2(1,:))
		Fs_2a = fft(s_2(:,k),2^M); 
        Fs_a(:,k) = Fs_2a.* Fs_a0 ;
        s_a = ifft(Fs_a);
    end
    purinto(s_a)
    %clearvars eta0 M eta tau
end 
%%
figure
temp = spectrogram(s_2(:,293),128,120,300,1000,'centered','yaxis');
x = [-3 3];
y = [-500 500];
imagesc(x, y, abs(temp))
 set(gca,'Ydir','normal')
set(gca,'FontSize',32,'Fontname','CMU Serif')
xlabel('$\eta$','Interpreter','latex')
ylabel('$f_\eta$','Interpreter','latex')
%%
for i = 1 : 5
	X_a(i,:) = frft(s_2(:,293,0.752 + 0.01*i))
	plot(abs(X_a(i,:)),'LineWidth', 2.5)
	hold on
end 
%axis([2700, 3300, 0, 7000])
set(gca,'Ydir','normal')   
xlabel('$u$','Interpreter','latex')
ylabel('$\left| X_a(u) \right|$','Interpreter','latex')
set(gca,'FontSize',28,'Fontname','CMU Serif')

