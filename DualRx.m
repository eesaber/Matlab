function DualRx()
% This function generate the SAR signal regarding to the input parameter.
% Usage: Gen_signal(vx, vy, ax, ay), 'vx' is range velocity, 'vy' is azimuth
% velocity, 'ax' is range accelaration and 'ay' is azimuth accelration. If no input parameters, the velocity in both direction are set
% to zero
% Noth that the parameters of SAR need to be modified in the function.

    delete parameter.mat
    if nargin == 0
        vx = 0; vy = 0; 
        ax = 0; ay = 0;
    end
    %% Parameters - C Band airborne SAR parameters
    % Target
    x_0 = 2160;
    y_0 = 400;
	
    d = 300; % Length of the target area 
    v_x = 0; a_x = 0; % rangecl
    v_y = 0; a_y = 0; % azimuth 
    % Platform
    h = 2200;
	R_0 = sqrt(x_0^2 + y_0^2 + h^2);
    d_a = 0.2;     
    v_p = 90; 
    f_0 = 9.6e9; c = 3e8 ; lambda = c/f_0 ;
    dur = 2 ; %
    PRF = 1000; %%2.35e3
    K_r = 5e15 ;
    T_p = 0.1e-6; % Pulse width
    B =  K_r*T_p;
   

    %% Signal 
    aa = 0;
    eta = linspace( aa - dur/2, aa + dur/2,PRF*dur) ;  % Slow time -6

    upran = sqrt( (x_0 + d)^2 + h^2) ; downran =  sqrt( (x_0 - d)^2 + h^2); %Set upper limit and down limit 
    tau =  2*(downran)/c: 1/4/B : 2*(upran)/c + T_p  ;  % fast time space

    R1 = sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
    R2 = sqrt(h^2 + (x_0 + v_x* eta + 0.5* a_x* eta.^2).^2 + (y_0 + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
   
    s1 = zeros(length(eta), length(tau));
	s2 = zeros(length(eta), length(tau));
    for k = 1 : length(eta)
        td1 = 2* R1(k)/ c ;
		td2 = (R1(k) + R2(k)) /c ;
        s1(k,:) = exp(-i* 4* pi*R1(k)/ lambda + i* pi* K_r* (tau - 2* R1(k)/ c).^2) .*(tau>=td1 & tau-td1<=T_p)  ;
		s2(k,:) = exp(-i* 2* pi*(R1(k) + R2(k))/ lambda + i* pi* K_r* (tau - (R1(k) + R2(k) )/ c).^2) .*(tau>=td2 & tau-td2<=T_p)  ;
	end
	abs()
end
