function [s] =  Gen_signal(vx, vy, ax, ay)
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
    x_n = 1000/1.41;
    y_n = 0;
	
    d = 50; % Length of the target area 
    v_x = vx; a_x = ax; % rangecl
    v_y = vy; a_y = ay; % azimuth 
    % Platform
    h = 1000/1.41;
	R_0 = sqrt(x_n^2 + y_n^2 + h^2);
    La = 2;     
    v_p = 100; 
    f_0 = 5e9; c = 3e8 ; lambda = c/f_0 ;
    dur = 2 ; %
    PRF = 2000; %%2.35e3
    K_r = 5e15 ;
    T_p = 0.1e-6; % Pulse width
    B =  K_r*T_p;
    T_a = 0.886 * (lambda * x_n)/(La *v_p *1);

    %% Signal 
    aa = 0;
    eta = linspace( aa - dur/2, aa + dur/2,PRF*dur) ;  % Slow time -6

    upran = sqrt( (x_n + d)^2 + h^2) ; downran =  sqrt( (x_n - d)^2 + h^2); %Set upper limit and down limit 
    tau =  2*(downran)/c: 1/4/B : 2*(upran)/c + T_p  ;  % fast time space

    R = sqrt(h^2 + (x_n + v_x* eta + 0.5* a_x* eta.^2).^2 + (y_n + v_y* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ; % range equation
    %R_0, ind_R_0] = min(R);
    save('parameter.mat'); %Save the parameters for the after processing.
    s = zeros(length(eta), length(tau));
    for k = 1 : length(eta)
        td = 2* R(k)/ c ;
        s(k,:) = exp(-i* 4* pi*R(k)/ lambda + i* pi* K_r* (tau - 2* R(k)/ c).^2) .*(tau>=td & tau-td<=T_p)  ;
    end
end