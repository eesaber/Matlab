function [s] =  Gen_signal(vx, vy)
    clear
    delete parameter.mat
    if nargin ~= 2
        vx = 0; vy = vx; 
    end
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
    save('parameter.mat'); %Save the parameters for the after processing.
     
    for k = 1 : length(eta)
        td = 2* R(k)/ c ;
        s(k,:) = exp(-i* 4* pi*R(k)/ lambda + i* pi* Kr* (tau - 2* R(k)/ c).^2) .*(tau>=td & tau-td<=Tp)  ;
    end
end