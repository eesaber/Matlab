%% HAF kernal of rect(t) 
s = ones(1,1000); % one point 0.001 s 
N = length(s) ;
s = [zeros(1,500) s zeros(1,500)];
t = linspace(-1, 1, length(s));

% K_1
K_1 = s ;

% K_2 
K_2 = zeros(N,length(s)); 
for xi = 1 : N / 2 
    K_2(xi + N/2,:) = [zeros(1,500 + xi) ones(1,1000 - 2 *xi) zeros(1,500 + xi)];
    K_2(1 -xi + N/2,:) = [zeros(1,500 + xi) ones(1,1000 - 2 *xi) zeros(1,500 + xi)];
end
%AF_2 = fft(K_2,2^next2pow)
% K_3 
K_3  = zeros(N, N, length(s));
for xi = 1 : N / 2
    K_3(xi + N / 2, :, :) = [K_2(:,xi:length(s)) zeros(N,xi-1)] .* [zeros(N,xi-1) K_2(:, xi:length(s)) ];
    K_3(1 - xi + N/2, :, :) = [K_2(:,xi:length(s)) zeros(N,xi-1)] .* [zeros(N,xi-1) K_2(:, xi:length(s))];
end 
close all
figure 
    plot(t, K_1,'k','linewidth',4)
figure
    imagesc(t, linspace(-1/2,1/2,1000), K_2)
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
    
    set(gcf,'color','w');
	set(gca,'Ydir','normal')   
	xlabel('t', 'Interpreter', 'latex')
	ylabel('$\xi_1$', 'Interpreter', 'latex')
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
	colormap('Jet')
	colorbar
    %export_fig K_2.jpg
 figure 
    surf(K_3(:,:,1000))

%% High-order Ambiguity function 
close all
clear
%Time space
t_mov = 1/6;
tot_point = 1024*6;
t_1 = linspace(0, 1, tot_point);
% a_m parameters
a_0 = 10 ;
a_1 = 100 ;
a_2 = 500 ;
a_3 = 200 ;
a_4 = 000 ;
% Siganl 
s = exp(j * 2 * pi * (a_0 + a_1 *t_1 + a_2 * t_1.^2 / 2 + a_3 *t_1.^3 / 6+ a_4 * t_1.^4 / 24)) ;
AF = GAF(s,3,3);
[~,qq] = max(abs(AF));
f = linspace(-tot_point/2, tot_point/2, 2^(nextpow2(length(s)-1)) );
xi = 1 / 3 / 2 ;
t_a_3 = f(qq) / 4 / xi^2;

%{
figure
    plot(f,abs(AF_2),'k','Linewidth',3.5)
    set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
    set(gcf,'color','w');
    xlabel('Hz', 'Interpreter', 'latex')
    xlim([-1000 0])
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
    %export_fig AF_2.jpg
figure
    plot(f,abs(AF_3),'k','Linewidth',3.5)
    set(gcf,'color','w');
    set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
    xlabel('Hz', 'Interpreter', 'latex')
    xlim([-100 100])
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
    %export_fig AF_3.jpg
	



%% Expand the High order Generalized Ambiguity function 
clc
clear
syms x xi_1 temp
for m = 1 : 5
    
    expand( (x + xi_1)^m )/ factorial(m) - expand( (x - xi_1)^m )/ factorial(m)
end

%% 
x_n = 5000;
y_n = 0;
d = 200; % Length of the target area 
v_x = -10; a_x = 0; % range
v_y = 0; a_y = 0; % azimuth 
% Platform
h = 5000;
La = 2;     
v_p = 100; 

eta = -50 : 0.1 : 50 ;
R1 = sqrt(h^2 + (x_n + 20* eta + 0.5* a_x* eta.^2).^2 + (y_n + -20* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ;
R2 = sqrt(h^2 + (x_n + 0* eta + 0.5* a_x* eta.^2).^2 + (y_n + 20* eta + 0.5* a_y * eta.^2 - v_p* eta).^2 ) ;
plot(R1)
hold
plot(R2)
%}