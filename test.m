close all
clear
%Time space
time_mov = 1/8; 
tot_point = 1024*8;
t_1 = linspace(0, 1, tot_point);
t_2 = linspace(time_mov * 2, 1, tot_point * (1 - time_mov * 2)) ;
% a_m parameters
a_0 = 100 ;
a_1 = 10 * 2*pi;
a_2 = 100 * 2*pi;
a_3 = 500 * 2*pi;
a_4 = 000 * 2*pi;
% Siganl 
s_1 = exp(j*(a_0 + a_1 *t_1 + a_2 * t_1.^2 + a_3 *t_1.^3 + a_4 * t_1.^4)) ;
s_2 = [exp(-j*(a_0 + a_1 *t_2 + a_2 * t_2.^2 + a_3 *t_2.^3 + a_4 * t_2.^4)), zeros(1, tot_point * time_mov * 2)];
s1D = s_1 .* s_2 ;

s = s1D(1: tot_point * 0.5) .* conj(s1D(tot_point * 0.25+1: tot_point * 0.75)) ;
figure
plot(real(s))
figure

spectrogram(s1D,256,250,256,tot_point,'yaxis')
figure
%spectrogram(s,256,250,256,tot_point,'yaxis')
plot(linspace(0,tot_point,length(s)),abs(fft(s)))  
    set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
	xlabel('Hz', 'Interpreter', 'latex')
	

%{
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