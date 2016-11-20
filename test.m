a = linspace(179/180*pi, 179.5/180*pi , 100);
figure
plot(1/sin(a))
hold
plot(1/(pi-a))
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