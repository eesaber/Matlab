clc
clear 
close all

%% Antenna Hw4 secret problem

freq = 75e6;
c = physconst('lightspeed');
lambda = c/freq;
N = 2;
dx= 2*lambda;

d = dipole('Length',1*lambda/10,'width',lambda/500);

dipole_array = linearArray;
dipole_array.Element = d;
dipole_array.NumElements = N;
dipole_array.ElementSpacing = dx;
hArray = figure;
show(dipole_array)
plot_para('Maximize',true,'Filename','hw4_sc1')
%axis tight

figure
pattern(dipole_array,freq)
plot_para('Maximize',true,'Filename','hw4_sc2')


%% Antenna Hw4 5-6-6

phi = linspace(0,2*pi,1000);
dis = 1/4;
E = 0;
for k = 0:11
	E = E + exp(j*(k*2*pi*dis*cos(phi)-k*2*pi*dis));
	%E = E + exp(j*(k*2*pi*dis*cos(phi)));
end
D = 2/sum((abs(E(1:numel(phi)/2))/abs(max(E(1:numel(phi)/2)))).^2.*sin(phi(1:numel(phi)/2))*(phi(2)-phi(1)));
fprintf('Directivity: %f\n', D)
figure
	polarplot(phi,abs(E)/abs(max(E)),'Linewidth',3) 
	set(gca,'GridAlpha',0.6)
	%title(['$D = $' num2str(D) ' or 95.534'],'Interpreter', 'latex')
	title(['$D = $' num2str(D)],'Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','hw4_565a')
D_guji = 28000/(2*(90-81.44))^2;

%% Antenna Hw4 5-6-4

E = 0;
for k = 0:7
	E = E + exp(j*(k*0.4*pi*cos(phi)-k*0.4*pi));
end
D = 2/sum((abs(E(1:numel(phi)/2))/abs(max(E(1:numel(phi)/2)))).^2.*sin(phi(1:numel(phi)/2))*(phi(2)-phi(1)));
fprintf('Directivity: %f\n', D)
figure
	polarplot(phi,abs(E)/abs(max(E)),'Linewidth',3)
	set(gca,'GridAlpha',0.6)
	title(['$D = $' num2str(D)],'Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','hw4_564a')
E = 0;
for k = 0:7
	E = E + exp(j*(k*0.4*pi*cos(phi)-k*0.4*pi - k*pi/8));
end
%D = 2/sum((abs(E)/abs(max(E)).^2)*(phi(2)-phi(1)));
D = 2/sum((abs(E(1:numel(phi)/2))/abs(max(E(1:numel(phi)/2)))).^2.*sin(phi(1:numel(phi)/2))*(phi(2)-phi(1)));
fprintf('Directivity: %f\n', D)
figure
	polarplot(phi,abs(E)/abs(max(E)),'Linewidth',3)
	set(gca,'GridAlpha',0.6)
	title(['$D = $' num2str(D)],'Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','hw4_564b')

%% Antenna Hw4 5-4-1

phi = linspace(0,2*pi,1000);
phase = pi/180*[121, -128];

E_1 = cos(pi*0.84*cos(phi) + phase(1));
E_2 = cos(pi*0.3*cos(phi) + phase(2));
figure(1)
polarplot(phi,abs(E_1),'Linewidth',3)
thetaticks(0:45:315)
set(gca,'GridAlpha',0.6)
plot_para('Maximize',true,'Filename','hw4_541_1')
figure(2)
polarplot(phi,abs(E_2),'Linewidth',3)
set(gca,'GridAlpha',0.6)
thetaticks(0:45:315)
plot_para('Maximize',true,'Filename','hw4_541_2')
figure(3)
polarplot(phi,abs(E_1.*E_2),'Linewidth',3)
set(gca,'GridAlpha',0.6)
thetaticks(0:45:315)
plot_para('Maximize',true,'Filename','hw4_541_s')


%% 5-3-1	

phi = linspace(0,2*pi,1000);
d = [1/16, 1/8, 1/4, 1];

for k = 1 : 4
	
E = cos(2*pi*d(k)*cos(phi));
P = (E>0)*0+(E<0)*-180;
figure
	subplot(2,1,1)
	plot(phi/pi*180,abs(E),'Linewidth',3)
	ylabel('$|E|$','Interpreter', 'latex')
	xlabel('$\theta^\circ$','Interpreter', 'latex')
	xlim([0,360])
	grid on 
	set(gca,'GridAlpha',0.6)
	plot_para()
	
	subplot(2,1,2)
	plot(phi/pi*180,P,'Linewidth',3)
	ylabel('Phase angle','Interpreter', 'latex')
	xlabel('$\theta^\circ$','Interpreter', 'latex')
	xlim([0,360])
	grid on
	set(gca,'GridAlpha',0.6)
	plot_para('Maximize',true,'Filename',['hw4_531_', num2str(k)])
end

%% 5-2-1
clear
close all
phi = linspace(0,2*pi,1000);
d = [1, 3/2, 4, 1/4];

for k = 1 : 4
E = cos(2*pi*d(k)/2*sin(phi));
figure
	polarplot(phi,E,'Linewidth',3)
	set(gca,'GridAlpha',0.6)
	plot_para('Maximize',true,'Filename',['hw4_521_', num2str(k)])
end

