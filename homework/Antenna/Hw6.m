%% 6-5-1
t = linspace(0,2*pi, 720);
E = abs(sin(t)./cos(t).*sin(pi/2*cos(t)));
figure
polarplot(t,E/max(E),'Linewidth',3)
set(gca,'GridAlpha',0.6)
title('$|E| / |E_{max}|$','Interpreter', 'latex')
plot_para('Maximize',true,'Filename','hw5_651')

%% 6-7-1
E = abs(cos(pi/2*sqrt(2)*sin(t))./sin(t-pi/4).*cos(pi/2*cos(t-pi/4)));
figure
polarplot(t,E/max(E),'Linewidth',3)
set(gca,'GridAlpha',0.6)
title('$|E| / |E_{max}|$','Interpreter', 'latex')
plot_para('Maximize',true,'Filename','hw5_671')

%% 7-8-2
E = abs(sin(2*pi/6*sin(t)));
E = E/max(E);
figure
    polarplot(t,E,'Linewidth',3)
    set(gca,'GridAlpha',0.6)
    title('$|E| / |E_{max}|$','Interpreter', 'latex')
    plot_para('Maximize',true,'Filename','hw5_7821')
E = abs(sin(2*pi*3/8*sin(t)));
figure
    polarplot(t,E/max(E),'Linewidth',3)
    set(gca,'GridAlpha',0.6)
    title('$|E| / |E_{max}|$','Interpreter', 'latex')
    plot_para('Maximize',true,'Filename','hw5_7822')
E = abs(sin(2*pi/1*sin(t)));
figure
    polarplot(t,E/max(E),'Linewidth',3)
    set(gca,'GridAlpha',0.6)
    title('$|E| / |E_{max}|$','Interpreter', 'latex')
    plot_para('Maximize',true,'Filename','hw5_7823')
    
%% 7-9-1
E = abs(sin(pi*cos(t)));
figure
    polarplot(t,E/max(E),'Linewidth',3)
    set(gca,'GridAlpha',0.6)
    title('$|E| / |E_{max}|$','Interpreter', 'latex')
    plot_para('Maximize',true,'Filename','hw5_791')

%% Bonus
close all
t = linspace(0,2*pi,360);
E = abs(2*cos(sqrt(2)*pi/2*cos(t)) - 2*cos(sqrt(2)*pi/2*sin(t)));
E = E.*((t < pi/4) + (t > 7*pi/4));
figure
    polarplot(t,E/max(E),'Linewidth',3)
    set(gca,'GridAlpha',0.6)
    title('$|E| / |E_{max}|$ H-plane','Interpreter', 'latex')
    thetaticks(0:45:315)
    plot_para('Maximize',true,'Filename','H')
    

array_reflect = (-2 + 2*cos(sqrt(2)*pi/2*cos(pi/2-t)));
dipole = (cos(pi*cos(t)/2)-cos(pi/2))./sin(t);
array_phase = 1/4*sin(4*5/2*pi/2*cos(t))./sin(5/2*pi/2*cos(t));
E = abs(array_reflect.*dipole.*array_phase);
E = E.*((t < pi));
figure
    polarplot(t,E/max(E),'Linewidth',3)
    set(gca,'GridAlpha',0.6)
    title('$|E| / |E_{max}|$ E-plane','Interpreter', 'latex')
    plot_para('Maximize',true,'Filename','E')
    

