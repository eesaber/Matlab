close all
theta = linspace(0, 2*pi, 1000);
%% HW 5-8-4
b = 5*pi/4*sin(theta);
E = 1 + exp(j*b) + exp(j*2*b) + exp(j*3*b);
figure
    polarplot(theta, abs(E)/max(abs(E)), 'Linewidth', 3)
    set(gca,'GridAlpha',0.6)
    plot_para('Maximize',true,'Filename','hw5_584a')

E = 1 + 3*exp(j*b) + 3*exp(j*2*b) + exp(j*3*b);
figure
    polarplot(theta, abs(E)/max(abs(E)), 'Linewidth', 3)
    set(gca,'GridAlpha',0.6)
    plot_para('Maximize',true,'Filename','hw5_584c')
%% HW 5-8-6
close all
clear all
theta = linspace(0, 2*pi, 2000);
b =  pi*sin(theta);
E = 0;
N = 5;
for k = 0 : N
    E = E + nchoosek(N,k)*exp(1j*k*b);
    %E = E + exp(1j*k*b);
end

D = 2/sum((abs(E(1:numel(theta)/2))/abs(max(E(1:numel(theta)/2)))).^2.*sin(theta(1:numel(theta)/2))*(theta(2)-theta(1)));
figure
    polarplot(theta, abs(E)/max(abs(E)), 'Linewidth', 3)
    set(gca,'GridAlpha',0.6)
    title(['$D = $' num2str(D)],'Interpreter', 'latex')
    plot_para('Maximize',true,'Filename','hw5_586')

%% HW 5-17-1
close all
E = 0;
theta = linspace(0.001,pi,1000);
E = 1/24*sin(12*pi*cos(theta)-6*pi)./sin(pi/2*cos(theta)-pi/2);
b_a = 2*pi*sum(abs(E).^2.*sin(theta)*(theta(2)-theta(1)));
b_M = 2*pi*sum(abs(E(theta<23/180*pi)).^2.*sin(theta(theta<23/180*pi))*(theta(2)-theta(1)));

fprintf('Beam area: %f, Main beam area: %f \n', b_a, b_M)
figure 
plot(theta*180/pi,abs(E))
figure 
polarplot(theta,abs(E))

%% 5-18-1
theta = linspace(0,2*pi,10000);
b = 4*pi*cos(theta);
E = 0;
E = 1 + exp(j*b);
figure 
    plot(theta*180/pi,abs(E)/max(abs(E)), 'Linewidth', 3)
    set(gca,'GridAlpha',0.6)
    plot_para('Maximize',true)

%% HW 5 extra

theta = linspace(-pi/2, pi/2, 2000);
b =  pi*sin(theta)-90.1/180*pi;
E = 0;
N = 10;
a = [1.0000    1.3555    1.9679    2.4787    2.7695    2.7695    2.4787    1.9679    1.3555    1.0000];
for k = 0 : N-1
    E = E + a(k+1)*exp(1j*k*b);
    %E = E + exp(1j*k*b);
end
figure
plot(theta*180/pi,20*log10(abs(E)/abs(max(E))), 'Linewidth', 3)
xlim([-90 90])
ylim([-40 0])
ylabel('$\frac{|E|}{\max |E|}$ (dB)','Interpreter', 'latex')
xlabel('$\theta$ (deg)','Interpreter', 'latex')
grid on
set(gca,'GridAlpha',0.6)
plot_para('Maximize',true,'Filename','hw5_extra')