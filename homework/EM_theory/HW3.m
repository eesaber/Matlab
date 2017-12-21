% EM theory HW4
eta = real(120*pi*[1 , 1/sqrt(4-1j)]) ;
t_in = linspace(0,pi/2,500);
s_t = 1/sqrt(4-1j)*sin(t_in);
c_t = sqrt(1-s_t.^2);
t_re = atan(sin(t_in)./real(sqrt((4-1j)^2-sin(t_in).^2)));

gamma_h = (eta(2)./c_t - eta(1)./cos(t_in))./(eta(2)./c_t + eta(1)./cos(t_in));
tau_h = 2*eta(2)./c_t./(eta(2)./c_t + eta(1)./cos(t_in));

gamma_v = (eta(2)*c_t - eta(1)*cos(t_in))./(eta(2)*c_t + eta(1)*cos(t_in));
tau_v = 2*eta(2)*cos(t_in)./(eta(2)*c_t + eta(1)*cos(t_in));

AR_r = 20*log10(abs( (abs(gamma_v+gamma_h) + abs(gamma_v-gamma_h))./(abs(gamma_v + gamma_h) - abs(gamma_v - gamma_h)) ));
AR_t = 20*log10(abs( (abs(tau_v - tau_h) + abs(tau_v + tau_h))./(abs(tau_v - tau_h) - abs(tau_v + tau_h)) ));

figure(1)
plot(t_in, abs(gamma_h),'k',t_in, abs(gamma_v),'k-.','Linewidth',3)
	xlim([0 pi/2])
	set(gca,'xtick',linspace(0,pi/2,4),'XTickLabel', linspace(0,90,4));
	xlabel('$\theta^\circ$','Interpreter', 'latex')
	ylabel('$|\Gamma(\theta)|$','Interpreter', 'latex')
	legend('|\Gamma_h(\theta)|','|\Gamma_v(\theta)|','Location','northwest')
	grid on 
	plot_para('Maximize',true,'Filename','hw4_1')
	
figure(2)
	plot(t_in, abs(tau_h),'k',t_in, abs(tau_v),'k-.','Linewidth',3)
	xlim([0 pi/2])
	set(gca,'xtick',linspace(0,pi/2,4),'XTickLabel', linspace(0,90,4));
	xlabel('$\theta^\circ$','Interpreter', 'latex')
	ylabel('$|\tau(\theta)|$','Interpreter', 'latex')
	legend('|\tau_h(\theta)|','|\tau_v(\theta)|','Location','northeast')
	grid on 
	plot_para('Maximize',true,'Filename','hw4_2')
	
figure(3)
	plot(t_in,AR_r,'k','Linewidth',3)
	hold on 
	plot(t_in, (abs(gamma_v +gamma_h)>abs(gamma_h-gamma_v))*80,'b','Linewidth',1)
	hold off
	xlim([0 pi/2])
	set(gca,'xtick',linspace(0,pi/2,4),'XTickLabel', linspace(0,90,4));
	xlabel('$\theta^\circ$','Interpreter', 'latex')
	ylabel('AR (dB)','Interpreter', 'latex')
	ylim([0 75])
	grid on 
	plot_para('Maximize',true,'Filename','hw4_3')
figure(4)
	plot(t_in,AR_t,'k', t_in,(abs(tau_v-tau_h)>abs(tau_h + tau_v))*10, 'k-.','Linewidth',3)
	xlim([0 pi/2])
	set(gca,'xtick',linspace(0,pi/2,4),'XTickLabel', linspace(0,90,4));
	xlabel('$\theta^\circ$','Interpreter', 'latex')
	ylabel('AR (dB)','Interpreter', 'latex')
	grid on 
	plot_para('Maximize',true,'Filename','hw4_4')