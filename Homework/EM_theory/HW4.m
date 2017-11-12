% EM theory HW4
theta = linspace(0,pi/2-0.00001,91);
eta = 120*pi;
eplison = [1, 2, 4-1j, 2, 1];
z = eta./sqrt(eplison);
l = [1, 1/4, 1/3, 1/4, 1]; 
co = ones(5,numel(theta));
si = co;
be = co;

si(1,:) = sin(theta);
co(1,:) = cos(theta);
for k = 2 : 5
	si(k,:) = sqrt(eplison(k-1)/eplison(k))*si(k-1,:);
	co(k,:) = sqrt(1-si(k,:).^2);
	be(k,:) = 2*pi*l(k)*sqrt(eplison(k))*co(k,:);
end
gamma_h = zeros(1, numel(theta));
gamma_v = gamma_h;
tau_h = gamma_h; 
tau_v = gamma_h;
for k = 1 : numel(theta)
	T = [cos(be(2,k)), 1j*z(2)/co(2,k)*sin(be(2,k)); 1j*sin(be(2,k))/z(2)*co(2,k), cos(be(2,k))]......
		*[cos(be(3,k)), 1j*z(3)/co(3,k)*sin(be(3,k)); 1j*sin(be(3,k))/z(3)*co(3,k), cos(be(3,k))]......
		*[cos(be(4,k)), 1j*z(4)/co(4,k)*sin(be(4,k)); 1j*sin(be(4,k))/z(4)*co(4,k), cos(be(4,k))];
	gamma_h(k) = (T(1,1) + T(1,2)/z(5)*co(5,k) - z(1)/co(1,k)*(T(2,1) + T(2,2)/z(5)*co(5,k)))/(T(1,1) + T(1,2)/z(5)*co(5,k) + z(1)/co(1,k)*(T(2,1) + T(2,2)/z(5)*co(5,k)));
	tau_h(k) = 2/(T(1,1) + T(1,2)/z(5)*co(5,k) + z(1)/co(1,k)*(T(2,1) + T(2,2)/z(5)*co(5,k)));
	
	T = [cos(be(2,k)), 1j*z(2)*co(2,k)*sin(be(2,k)); 1j*sin(be(2,k))/z(2)/co(2,k), cos(be(2,k))]......
		*[cos(be(3,k)), 1j*z(3)*co(3,k)*sin(be(3,k)); 1j*sin(be(3,k))/z(3)/co(3,k), cos(be(3,k))]......
		*[cos(be(4,k)), 1j*z(4)*co(4,k)*sin(be(4,k)); 1j*sin(be(4,k))/z(4)/co(4,k), cos(be(4,k))];
	gamma_v(k) = (T(2,1)*z(5)*co(5,k) + T(2,2) - (T(1,1)*z(5)*co(5,k) + T(1,2))/z(1)/co(1,k))/(T(2,1)*z(5)*co(5,k) + T(2,2) + (T(1,1)*z(5)*co(5,k) + T(1,2))/z(1)/co(5,k));
	tau_v(k) = 2/(T(2,1)*z(5)*co(5,k) + T(2,2) + (T(1,1)*z(5)*co(5,k) + T(1,2))/z(1)/co(5,k));
end




figure(1)
plot(theta, abs(gamma_h),'k', theta, abs(gamma_v),'k-.','Linewidth',3)
	xlim([0 pi/2])
	set(gca,'xtick',linspace(0,pi/2,6),'XTickLabel', linspace(0,90,6));
	xlabel('$\theta^\circ_i$','Interpreter', 'latex')
	ylabel('$|\Gamma(\theta)|$','Interpreter', 'latex')
	legend('|\Gamma_h(\theta)|','|\Gamma_v(\theta)|','Location','west')
	grid on 
	plot_para('Maximize',true,'Filename','hw5_1')
figure(2)
yyaxis left
	plot(theta, abs(tau_h),theta, abs(tau_v),'-.','Linewidth',3)
	xlim([0 pi/2])
	set(gca,'xtick',linspace(0,pi/2,6),'XTickLabel', linspace(0,90,6));
	xlabel('$\theta^\circ_i$','Interpreter', 'latex')
	ylabel('$|\tau(\theta)|$','Interpreter', 'latex')
	legend('|\tau_h(\theta)|','|\tau_v(\theta)|','Location','west')
yyaxis right
	plot(theta, asind(si(5,:)),'Linewidth',3)
	ylabel('$\theta_t$','Interpreter', 'latex')
	ylim([0 90])
	grid on 
	plot_para('Maximize',true,'Filename','hw5_2')
figure(3)
	plot(theta,abs(gamma_h).^2 + abs(tau_h).^2*z(1)./co(1,:)/z(5).*co(5,:),'k','Linewidth',3)
	hold on 
	plot(theta,abs(gamma_v).^2 + abs(tau_v).^2*z(5).*co(5,:)/z(5)./co(5,:),'k-.','Linewidth',3)
	xlim([0 pi/2])
	set(gca,'xtick',linspace(0,pi/2,6),'XTickLabel', linspace(0,90,6));
	%ylim([0.99 1.01])
	xlabel('$\theta^\circ_i$','Interpreter', 'latex')
	ylabel('Remain power','Interpreter', 'latex')
	legend('h-pol','v-pol','Location','northwest')
	plot_para('Maximize',true,'Filename','hw5_3')