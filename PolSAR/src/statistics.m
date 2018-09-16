%% Distribution
clear Y
x = linspace(0,5,5000);
a = [0.001, 0.003, 0.007, 0.01, 0.03, 0.07, 0.1, 0.2, 0.3, 0.7, 1,3];
b = 2;
type = 5; 
for n = 1 : numel(a)
	switch type
		case 1
			Y(n,:) = pdf('Beta',x, a(n), b);
		case 2
			Y(n,:) = pdf('Rayleigh',x, a(n), b);
		case 3
			Y(n,:) = pdf('InverseGaussian',x, a(n), b);
		case 4
			Y(n,:) = pdf('Gamma',x, a(n), b);
        case 5
            Y(n,:) = pdf('Chisquare',x, n);
	end
end
figure
plot(x,Y,'Linewidth',2)
legend(num2str(a'))
grid on
ylim([0 10])
plot_para('Maximize',true)
%% Log-cumulants diagram
x = linspace(0,5,2000);
plot(psi(3,x), -psi(2,x), 'k', -psi(3,x), -psi(2,x), 'k--', 'Linewidth',2.5)
set(gca, 'Xtick', -5:5,'xlim',[-5 5])
text(-2,0.4,'A','Fontsize',40, 'Fontname','CMU serif')
text(2,0.4,'B','Fontsize',40, 'Fontname','CMU serif')
text(-0.2,1,'C','Fontsize',40, 'Fontname','CMU serif')
xlabel('$\kappa_3$','interpreter','latex')
ylabel('$\kappa_2$','interpreter','latex')
plot_para('Maximize',true,'Filename','log_cumulant_diagram','Ratio',[4 3 1])