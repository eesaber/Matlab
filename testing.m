%% DATA generate
f1 = 9.91e10;
f2 = 8.8e10;
f3 = 1.27e10;
%f_ = -sort(-[f1, f2, f3]) ;
f_ = -sort(-[f1, f2, f1-f2, f1+f2]);
lamda_ = c./f_;
N = -300;
space = 2*pi*(0:N/abs(N):N);
clear k k_bar temp 
temp = mod(4*pi./lamda_*(d_a*y_0/R_0), 2*pi) .* lamda_;
k = repmat(temp.',1,length(space)) + lamda_.' * space ;
k_bar = [k(:,1),  diff(k,1,2)];

%% Plot the 
x_l = 3;
figure 
subplot('position',[0.1 0.5 0.8 0.15]);
%subplot('position',[0.1 0.3375 0.8 0.1125]);
barh(k_bar,1,'stack','FaceColor',[255,197,108]/256,'Linewidth',1.25,'EdgeColor',[204, 158, 86]/256, 'ShowBaseLine', 'off');
set(gca,'TickLength',[0 0])
ylim([0.5 numel(k)/length(k)+0.5])
hold all
if y_0 >= 0 
	annotation('arrow', (0.1+ 0.8*(4*pi*d_a*y_0/R_0 / abs(x_l - 0)))*ones(1,2), [.66 .65],'Color',[214,57,50]/256)
else 
	annotation('arrow', (0.9+ 0.8*(4*pi*d_a*y_0/R_0 / abs(x_l - 0)))*ones(1,2), [.66 .65],'Color',[214,57,50]/256)
end
%plot(4*pi*d_a*y_0/R_0*[1 1],[.5 4.5],'Color', [214,57,50]/256, 'Linewidth',3)
x_l = 2*x_l*(N>0) - x_l;
xlim(sort([0,x_l]))
if N > 0
	xlabel(['$\lambda_\alpha b_{am} + 2 N \lambda_\alpha \pi$, $0 \leq N \leq $' int2str(N)],'Interpreter', 'latex')
else 
	xlabel(['$\lambda_\alpha b_{am} + 2 N \lambda_\alpha \pi$, ', int2str(N) ' $\leq N \leq 0$'],'Interpreter', 'latex')
end


[temp, ind_] = sort(reshape(k,1,numel(k)));
row_ = reshape( (1:numel(k)/length(k)).' *ones(1,length(space)), 1, numel(k)); 
row_ = row_(ind_);
for i = 1 : numel(k)- numel(k)/length(k)
	if temp(i+numel(k)/length(k)-1) - temp(i) < 0.01 && sum(diff(sort(row_(i: i + numel(k)/length(k) - 1))) == 0) == 0
		for j = 0 : 3
			%animatedline([temp(i+j) temp(i+j)],[row_(i+j)-0.5 row_(i+j)+0.5],'Color',[ 2,65,226]/256,'Linewidth',2); 
			%drawnow
			plot([temp(i+j) temp(i+j)],[row_(i+j)-0.5 row_(i+j)+0.5],'Color',[ 2,65,226]/256,'Linewidth',2); 
			
		end
	end
end

plot_para('Maximize',true,'Filename','1')


%{
track = 3;
fprintf ('Real phase: %f\n', 4*pi/lamda_(track)*(d_a*y_0/R_0) )
fprintf ('F band: %i, Amb number: %i, Real phase: %f\n', track, 2, mod(4*pi/lamda_(track)*(d_a*y_0/R_0), 2*pi) + 2 * 2*pi )
fprintf ('For each interval, y_0 will differ %f\n', 10)

figure
plot(k(1,:), ones(length(space)),'ko', 'linewidth',1.5 )
hold on 
plot(k(2,:), 1.25*ones(length(space)),'bo','linewidth',1.5 )
%ylim([0.75 1.5])
%xlim([-0.45 0])
%xlim([0 0.45])
xlabel('$b_{am} + 2 N \pi$, $0 \leq N \leq 20$ ','Interpreter', 'latex')
grid on
grid minor
%plot_para('Filename','1','Maximize',true)
%}

%{
xi = -2:0.01:2;
real_12 = 4*pi*f_0/c*(d_a*y_0/1.5/R_0 - d_a*(v_p-v_y)/2/R_0*xi);
real_13 = 4*pi*f_0/c*(d_a*y_0/R_0 - d_a*(v_p-v_y)/R_0*xi);
close all
figure
plot(xi,real_12/pi - 2,xi,real_13/pi - 4)
figure
plot(xi,unwrap(angle(exp(j*real_12)))/pi,xi,unwrap(angle(exp(j*real_13)))/pi)
figure
plot(xi,real_12/pi,xi,real_13/pi)


x = 0: 0.01:1;
y = 1 - x;
ref = 1 - x.^2 - y.^2;
b = x.*(1 - x + y).^2 + y.*(-1-x+y).^2;
c = -x .* log(x) - y.*log(y);


plot(ref/max(ref))
figure
plot(ref/max(ref) - b/max(b))
figure
plot(ref/max(ref) - c/max(c))


syms t x_0 y_0 v_x v_y v_z v_p a_x a_y h 
	f = sqrt( (x_0 + v_x * t + a_x* t^2 / 2)^2 + h^2 + ( y_0 + (v_y - v_p) * t + a_y * t^2 / 2 )^2  );
%f = (y_0 + (v_y - v_p) * t + a_y * t^2 / 2) / sqrt( (x_0 + v_x * t + a_x* t^2 / 2)^2 + h^2 ...
%		+ ( y_0 + (v_y - v_p) * t + a_y * t^2 / 2 )^2  );

	taylor(f, t,'ExpansionPoint',0,'order', 2)
	taylor(f, t,'ExpansionPoint',0,'order', 3)
	taylor(f, t,'ExpansionPoint',0,'order', 4)



(h^2 + x_0^2 + y_0^2)^(1/2) 
% t 
((2*v_x*x_0 - 2*y_0*(v_p - v_y)))/(2*(h^2 + x_0^2 + y_0^2)^(1/2))
% t^2
- (t^2*((2*v_x*x_0 - 2*y_0*(v_p - v_y))^2/(4*(h^2 + x_0^2 + y_0^2)^(1/2)) ...
	- (h^2 + x_0^2 + y_0^2)^(1/2)*((v_p - v_y)^2 + a_x*x_0 + a_y*y_0 + v_x^2)))/(2*(h^2 + x_0^2 + y_0^2))
% t^3	
(t^3*((3*(a_x*v_x - a_y*(v_p - v_y))*(h^2 + x_0^2 + y_0^2)^(1/2))/2 + ...
	(3*(2*v_x*x_0 - 2*y_0*(v_p - v_y))*((2*v_x*x_0 - 2*y_0*(v_p - v_y))^2/(4*(h^2 + x_0^2 + y_0^2)^(1/2)) ...
	- (h^2 + x_0^2 + y_0^2)^(1/2)*((v_p - v_y)^2 + a_x*x_0 + a_y*y_0 + v_x^2)))/(4*(h^2 + x_0^2 + y_0^2))))/(3*(h^2 + x_0^2 + y_0^2))
 
%}
%{
n = 0:160;
la = 0.98;
-en = sum(exp(-la) * la.^n ./  factorial(n) .* log((exp(-la) * la.^n ./  factorial(n))))
10000*en/log(2)
%}