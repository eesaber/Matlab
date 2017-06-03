f1 = 1.57e10;
f2 = 1.23e10;
f3 = 1.3e10;
%f_ = -sort(-[f1, f2, f3]) ;
f_ = -sort(-[f1, f2, f1-f2, f1+f2]);
lamda_ = c./f_;
N = -100;
space = 2*pi*(0:N/abs(N):N);
clear k k_bar temp b 
temp = mod(4*pi./lamda_*(d_a*y_0/R_0), 2*pi) .* lamda_;
k = repmat(temp.',1,length(space)) + lamda_.' * space ;
k_bar = [k(:,1),  diff(k,1,2)];

close all
figure 
%barh(k_bar,1,'stack','FaceColor','none','Linewidth',1);
b = barh(k_bar(1:4,:),1,'stack','FaceColor','none','Linewidth',1.25,'EdgeColor','k');
yes_color = 0;
if yes_color 
	c_num = 15;
	c_map = colormap(copper(c_num));
	for i = 2 : length(space)
		b(i).FaceColor = c_map(mod(i,c_num)+1,:);
	end
end
if N > 0
	xlabel(['$\lambda_\alpha b_{am} + 2 N \lambda_\alpha \pi$, $0 \leq N \leq $' int2str(N)],'Interpreter', 'latex')
	xlim([0,7])
else 
	xlabel(['$\lambda_\alpha b_{am} + 2 N \lambda_\alpha \pi$, ', int2str(N) ' $\leq N \leq 0$'],'Interpreter', 'latex')
	xlim([-7,0])
end
plot_para('Maximize',true,'Filename','1')

track = 3;
fprintf ('Real phase: %f\n', 4*pi/lamda_(track)*(d_a*y_0/R_0) )
fprintf ('F band: %i, Amb number: %i, Real phase: %f\n', track, 2, mod(4*pi/lamda_(track)*(d_a*y_0/R_0), 2*pi) + 2 * 2*pi )
fprintf ('For each interval, y_0 will differ %f\n', 10)
%{
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