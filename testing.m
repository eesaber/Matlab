%% Plot Each Region onto H/alpha diagram
% bare soil
fan = [343,461;214,368]; % [row, row; col, col]
size_f = (fan(1,2)-fan(1,1)+1)*(fan(2,2)-fan(2,1)+1);
temp = cat(1,cat(2, reshape(T_11(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])), ......
             cat(2,reshape(conj(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_22(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])),.......
             cat(2,reshape(conj(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(conj(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_33(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])));
ind_ = randperm(size_f);
figure
H_Alpha(temp(:,:,ind_(1:500)),'Pointcolor','b')
% building 
fan = [176,273; 1,212]; % [row, row; col, col]
size_f = (fan(1,2)-fan(1,1)+1)*(fan(2,2)-fan(2,1)+1);
temp = cat(1,cat(2, reshape(T_11(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])), ......
             cat(2,reshape(conj(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_22(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])),.......
             cat(2,reshape(conj(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(conj(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_33(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])));
ind_ = randperm(size_f);
H_Alpha(temp(:,:,ind_(1:500)),'Pointcolor','r')
% canopy
fan = [214,285; 661,801]; % [row, row; col, col]
size_f = (fan(1,2)-fan(1,1)+1)*(fan(2,2)-fan(2,1)+1);
temp = cat(1,cat(2, reshape(T_11(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])), ......
             cat(2,reshape(conj(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_22(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])),.......
             cat(2,reshape(conj(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(conj(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_33(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])));
ind_ = randperm(size_f);
H_Alpha(temp(:,:,ind_(1:500)),'Pointcolor','k')

%%
syms x y z a b c e f g r 
r = 1;
eqn1 = r*e*x + y + z == a;
eqn2 = r*f*x + y - z == b;
eqn3 = r*x + 1j*g*y + b*z == c;

sol = solve([eqn1, eqn2, eqn3], [x, y, z]);
xSol = sol.x
ySol = sol.y
zSol = sol.z
%%
clear
syms U_x U_y U_z D_x D_y D_z 
syms k_x k_y k_z kp_z k kp 
syms h h_x h_y xi
syms T R c_t c_i s_t x
syms KEEP %  1: Keep partial derivative, 0: Drop it
KEEP = 0;

eqn1 = k_x*U_x + k_y*U_y + k_z*U_z == 0;
eqn2 = -k_x*D_x - k_y*D_y + kp_z*D_z == 0;
eqn3 = (1-1j*k_z*h)*U_x + KEEP*h_x*(1-1j*k_z*h)*U_z - (1+1j*kp_z*h)*D_x - KEEP*h_x*(1+1j*kp_z*h)*D_z == 0;
eqn4 = (1-1j*k_z*h)*U_y + KEEP*h_y*(1-1j*k_z*h)*U_z - (1+1j*kp_z*h)*D_y - KEEP*h_y*(1+1j*kp_z*h)*D_z == 0;

eqn5 = -1j*k_y*(1-1j*k_z*h)*KEEP*h_x*U_x + ( -1j*k_z - k_z^2*h + 1j*k_x*(1-1j*k_z*h)* KEEP*h_x)*U_y + 1j*k_y*(1-1j*k_z*h)*U_z ......
		+ 1j*k_y*(1+1j*kp_z*h)*KEEP*h_x*D_x + (-1j*kp_z + kp_z^2*h - 1j*k_x*(1+1j*kp_z*h)*KEEP*h_x)*D_y - 1j*k_y*(1+1j*kp_z*h)*D_z ......
        == -(T*kp^2*c_t^2 - (T-1)*k^2*c_i^2)*xi;
%		== -(T*kp^2*c_t^2 - R*k^2*c_i^2)*xi;
    
eqn6 = (-1j*k_z - k_z^2*h - 1j*k_y*(1-1j*k_z*h)*KEEP*h_y)*U_x + 1j*k_x*(1-1j*k_z*h)*KEEP*h_y*U_y + 1j*k_x*(1-1j*k_z*h)*U_z ......
        + (-1j*kp_z + kp_z^2*h + 1j*k_y*(1+1j*kp_z*h)*KEEP*h_y)*D_x - 1j*k_x*(1+1j*kp_z*h)*KEEP*h_y*D_y - 1j*k_x*(1+1j*kp_z*h)*D_z == 0;
    
sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [U_x, U_y, U_z, D_x, D_y, D_z]);

sol.U_x
sol.U_y



%%
%{
s_p = 5;
[x, y] = meshgrid(linspace(-1, 1, s_p), linspace(-1, 1, s_p))
[a,b,c]=sphere(4)
x(sqrt(x.^2 + y.^2)>1) = nan;
y(sqrt(x.^2 + y.^2)>1) = nan;
z = real(sqrt(1-x.^2-y.^2));
surf(z,'FaceColor','w')

syms U_x U_y U_z D_x D_y D_z 
syms k_x k_y k_z kp_z k kp 

syms h h_x h_y xi
syms T R c_t c_i s_t


eqn1 = k_x*U_x + k_y*U_y + k_z*U_z == 0;
eqn2 = -k_x*D_x - k_y*D_y + kp_z*D_z == 0;
eqn3 = (1-1j*k_z*h)*U_x + h_x*(1-1j*k_z*h)*U_z - (1+1j*kp_z*h)*D_x - h_x*(1+1j*kp_z*h)*D_z == 0;
eqn4 = (1-1j*k_z*h)*U_y + h_x*(1-1j*k_z*h)*U_z - (1+1j*kp_z*h)*D_y - h_y*(1+1j*kp_z*h)*D_z == 0;
eqn5 = -1j*k_y*(1-1j*k_z*h)*h_x*U_x + [(-1j*k_z - k_z^2*h) + 1j*k_x*(1-1j*k_z*h)*h_x]*U_y + 1j*k_y*(1-1j*k_z*h)*U_z ......
	+ 1j*k_y*(1+1j*kp_z*h)*h_x*D_x +[-1j*kp_z + kp_z^2*h - 1j*k_x*(1+1j*kp_z*h)*h_x]*D_y - 1j*k_y*(1+1j*kp_z*h)*D_z == 0;
eqn6 = [-1j*k_z - k_z^2*h -1j*k_y*(1-1j*k_z*h)*h_y]*U_x + 1j*k_x*(1-1j*k_z*h)*h_x*U_y + 1j*k_x*(1-1j*k_z*h)*U_z ......
	+ [-1j*kp_z + kp_z^2*h + 1j*k_y*(1+1j*kp_z*h)*h_y]*D_x - 1j*k_x*(1+1j*kp_z*h)*h_y*D_y - 1j*k_y*(1+1j*kp_z*h)*D_z == 0;

sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [U_x, U_y, U_z, D_x, D_y, D_z]);
sol.U_x
sol.U_y

%}
%{
phi = linspace(-pi,pi,1000);
E = cos(0.6*pi*cos(phi))-cos(0.6*pi*sin(phi));
E = abs(E)/max(abs(E));
figure(1)
plot(phi,E)

L = linspace(0.1,30);
a = 2*L.*tan(acos(L./(L+0.25))); 
figure(2)
	plot(L,a,'k','Linewidth',3)
	xlabel('$L/\lambda$','Interpreter', 'latex')
	ylabel('$a/\lambda$','Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','ant1')
figure(3)
	plot(L,10*log10(7.5*a.^2),'k','Linewidth',3)
	xlabel('$L/\lambda$','Interpreter', 'latex')
	ylabel('$D$ (dB)','Interpreter', 'latex')
	plot_para('Maximize',true,'Filename','ant2')
%}

%{
syms a b c d co si;
R_2 = [co, si; -si, co];
R_3 = [1+co, sqrt(2)*si, 1-co; -sqrt(2)*si, 2*co, sqrt(2)*si; 1-co, -sqrt(2)*si, 1+co];
S = [a c; c b];
S = R_2.'*S*R_2
rho = [co*(a*co + b*si) + si*(b*co + a*si); sqrt(2)*(co*(b*co + a*si) - si*(a*co + b*si)); co*(a*co - b*si) - si*(b*co - a*si)];
rho*conj(rho).';
%C = 2*[a^2, 0, a*conj(d); 0, 2*b^2, 0; conj(a)*d, 0, d^2];
syms a11 a12 a13 a21 a22 a23 a31 a32 a33
C = [a11, a12, a13; a21, a22, a23; a31, a32, a33];
C = 1/4*R_3*C*R_3.';

%}

%{
%% DATA generate
f1 = 9.6e9;
f2 = 1.9e9;
f3 = 1.27e9;
%f_ = -sort(-[f1, f2, f3]) ;
f_ = -sort(-[f1, f2, f1-f2, f1+f2]);
lamda_ = c./f_;	
N=20;
T = 100;
N_down = min(floor(4*pi./lamda_*d_a/2/R_0*(-T*(v_p - v_y) /2 + y_0) / 2 /pi));
N_up = max(floor(4*pi./lamda_*d_a/2/R_0*(T*(v_p - v_y) /2 + y_0) / 2 /pi));
space = 2*pi*(N_down: N_up);
clear k k_bar temp 
temp = mod(4*pi./lamda_*(d_a*y_0/2/R_0), 2*pi) .* lamda_;
k = repmat(temp.',1,length(space)) + lamda_.' * space ;
k_bar = [k(:,1),  diff(k,1,2)];

close all
%% Plot the 
x_l = 6;
figure 
subplot('position',[0.1 0.5 0.8 0.2]);
%subplot('position',[0.1 0.3375 0.8 0.1125]);
barh(k_bar,1,'stack','FaceColor','none','Linewidth',1.25,'EdgeColor','k' ...
	,'BaseValue', -100, 'ShowBaseLine', 'on');
set(gca,'TickLength',[0 0])
ylim([0.5 numel(k)/length(k)+0.5])

%yticklabels({'\lambda_1 ', '\lambda_2','\lambda_3'})
%yticks([0.75 2 3 4.2])
yticklabels({'\lambda_4 ', '\lambda_3', '\lambda_2', '\lambda_1'})
hold all
x_l = 2*x_l*(N>0) - x_l;
xlim([-x_l,x_l])
%xlabel(['$( \beta_k \hat{y}_k + 2 N \pi)\lambda_k$, ', int2str(N_down), '$ \leq N \leq $', int2str(N_up)],'Interpreter', 'latex')
xlabel('$(\beta_\alpha \tilde{y}_\alpha + 2 N \pi)\lambda_\alpha$','Interpreter', 'latex')
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

annotation('arrow', (0.5+ 0.8*(4*pi*d_a*y_0/2/R_0 / 2 / x_l ))*ones(1,2), [.69 .7],'Color','k')
plot(4*pi*d_a*y_0/2/R_0*[1 1],[.5 4.5],'Color', 'k', 'Linewidth',2.5)
plot_para('Maximize',true,'Filename','AmbiGene1')



track = 3;
fprintf ('Real phase: %f\n', 4*pi/lamda_(track)*(d_a*y_0/2/R_0) )
fprintf ('F band: %i, Amb number: %i, Real phase: %f\n', track, 2, mod(4*pi/lamda_(track)*(d_a*y_0/2/R_0), 2*pi) + 2 * 2*pi )
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