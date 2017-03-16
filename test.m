%% High-order Ambiguity function 
close all
clear
%Time space
t_mov = 1/6;
tot_point = 1024*6;
t_1 = linspace(0, 1, tot_point);
% a_m parameters
a_0 = 10 ;
a_1 = 100 ;
a_2 = 500 ;
a_3 = linspace(1,100,100) ;
a_4 = 000 ;
% Siganl 
f = linspace(-tot_point/2, tot_point/2, 2^(nextpow2(length(t_1)-1)) );
t_a_3 = zeros(1,length(a_3));
for i = 1 : length(a_3)
	s = exp(j * 2 * pi * (a_0 + a_1 *t_1 + a_2 * t_1.^2 / 2 + a_3(i) *t_1.^3 / 6+ a_4 * t_1.^4 / 24)) ;
	AF = GAF(s,3,3);
	[~,qq] = max(abs(AF));
	xi = 1 / 3 / 2 ;
	t_a_3(i) = f(qq) / 4 / xi^2;
	fprintf('a_3 = %f, Estimated a_3 = %f \n', a_3(i), t_a_3(i))
end

%% the change of c_3 
clear
% parameter 
eta = linspace(-10,10,10000);
x_0 = 1000/1.41; y_0 = 0; h = 1000/1.41; 
v_r = linspace(-25,25,52); v_y = linspace(-10,10,81);
v_p = 100;
f_0 = 5e9; c = 3e8 ; lambda = c/f_0 ;
co_3 = zeros(length(v_r),length(v_y));
R_0 = zeros(length(v_r),length(v_y));
% generate c_3
for h = 1 : length(v_r)
	for k = 1 : length(v_y)
		R = sqrt(h^2 + (x_0 + v_r(h)* eta ).^2 + (y_0 + v_y(k)* eta - v_p* eta).^2 );
		[R_0(h,k), ~] = min(R);
		co_3(h,k) = -1 / lambda * 6* v_r(h) * (v_y(k) - v_p)^2 / R_0(h,k)^2 ;
	end
end
%%
close all
figure
	imagesc(co_3)
	xlabel('$v_r$', 'Interpreter', 'latex')
	ylabel('$v_y$', 'Interpreter', 'latex')
	%xlim([-10 10])
	%ylim([-10 10])
	hold on
	[C,h] = contour(co_3,10,'k--','LineWidth',3,'ShowText','on');
	clabel(C,h,'FontSize',40,'Color','black','LabelSpacing',1000)
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
	set(gcf,'color','w');
	colormap('Jet')
	colorbar
	pause(0.00001);
	frame_h = get(handle(gcf),'JavaFrame');
	set(frame_h,'Maximized',1); 
	export_fig c3Contour.jpg
	
	
figure
	imagesc(v_r, v_y, co_3)
	xlabel('$v_r$', 'Interpreter', 'latex')
	ylabel('$v_y$', 'Interpreter', 'latex')
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
	set(gcf,'color','w');
	pause(0.00001);
	frame_h = get(handle(gcf),'JavaFrame');
	set(frame_h,'Maximized',1); 
	export_fig c3Contour1.jpg
%set(gca,'xtick',[-25:25],'ytick',[-10:10]);
%{
figure
    plot(f,abs(AF_2),'k','Linewidth',3.5)
    set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
    set(gcf,'color','w');
    xlabel('Hz', 'Interpreter', 'latex')
    xlim([-1000 0])
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
    %export_fig AF_2.jpg
figure
    plot(f,abs(AF_3),'k','Linewidth',3.5)
    set(gcf,'color','w');
    set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
    xlabel('Hz', 'Interpreter', 'latex')
    xlim([-100 100])
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
    %export_fig AF_3.jpg
	



%% Expand the High order Generalized Ambiguity function 
clc
clear
syms x xi_1 temp
for m = 1 : 5
    
    expand( (x + xi_1)^m )/ factorial(m) - expand( (x - xi_1)^m )/ factorial(m)
end

%% 
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