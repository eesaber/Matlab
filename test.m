%% Aim:  Verify: High-order Ambiguity function
	close all
	clear
	%Time space
	tot_point = 1024*6;
	t_1 = linspace(0, 1, tot_point);
	% a_m parameters
	c_0 = 10 ;
	c_1 = 100 ;
	c_2 = 500 ;
	c_3 = -100;
	c_4 = 0 ;
	m = 3 ;
	% Siganl 
	
	s = exp(j * 2 * pi * (c_0 + c_1 *t_1 + c_2 * t_1.^2 / 2 + c_3 *t_1.^3 / 6+ c_4 * t_1.^4 / 24));
	AF = GAF(s,3,3,2);
	[~,qq] = max(abs(AF));
	f = linspace(-tot_point/2, tot_point/2, length(AF));
	xi = 1 / 3 / 2 ;
	t_c_3 = f(qq) / 4 / xi^2;
	fprintf('c_3 = %f, Estimated c_3 = %f \n', c_3, t_c_3)
	figure
		plot(f,abs(AF),'k','Linewidth',3.5)
		set(gcf,'color','w');
		set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
		xlabel('Hz', 'Interpreter', 'latex')
		xlim([-100 100])
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1); 
		%export_fig AF_3.jpg
	s = s .* exp(-j * 2 * pi * c_3 * t_1.^3 / 6) ;
	AF = GAF(s,2,2,2);
	[~,qq] = max(abs(AF));
	xi = 1 / 2 / 2;
	t_c_2 = f(qq) / 2/ xi;
	fprintf('c_2 = %f, Estimated c_2 = %f \n', c_2, t_c_2)

	figure
		plot(f,abs(AF),'k','Linewidth',3.5)
		set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
		set(gcf,'color','w');
		xlabel('Hz', 'Interpreter', 'latex')
		xlim([100 400])
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1); 
		%export_fig AF_2.jpg
%% Aim: Total points will not makes resolution better
% We find GAF has a poor resolution, and want to chekc if add more points
% can improve the resolution. 
	clear AF
	close all
	tot_point = 1024*2;
	t_1 = linspace(0, 1, tot_point) ;
	c_3 = linspace(1,25,50);
	c_2 = linspace(-20,20,100);
	t_c_3 = zeros(5,length(c_3));
	t_c_2 = zeros(5,length(c_3),length(c_2));
	for k = 0 : 4
		f = linspace(-tot_point/2, tot_point/2, tot_point*2^k); 
		for n = 1 : length(c_3)
			for g = 1 : length(c_2)
				s = exp(j * 2 * pi * (c_0 + c_1 *t_1 + c_2(g) * t_1.^2 / 2 + c_3(n) *t_1.^3 / 6+ c_4 * t_1.^4 / 24));
				AF = GAF(s,3,3,k);
				[~,qq] = max(abs(AF));
				xi = 1 / 3 / 2 ;
				t_c_3(k+1,n) = f(qq) / 4 / xi^2;
				s = s .* exp(-j * 2 * pi * t_c_3(k+1,n) * t_1.^3 / 6) ;
				AF = GAF(s,2,2,k);
				[~,qq] = max(abs(AF));
				xi = 1 / 2 / 2;
				t_c_2(k+1,n,g) = f(qq) / 2/ xi;
			end
		end
	end

	% Fig.1 Quantization error of c_3
	%{
		left_color = [1 0.6 0.2];
		right_color = [0 0 0];
		set(figure,'defaultAxesColorOrder',[left_color; right_color]);
		yyaxis left
		plot(c_3, t_c_3(1,:) - c_3,'--','Linewidth',3.5)
		yyaxis right
		plot(c_3, t_c_3(4,:) - c_3, c_3, t_c_3(5,:) - c_3,'-.','Linewidth',3.5)
	%}
	% Fig.2 Show three lines of "c_3" in different total point
		plot(c_3, t_c_3(1,:),'--k', c_3, t_c_3(4,:) ,'k', c_3, t_c_3(5,:) ,'k-.','Linewidth',3.5);
		set(gca,'FontSize',40,'Fontname','CMU Serif Roman','Linewidth',2)
		set(gcf,'color','w');
		xlabel('$c_3$', 'Interpreter', 'latex')
		ylabel('$\tilde{c}_3 - c_3$', 'Interpreter', 'latex')
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1); 
		%export_fig c3ResoDiffN.jpg
	%% Subaim: The effect of cccuracy of c_3 when estimate c_2
	figure
		%{surf(c_2,c_3,squeeze(t_c_2(5,:,:)) - repmat(c_2,length(c_3),1) )
		xlabel('$c_2$', 'Interpreter', 'latex')
		ylabel('$c_3$', 'Interpreter', 'latex')
		zlabel('$\tilde{c}_2 - c_2$', 'Interpreter', 'latex')
		%}
		imagesc(c_2,c_3,squeeze(t_c_2(5,:,:)) - repmat(c_2,length(c_3),1) )
		xlabel('$c_2$', 'Interpreter', 'latex')
		ylabel('$c_2 - \tilde{c}_2$', 'Interpreter', 'latex')
		colorbar
		colormap('Jet')

		set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
		set(gcf,'color','w');
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1); 
		export_fig c2Contour.jpg
%% Aim: Plot the contour of c_3. 
%I want to show that c_3 changes little when v_y changes a lot 
	clear
	% parameter 
	T_a = 1;
	eta = linspace(-T_a,T_a,T_a * 2000);
	x_0 = 1000/1.41; y_0 = 0; h = 1000/1.41; 
	v_r = linspace(-25,25,52); v_y = linspace(-20,20,81);
	v_p = 100;
	f_0 = 5e9; c = 3e8 ; lambda = c/f_0 ;
	% vector space
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
	figure
		imagesc(co_3)
		xlabel('$v_y$', 'Interpreter', 'latex')
		ylabel('$v_r$', 'Interpreter', 'latex')
		%xlim([-10 10])
		%ylim([-10 10])
		hold on
		[C,h] = contour(co_3,10,'k--','LineWidth',3,'ShowText','on');
		clabel(C,h,'FontSize',40,'Color','black','LabelSpacing',1000)
		set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
		set(gcf,'color','w');
		set(gca,'linewidth',2.5)
		colormap('Jet')
		colorbar
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1); 
		export_fig c3Contour.jpg
	figure
		imagesc(v_y, v_r, co_3)
		xlabel('$v_y$', 'Interpreter', 'latex')
		ylabel('$v_r$', 'Interpreter', 'latex')
		set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
		set(gcf,'color','w');
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1); 
		export_fig c3Contour1.jpg
		set(gca,'xtick',[-20:5:20],'ytick',[-25:5:25]);
		set(gca,'linewidth',2.5)
		colorbar
		export_fig c3Contour1.jpg
%% Range equation expansion 
	syms t x_0 y_0 v_x v_y v_p a_x a_y h 
	f = sqrt( (x_0 + v_x * t + a_x* t^2 / 2)^2 + h^2 ...
		+ (y_0 + (v_y - v_p) * t + a_y * t^2 / 2 )^2  );
	taylor(f, t, 'order', 2)
	taylor(f, t, 'order', 3)
	taylor(f, t, 'order', 4) 

