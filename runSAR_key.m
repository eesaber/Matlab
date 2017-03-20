v_space = linspace(-20,20,41);
if ispc 
	cd D:\Code\Simu
else
	cd ~/Code/Matlab 
end
%% v_y change , v_x = 3
v_yt = zeros(1,41); v_rt = zeros(1,41); a_rt = zeros(1,41);
co_3 = zeros(1,41);
c_3 = zeros(1,41);
for m = 1 : length(v_space)
   [co_3(m), c_3(m), v_rt(m), v_yt(m), a_rt(m)] = SAR_key(5, v_space(m),10,0);
end

%%
close all
figure
	plot(v_space, v_yt  ,'k','Linewidth',4)
	xlabel('$v_y$ (m/s)','Interpreter','latex')
	ylabel('$\tilde{v}_y - v_y$ (m/s)','Interpreter','latex')
	set(gcf,'color','w');
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman','Linewidth',3)
	pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
	export_fig vyError.jpg 
figure
	plot(v_space, a_rt,'k','Linewidth',4)
	xlabel('$v_y$ (m/s)','Interpreter','latex')
	ylabel('$\tilde{a}_r$ (m/s)','Interpreter','latex')
	set(gcf,'color','w');
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman','Linewidth',3)
	pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
	export_fig ayError.jpg
%{
%%
my_plot(v_rt - v_r, v_space, '$\tilde{v}_r - v_r$(m/s)','ErrKey_vx')

%% v_x = 0 , v_y change
v_yt = zeros(1,41);
for m = 1 : 1
   [v_r(m), v_rt(m), v_yt ] = SAR_key(v_space(m), v_space(m), 0, 0);
end
%%
my_plot(v_yt - v_space, v_space, '$\tilde{v}_y - v_y$(m/s)','ErrKey_vy')
%}