% Verify the fact that v_y is very hard to estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Result								%
% <a_y, v_y> == <0, change>				%
% The error is small					%
% <a_y, v_y> == <3, change>				%
% The error is very large				%
% Conclusion: v_y needs new approach	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_space = linspace(-20,20,41);
if ispc 
	cd D:\Code\Simu
else
	cd ~/Code/Matlab 
end
%% v_y change , v_x = 3

v_yt = zeros(1,41); v_rt = zeros(1,41); a_rt = zeros(1,	41);
co_3 = zeros(1,41);
c_3 = zeros(1,41);
for m = 1 : length(v_space)	
   [co_3(m), c_3(m), v_rt(m), v_yt(m), a_rt(m)] = SAR_key(10, v_space(m),3,0);
end
close all
figure
	plot(v_space, v_yt - v_space  ,'k','Linewidth',4)
	xlabel('$v_y$ (m/s)','Interpreter','latex')
	ylabel('$\tilde{v}_y - v_y$ (m/s)','Interpreter','latex')
	set(gcf,'color','w');
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman','Linewidth',3)
	pbaspect([7 5 1])
	pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
	export_fig vyError1.jpg 
figure
	plot(v_space, a_rt-3,'k','Linewidth',4)
	xlabel('$v_y$ (m/s)','Interpreter','latex')
	ylabel('$\tilde{a}_r - a_r$ (m/s$^2$)','Interpreter','latex')
	set(gcf,'color','w');
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman','Linewidth',3)
	pbaspect([7 5 1])
	pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
	export_fig arError1.jpg
