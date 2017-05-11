v = -20:2:20;
v_Est_Vy_ = zeros(length(v),length(v));
v_Est_Err_Vy_ = v_Est_Vy_;

%{
v_Est_Vy_wvd = zeros(length(v),length(v));
v_Est_Vy_crr = v_Est_Vy_wvd;
v_Est_Err_Vy_wvd = zeros(length(v),length(v));
v_Est_Err_Vy_crr = v_Est_Vy_wvd;
%}

for vx = 1 : length(v)
	for vy = 1 : length(v)
		%v_Est_Vy_(vx,vy) = DualRx(v(vx),v(vy),0,0);
		%v_Est_Err_Vy_(vx,vy) =  v_Est_Vy_(vx,vy) - v(vy);
		[v_Est_Vy_wvd(vx,vy), v_Est_Vy_crr(vx,vy)]= DualRx(v(vx),v(vy),0,0);
		v_Est_Err_Vy_crr(vx,vy) = v_Est_Vy_crr(vx,vy) - v(vy);
		v_Est_Err_Vy_wvd(vx,vy) = v_Est_Vy_wvd(vx,vy) - v(vy);
	end
	fprintf('.')
	if mod(vx,10) == 0
		fprintf('\n')
	end
end
%{
figure
imagesc(v,v,v_Est_Err_Vy_wvd)
	xlabel('$v_y$', 'Interpreter', 'latex')
	ylabel('$v_x$', 'Interpreter', 'latex')
	set(gca,'Ydir','normal'),colorbar, colormap('Jet')
	%caxis([-1 1])
	plot_para('Maximize',true,'Filename','errmap')

imagesc(v,v,v_Est_Err_Vy_crr)
	xlabel('$v_y$', 'Interpreter', 'latex')
	ylabel('$v_x$', 'Interpreter', 'latex')
	set(gca,'Ydir','normal'),colorbar, colormap('Jet')
	%caxis([-1 1])
	plot_para('Maximize',true,'Filename','errmap')
%}
	

figure
imagesc(v,v,v_Est_Err_Vy_)
	xlabel('$v_y$', 'Interpreter', 'latex')
	ylabel('$v_x$', 'Interpreter', 'latex')
	set(gca,'Ydir','normal'),colorbar, colormap('Jet')
	caxis([-1 1])
	plot_para('Maximize',true,'Filename','errmap')