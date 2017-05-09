v = -20:2:20;
v_Est_Vy_ = zeros(length(v),length(v));
v_Est_Err_Vy_ = v_Est_Vy_;

for vx = 1 : length(v)
	for vy = 1 : length(v)
		v_Est_Vy_(vx,vy) = DualRx(v(vx),v(vy),0,0);
		v_Est_Err_Vy_(vx,vy) =  v_Est_Vy_(vx,vy) - v(vy);
	end
	fprintf('.')
	if mod(vx,10) == 0
		fprintf('\n')
	end
end
%%
figure
imagesc(v,v,v_Est_Err_Vy_)
	xlabel('$v_y$', 'Interpreter', 'latex')
	ylabel('$v_x$', 'Interpreter', 'latex')
	set(gca,'Ydir','normal'),colorbar, colormap('Jet')
	caxis([-1 1])
%hold on 
%[C,h] = contour(v,v,v_Est_Err_Vy_,'-k','linewidth',1.5,'ShowText','on');
%clabel(C,h,'FontSize',40,'Color','black','LabelSpacing',1000)
	plot_para('Maximize',true,'Filename','errmap')

%%
find (abs(v_Est_Err_Vy_) > 10);
temp = reshape(v_Est_Err_Vy_,1,41^2);
temp(find(abs(v_Est_Err_Vy_) > 10)) = 1;
temp = reshape(temp,41,41);
figure
imagesc(v,v,temp)
	xlabel('$v_y$', 'Interpreter', 'latex')
	ylabel('$v_x$', 'Interpreter', 'latex')
	set(gca,'Ydir','normal')
	colormap('Jet')
	caxis([-10, 2])
	colorbar
	plot_para(1,1,'VaziErr')