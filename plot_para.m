function plot_para(Spec, Maximize,Filename)
	
	if nargin == 1 
		Maximize = 0 ;
	end
	if nargin ~= 3
		Filename = 'Default';
	end
	%%
	if Spec == 1 
		set(gca,'FontSize',30,'Fontname','CMU Serif Roman','Linewidth',2)
		set(gcf,'color','w');
		pbaspect([7 5 1])
	end
	
	if Maximize
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1);
	end
	if strcmp(Filename,'Default') == 0
		a = [Filename, '.jpg'];
		export_fig(a)
		clear a
	end
end