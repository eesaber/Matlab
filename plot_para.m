function plot_para(Maximize,Filename)
	parser = inputParser;

	if nargin == 0
		Maximize = 0 ;
	end
	if nargin ~= 2
		Filename = 'Default';
	end
	%%
	get_size = get(0,'screensize');
	if get_size(3:4) == [1920 1080]
		font_size = 40;
	else
		font_size = 28;
	end
	set(gca,'FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
	set(gcf,'color','w');
	pbaspect([7 5 1])
	
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