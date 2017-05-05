function plot_para(varargin)
	
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'char'},{'nonempty'}); 
	validationFcn_2_ = @(x) validateattributes(x,{'logical'},{'scalar'}); 
	addParameter(parse_,'Filename','None',validationFcn_1_);
	addParameter(parse_,'Maximize',0,validationFcn_2_);
	parse(parse_,varargin{:})

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

	if parse_.Results.Maximize
		pause(0.00001);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1);
	end
	if strcmp(parse_.Results.Filename,'None') == 0
		a = [parse_.Results.Filename, '.jpg'];
		export_fig(a)
		clear a
	end
end