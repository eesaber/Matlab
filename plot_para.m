function plot_para(varargin)
	
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'char'},{'nonempty'}); 
	validationFcn_2_ = @(x) validateattributes(x,{'logical'},{'scalar'});
	validationFcn_3_ = @(x) validateattributes(x,{'numeric'},{'positive'});
    validationFcn_4_ = @(x) validateattributes(x,{'numeric'},{'positive'});
	addParameter(parse_,'Filename','None',validationFcn_1_);
	addParameter(parse_,'Maximize',0,validationFcn_2_);
	addParameter(parse_,'Ratio',[0, 0, 0],validationFcn_3_);
    addParameter(parse_,'Fontsize',0,validationFcn_4_);
	parse(parse_,varargin{:})

	%%
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
    set(gcf,'color','w');
	get_size = get(0,'screensize');
    if parse_.Results.Fontsize 
        font_size = parse_.Results.Fontsize ;
    else
        if get_size(3:4) == [1920 1080]
            font_size = 40;
        else
            font_size = 28;
        end
    end
    set(gca,'FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
    %set(gca,'FontSize',font_size,'Fontname','CMU Sans Serif','Linewidth',2)
    
	if parse_.Results.Ratio(1)
		pbaspect(parse_.Results.Ratio)
	end

	if parse_.Results.Maximize
        pause(0.01);
		frame_h = get(handle(gcf),'JavaFrame');
		set(frame_h,'Maximized',1);
    end
    pause(0.5);
	if strcmp(parse_.Results.Filename,'None') == 0
		a = [parse_.Results.Filename, '.jpg'];
		export_fig(a)
		clear a
	end
end