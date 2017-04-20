function purinto(s, label_x, label_y , x, y)
% Usage: purinto(s, label_x, label_y), where "s" is the 2D matrix,
% "label_x" is the label of X-axis and "label_y" is label of Y-axis. If no
% label input, it will be default to print "range" for X-axis and "Azimuth"
% for Y-axis.

   
	get_size = get(0,'screensize');
	if get_size(3:4) == [1920 1080]
		font_size = 40;
	else
		font_size = 28;
	end
%%It will print the magnitude of the input signal 
	figure
	imagesc(abs(s))
    if nargin == 5
        imagesc(x, y, abs(s))
        %axis off 
    end
    % Maximize the figure
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
    set(gcf,'color','w');
	set(gca,'Ydir','normal')   
	set(gca,'FontSize',font_size,'Fontname','CMU Serif Roman')
	colormap('Jet')
	colorbar
	pbaspect([7 5 1])
    
    
end