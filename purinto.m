function purinto(s, label_x, label_y , x, y)
% Usage: purinto(s, label_x, label_y), where "s" is the 2D matrix,
% "label_x" is the label of X-axis and "label_y" is label of Y-axis. If no
% label input, it will be default to print "range" for X-axis and "Azimuth"
% for Y-axis.

    if nargin == 1
        label_x = 'Range $(\Delta \tau)$'; 
        label_y = 'Azimuth $(\Delta t)$'; 
    end 
%%It will print the magnitude of the input signal 
	figure
    if nargin ~= 3
        imagesc(abs(s))
    end
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
	xlabel(label_x, 'Interpreter', 'latex')
	ylabel(label_y, 'Interpreter', 'latex')
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
	colormap('Jet')
	colorbar
	%pbaspect([4 3 1])
    
    
end