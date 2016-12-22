function my_plot(s1, s2, x_s, y_la, ti)


%%It will print the magnitude of the input signal 
	figure
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
	plot(x_s, s1,'k-.' ,'Linewidth',4)
    if s2 ~= 0
        hold
        plot(x_s, s2,'k' ,'Linewidth',4)
    end
    ymin = floor(min( [min(s1) min(s2)]));
    ymax = ceil(max([max(s1) max(s2)] ) );
    axis([0 (length(s1) - 1) ymin ymax ])
	xlabel('$v$ (m/s)','Interpreter','latex')
    ylabel(y_la,'Interpreter','latex')
    set(gca,'linewidth',2.5)
	set(gca,'FontSize',42,'Fontname','CMU Serif Roman')
	pbaspect([4 3 1])
    saveas(gcf, [ti '.jpg'])
    
end