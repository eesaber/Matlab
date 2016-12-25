function my_plotyy (s1, s2, x_s, l_y_la, r_y_la,ti)
	figure
    [ax,hLine1,hLine2] = plotyy(x_s, s1, x_s, s2);
    qq = [hLine1,hLine2] ;
    hLine2.LineStyle = '-.';
	xlabel('$v$ (m/s)','Interpreter','latex')
    ylabel (ax(1), l_y_la,'Interpreter','latex')
    ylabel (ax(2), r_y_la,'Interpreter','latex')
    set(qq,'color','k')
    set(qq,'linewidth',2.5)
	%set(ax,'FontSize',42,'Fontname','CMU Serif Roman') % 4:3 moniter at EE527
    set(ax,'FontSize',36,'Fontname','CMU Serif Roman') % 4:3 moniter at laptop    
	%pbaspect([4 3 1])
    saveas(gcf, [ti '.jpg'])
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
end

