function my_plot(s, x_s, y_la, ti)
print = true ;
if(print) 
%%It will print the magnitude of the input signal 
	figure
	plot(x_s, s, 'Linewidth',4)
    axis([0 (length(s) - 1) -inf inf ])
	xlabel('$v_x$','Interpreter','latex')
    ylabel(y_la,'Interpreter','latex')
	set(gca,'FontSize',26,'Fontname','CMU Serif Roman')
	pbaspect([4 3 1])
    saveas(gcf, [ti '.jpg'])
end

end