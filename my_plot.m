function my_plot(s, x_s, y_la)
print = true ;
if(print) 
%%It will print the magnitude of the input signal 
	figure
	plot(x_s, s, 'Linewidth',4)
	xlabel('$v_x$','Interpreter','latex')
    ylabel(y_la,'Interpreter','latex')
	set(gca,'FontSize',32,'Fontname','CMU Serif Roman')
	pbaspect([4 3 1])
end

end