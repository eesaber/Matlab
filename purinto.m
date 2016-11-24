function purinto(s)
print = true ;
if(print) 
%%It will print the magnitude of the input signal 
    figure
    imagesc(abs(s))
    set(gca,'Ydir','normal')   
    xlabel('Range $(\Delta \tau)$','Interpreter','latex')
    ylabel('Azimuth $(\Delta \eta)$','Interpreter','latex')
    set(gca,'FontSize',26,'Fontname','CMU Serif')
    colorbar
    pbaspect([4 3 1])
end