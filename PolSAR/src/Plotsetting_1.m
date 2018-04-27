function Plotsetting_1(Clim)
    n = 4;
    set(gca,'Ydir','normal','Clim',Clim)
    switch n
        case 1
        set(gca,'View',[90 90], 'XTick',1:900:4501)
        xt=arrayfun(@num2str,get(gca,'xtick')+500-1,'un',0);
        yt=arrayfun(@num2str,get(gca,'ytick')+30500,'un',0);
        set(gca,'xticklabel',xt,'yticklabel',yt)
        case 2 
        set(gca, 'View',[90 90], 'XTick',1:500:2501)
        xt=arrayfun(@num2str,get(gca,'xtick')+2000-1,'un',0);
        yt=arrayfun(@num2str,get(gca,'ytick')+22000,'un',0);
        set(gca,'xticklabel',xt,'yticklabel',yt)
        case 3
        set(gca, 'View',[90 90], 'XTick',1:300:1401)
        xt=arrayfun(@num2str,get(gca,'xtick')+600-1,'un',0);
        yt=arrayfun(@num2str,get(gca,'ytick')+26000,'un',0);
        set(gca,'xticklabel',xt,'yticklabel',yt)
        otherwise
            1;
    end
    colormap jet; colorbar    
end