function Plotsetting_1(Clim, varargin)
    parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'numeric'},{});
	addParameter(parse_,'Colorbar_unit',[],validationFcn_1_);
	parse(parse_,varargin{:})
    
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
    set(gca,'XTick',0:5000:30000, 'Xticklabel',cellstr(int2str((0:5000*7:30000*7)'/1000))')
    set(gca,'YTick',0:600:3300, 'Yticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
    
    
    colormap jet; colorbar
    if parse_.Results.Colorbar_unit
        title(colorbar,'(dB)','Position', parse_.Results.Colorbar_unit)
    end
end