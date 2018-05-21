function Plotsetting_GOM1(Clim, varargin)
    parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'numeric'},{});
	addParameter(parse_,'Colorbar_unit',[],validationFcn_1_);
	parse(parse_,varargin{:})  
    set(gca,'Ydir','normal','Clim',Clim)
    set(gca,'XTick',0:5000:30000, 'Xticklabel',cellstr(int2str((0:5000*7:30000*7)'/1000))')
    set(gca,'YTick',0:600:3300, 'Yticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
    colormap jet; colorbar
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    
    if parse_.Results.Colorbar_unit
        title(colorbar,'(dB)','Position', parse_.Results.Colorbar_unit)
    end
end