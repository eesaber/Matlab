function Plotsetting_NorthSea2(Clim, varargin)
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'string'},{});
    validationFcn_2_ = @(x) validateattributes(x,{'logical'},{});
    addParameter(parse_,'Colorbar_unit','',validationFcn_1_);
    addParameter(parse_,'Submap',0,validationFcn_2_);
    parse(parse_,varargin{:})
    if parse_.Results.Submap == 0
        set(gca,'XTick',0:1000:7000, 'Xticklabel',cellstr(int2str((0:1000*7:7000*7)'/1000))')
        set(gca,'YTick',0:600:3300, 'Yticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
    else
        %set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))')
        set(gca,'YTick',1:600:3300, 'Yticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
    end
    set(gca,'Ydir','normal','Clim',Clim)
    colormap jet; colorbar
    xlabel('Azimuth (km)', 'Interpreter', 'latex')
    ylabel('Range (km)', 'Interpreter', 'latex')
    if ~strcmp(parse_.Results.Colorbar_unit, '')
        title(colorbar,parse_.Results.Colorbar_unit,'Position', [40 -70])
    end
end