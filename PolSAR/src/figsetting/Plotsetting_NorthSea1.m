function Plotsetting_NorthSea1(Clim, varargin)
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'string'},{});
    validationFcn_2_ = @(x) validateattributes(x,{'logical'},{});
    addParameter(parse_,'Colorbar_unit','',validationFcn_1_);
    addParameter(parse_,'Submap',0,validationFcn_2_);
    parse(parse_,varargin{:})
    
    set(gca,'Ydir','normal','Clim',Clim)
    colormap jet; colorbar
    xlabel('Azimuth (km)', 'Interpreter', 'latex')
    ylabel('Range (km)', 'Interpreter', 'latex')
    if ~strcmp(parse_.Results.Colorbar_unit, '')
        title(colorbar,parse_.Results.Colorbar_unit,'Position', [40 -70])
    end
end