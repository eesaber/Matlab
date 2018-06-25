function Plotsetting_NorthSea1(Clim, sub_map, varargin)
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'string'},{});
    addParameter(parse_,'Colorbar_unit','',validationFcn_1_);
    parse(parse_,varargin{:})
    
    set(gca,'Ydir','normal','Clim',Clim)
    colormap jet; colorbar
    xlabel('Azimuth (km)', 'Interpreter', 'latex')
    ylabel('Range (km)', 'Interpreter', 'latex')
    if ~strcmp(parse_.Results.Colorbar_unit, '')
        title(colorbar,parse_.Results.Colorbar_unit,'Position', [40 -70])
    end
end