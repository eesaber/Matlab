function Plotsetting_dummy(Clim, varargin)
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'string'},{});
    validationFcn_2_ = @(x) validateattributes(x,{'logical'},{});
    addParameter(parse_,'Colorbar_unit','',validationFcn_1_);
    addParameter(parse_,'Submap',0,validationFcn_2_);
    parse(parse_,varargin{:})    
    set(gca,'Ydir','normal','Clim',Clim)
    colormap jet; colorbar
   
end