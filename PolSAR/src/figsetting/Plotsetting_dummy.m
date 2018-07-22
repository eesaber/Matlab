function Plotseeting_dummy(Clim, sub_map, varargin)
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'string'},{});
    addParameter(parse_,'Colorbar_unit','',validationFcn_1_);
    parse(parse_,varargin{:})    
    set(gca,'Ydir','normal','Clim',Clim)
    colormap jet; colorbar
    
end