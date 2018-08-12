function RCS(obj, varargin)
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'logical'},{});
    addParameter(parse_, 'Span', 0, validationFcn_1_);
    addParameter(parse_, 'vv',   0, validationFcn_1_);
    addParameter(parse_, 'hh',   0, validationFcn_1_);
    addParameter(parse_, 'hv',   0, validationFcn_1_);
    parse(parse_,varargin{:})
    if parse_.Results.Span 
        figure
        imagesc(10*log10(obj.vv_vv + obj.hh_hh + 2*obj.hv_hv))
        obj.plotSetting(obj.POW_RANGE)
        colormap gray
        plot_para('Maximize',true,'Filename','Span')
    end
    if parse_.Results.vv
        figure
        imagesc(10*log10(obj.vv_vv))
        obj.plotSetting(obj.POW_RANGE)
        colormap gray
        plot_para('Maximize',true,'Filename','sigma_vv')
    end
    if parse_.Results.hh
        figure
        imagesc(10*log10(obj.hh_hh))
        obj.plotSetting(obj.POW_RANGE)
        colormap gray
        plot_para('Maximize',true,'Filename','sigma_hh')
    end
    if parse_.Results.hv
        figure
        imagesc(10*log10(obj.hv_hv))
        obj.plotSetting([-30 -15])
        colormap gray
        plot_para('Maximize',true,'Filename','sigma_hv')
    end
end