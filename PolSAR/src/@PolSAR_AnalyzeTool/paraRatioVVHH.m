function paraRatioVVHH(obj, x , y)
    %% σ_vv/ σ_hh ratio
    switch nargin
        case 2
            mask = ones(x,y)/(x*y);
        case 1
            mask = ones(x,3)/(x*3);
        otherwise
            mask = ones(9,9)/(9*9);
    end
    if numel(obj.vv_vv)==0 || numel(obj.hh_hh)==0
        error('vv_vv can not be empty')
    end

    disp('σ_vv/ σ_hh')
    f_name = 'para_hvratio_g';
    figure
    imagesc(10*log10(conv2(obj.vv_vv./ obj.hh_hh, mask, 'same')))
    obj.plotSetting([-3 3])
    plot_para('Maximize',true,'Filename',f_name)
end