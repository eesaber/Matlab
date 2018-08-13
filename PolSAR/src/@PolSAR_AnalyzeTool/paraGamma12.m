function paraGamma12(obj, x, y)
    %% correlation ratio |T_11|/sqrt{T_11 T_22}
    switch nargin
    case 2
        mask = ones(x,y)/(x*y);
    case 1
        mask = ones(x,3)/(x*3);
    otherwise
        mask = ones(9,9)/(9*9);
    end
    if numel(obj.T_11)==0 || numel(obj.T_22)==0 || numel(obj.T_12)==0
        error('T_11 can not be empty')
    end
    f_name = 'para_gamma12';
    figure
    imagesc(conv2(abs(obj.T_12)./sqrt(obj.T_11.* obj.T_22), mask, 'same'))
    obj.plotSetting([0 1])
    plot_para('Maximize',true,'Filename',f_name)
end