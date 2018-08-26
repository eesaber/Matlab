function paraMoisture(obj, x, y)
    %% correlation ratio (T_22 + T_33)/ T_11
    disp('Soil moisture by (T_22 + T_33)/ T_11')
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
    f_name = 'para_moisture';
    figure
    imagesc(conv2((obj.T_22 + obj.T_33)./obj.T_11, mask, 'same'))
    obj.plotSetting([0 1])
    plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH '/' f_name])
end