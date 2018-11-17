function paraRoughness1(obj, x, y)
    % PARAROUGHNESS1 Correlation ratio |T_11|/sqrt{T_11 T_22}
    %
    % Syntax: PARAROUGHNESS1(x,y)
    %
    % Inputs:
    %    x  - row number of boxcar filter with default value 3.
    %    y  - column number of boxcar filter with default value 3.
    %
    % Other m-files required: plot_para
    % Subfunctions: none
    % MAT-files required: none
    %
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    disp('Surface roughness by \frac{T_11}{\sqrt{T_11 T_22}}')
    switch nargin
    case 3
        mask = ones(x,y)/(x*y);
    case 2
        mask = ones(x,3)/(x*3);
    otherwise
        mask = ones(3,3)/(3*3);
    end
    if numel(obj.T_11)==0 || numel(obj.T_22)==0 || numel(obj.T_12)==0
        error('T_11 can not be empty')
    end
    f_name = 'para_roughness1';
    figure
    imagesc(conv2(abs(obj.T_12)./sqrt(obj.T_11.* obj.T_22), mask, 'same'))
    obj.plotSetting([0 1])
    plot_para('Maximize',true,'Ratio', [4 3 1], 'Filename',[obj.OUTPUT_PATH '/' f_name])
end