    function paraRatioVVHH(obj, x, y)
    % PARARATIOVVHH  Co-polarization ratio σ_vv/ σ_hh 
    %
    % Syntax: PARARATIOVVHH(x,y)
    %
    % Inputs:
	%    x  - row number of boxcar filter with default value 3.
	%    y  - column number of boxcar filter with default value 3.
    %
    % Other m-files required: plot_para
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    switch nargin
    case 3
        mask = ones(x,y)/(x*y);
    case 2
        mask = ones(x,3)/(x*3);
    otherwise
        mask = ones(3,3)/(3*3);
    end
    if numel(obj.vv_vv)==0 || numel(obj.hh_hh)==0
        error('vv_vv can not be empty')
    end

    disp('σ_vv/ σ_hh')
    f_name = 'para_hvratio';
    figure
    imagesc(10*log10(conv2(obj.vv_vv./ obj.hh_hh, mask, 'same')))
    
    obj.plotSetting([-3 3])
    plot_para('Maximize',true, 'Filename',[obj.OUTPUT_PATH '/' f_name])
end