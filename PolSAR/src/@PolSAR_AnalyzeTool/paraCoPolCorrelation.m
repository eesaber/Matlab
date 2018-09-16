function paraCoPolCorrelation(obj, x, y)
    % PARACOPOLCORRELATION Co-polarization correlation
    %
    % Syntax:  PARACOPOLCORRELATION(x,y)
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

    disp('S_{hh} S_{vv}^* / \sqrt{σ_vv σ_hh}')
    f_name = 'para_copolcorrelation';
    im = (angle(conv2(obj.hh_vv./sqrt(obj.hh_hh.*obj.vv_vv), mask, 'same')));
    %im = 180/pi*angle(conv2(obj.hh_vv./sqrt(obj.hh_hh.*obj.vv_vv), mask, 'same'));
    %im = 10*log10(obj.hh_hh);
    im = rescale(im, 0, 1);
    im = histeq(im);
    im = adapthisteq(im);
    
    im = entropyfilt(im);
    im = imgaussfilt(im, 8);
    
    im = rescale(im, 0, 1);
    im = histeq(im);
    %im = imsharpen(im);
    
    figure
    imagesc(im)
    %obj.plotSetting([0 1])
    plot_para('Maximize',true, 'Ratio', [4 3 1], 'Filename',[obj.OUTPUT_PATH '/' f_name])

end