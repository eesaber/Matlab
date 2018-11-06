function isSNR(obj, NESZ, varargin)
    % ISSNR is used to check if the signal is above the signal-to-noise
    % ratio.
    %
    % Syntax: eigenDecomp(NESZ, NAME, VALUE)
    % Specifies function properties using one or more Name,Value pair arguments. 
    % Name-value pair settings apply to all the lines plotted. 
    %
    % Inputs:
    %    NESZ      - It is noise equivalent sigam zero (NESZ). NESZ should be 
    %                scalar or a matrix with same size as image.
    % Name-Value Pair Arguments:
    %    ('vv', s) - The pixel above SNR thershold is marked. If s is ture, 
    %                function will show the figure.
    %    ('hh', s) - The pixel above SNR thershold is marked. If s is ture, 
    %                function will show the figure.
    %    ('hv', s) - The pixel above SNR thershold is marked. If s is ture, 
    %                function will show the figure.
    %
    % Example: 
    %    1. ISSNR(NESZ, 'vv', true)
    %       It will plot a 2-D image which mark the pixel satisfy SNR with white point
    %       and mark the pixel below SNR 
    %
    % Other m-files required: plot_para
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'logical'},{});
    addParameter(parse_,'vv', 0,validationFcn_1_);
    addParameter(parse_,'hh', 0,validationFcn_1_);
    addParameter(parse_,'hv', 0,validationFcn_1_);
    parse(parse_,varargin{:})

    if parse_.Results.vv
        figure
        imagesc(10*log10(obj.vv_vv)-NESZ>6)
        %title('\sigma_{vv}')
        obj.plotSetting([0,1])
        colormap gray
        colorbar off
        set(gca,'Visible','on')
        plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH, '/SNR_sigma_vv'])
    end
    if parse_.Results.hh
        figure
        imagesc(10*log10(obj.hh_hh)-NESZ>6)
        %title('\sigma_{hh}')
        colormap gray
        plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH, '/SNR_sigma_hh'])
    end
    if parse_.Results.hv
        figure
        imagesc(10*log10(obj.hv_hv)-NESZ>6)
        %title('\sigma_{hv}')
        colormap gray
        plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH, '/SNR_sigma_hv'])
    end
    
end