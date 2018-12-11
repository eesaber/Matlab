classdef PolSAR_AnalyzeTool < handle
    % PolSAR_AnalyzeTool is the class which wraps the polarimetric
    % analysis method.
    %
    % Syntax:
    %	* obj = PolSAR_AnalyzeTool(image_size, carrier, plotSetting, varargin)
    %
    % Inputs:
    %	* image_size  : A 2x1 array which specify the dimension of image.
    %	* carrier     : A char array which specify the platform name, including 
    %                   'ERS-2', 'ALOS PALSAR' and 'UAVSAR'.
    %	* plotSetting : A function pointer.
    %
    %  Name-Value Pair Arguments:
    %	* 'inputDataDir'  - directory of input file
    %	* 'outputDataDir' - directory of output file
    %
    % Outputs:
    %   Return a PolSAR_AnalyzeTool object.
    %
    % Other m-files required:
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    properties
        OUTPUT_PATH;
        INPUT_PATH;
        IMAGE_SIZE;
        PLATFORM;
        POW_RANGE = [-25 -5];
        plotSetting;
        IS_BIGFILE = false;
        hh_hh; hv_hv; vv_vv; hh_hv; hh_vv; hv_vv;
        T_11; T_22; T_33; T_12; T_13; T_23;
    end
    methods
        function obj = PolSAR_AnalyzeTool(image_size, carrier, plotSetting, varargin)
            parse_ = inputParser;
            validationFcn_1_ = @(x) validateattributes(x,{'char'},{'nonempty'});
            validationFcn_2_ = @(x) validateattributes(x,{'logical'},{});
            addParameter(parse_,'inputDataDir', '',validationFcn_1_);
            addParameter(parse_,'outputDataDir','',validationFcn_1_);
            addParameter(parse_,'null', false, validationFcn_2_);
            parse(parse_,varargin{:})
            if parse_.Results.null 
                return;
            end

            if strcmp(parse_.Results.inputDataDir,'')
                obj.INPUT_PATH = uigetdir('/media/akb/2026EF9426EF696C/raw_data',...
                    'Select the folder which contain the polarimetric file');
            else
                obj.INPUT_PATH = parse_.Results.inputDataDir;
            end
            if strcmp(parse_.Results.outputDataDir,'')
                obj.OUTPUT_PATH = uigetdir('/home/akb/Code/Matlab/PolSAR',...
                    'Select the folder for saving output files');
            else
                obj.OUTPUT_PATH = parse_.Results.outputDataDir;
                if exist(obj.OUTPUT_PATH) == 0
                    mkdir(obj.OUTPUT_PATH)
                end
            end
            obj.plotSetting = plotSetting;
            obj.IMAGE_SIZE = image_size;
            obj.PLATFORM =  carrier;
            obj.readPolsarData()
            
        end
        % class function declearation
        readPolsarData(obj);
        %
        varargout = logCumulant(obj)
        logCumulantDiagram(obj, kai_2, kai_3)
        [X] = HMRF(obj, im, labels, nColors)
        [labels] = imageKCluster(obj, im, nColors)
        showLabels(obj, labels, nColors)
        varargout = myFCM(obj, x, varargin);
        varargout = mySVM(obj, x, y);
        y = myMaxPooling(obj, x, r_size, c_size);
        %
        [H, alpha_bar] = eigenDecomposition(obj, Calculate, Filename, varargin)
        pauliDecomposition(obj, varargin)
        fourComponentDecomposition(obj, varargin)
        HAlphaDiagram(obj, H, alpha_bar, varargin)
        dist_ = wishartDist(obj, x, y)
        %
        coh2cov(obj)
        cov2coh(obj)
        setCov2Zero(obj)
        setCoh2Zero(obj)
        theta = getPolAngle(obj, ang_range)
        %
        RCS(obj, varargin)
        isSNR(obj, NSZE, varargin)
        paraMoisture(obj, x, y)
        paraRatioVVHH(obj, x, y)
        paraRoughness1(obj, x, y)
        paraRoughness2(obj, x, y)
        co_pol = paraCoPolCorrelation(obj, x, y)
        %
        geocode(obj)
    end
end