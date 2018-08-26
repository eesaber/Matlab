classdef PolSAR_AnalyzeTool < handle
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
            addParameter(parse_,'inputDataDir','',validationFcn_1_);
            addParameter(parse_,'outputDataDir','',validationFcn_1_);
            addParameter(parse_,'null', false,validationFcn_2_);
            parse(parse_,varargin{:})
            if parse_.Results.null 
                return
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
            cd(obj.OUTPUT_PATH)
        end
        % class function declearation
        readPolsarData(obj);
        [H, alpha_bar] = eigenDecomposition(obj, Calculate, Filename, varargin)
        [kai_1, kai_2, kai_3] = logCumulant(obj)
        logCumulantDiagram(obj)
        segmentation(obj)
        pauliDecomposition(obj, varargin)
        fourComponentDecomposition(obj, varargin)
        HAlphaDiagram(obj, H, alpha_bar, varargin)
        %
        RCS(obj, varargin)
        isSNR(obj, NSZE, varargin)
        %
        coh2cov(obj)
        cov2coh(obj)
        setCov2Zero(obj)
        setCoh2Zero(obj)
        getPolAngle(obj, T, plot_set, ang_range)
        %
        paraMoisture(obj, x, y)
        paraRatioVVHH(obj, x, y)
        paraRoughness1(obj, x, y)
        %
        geocode(obj)
   end
end