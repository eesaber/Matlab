classdef SeaIce < PolSAR_AnalyzeTool
    % SEAICE is inherited from PolSAR_AnalyzeTool.
    %
    % Syntax:
	%	* obj = SEAICE(image_size, carrier, plotSetting, varargin)
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
    %   Return a SEAICE object.
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
        H, alpha_bar;
        kai_1; kai_2; kai_3; kai_4;% log-cumulant;
        im; % log-cumulant plot
        groundtruth;
        y_hat;
    end
    methods
        % Constructor
        function obj = SeaIce(image_size, carrier, plotSetting, varargin)
            obj = obj@PolSAR_AnalyzeTool(image_size, carrier, plotSetting,...
                   varargin{:});
        end

        % Extend the base class function 
        function readPolsarData(obj)
            % READPOLSARDATA is extended from PolSAR_AnalyzeTool.
            % Load data as the base class does, and clear the coherent matrix.
            readPolsarData@PolSAR_AnalyzeTool(obj);
            obj.setCoh2Zero();
        end

        function logCumulant(obj)
            [obj.kai_1, obj.kai_2, obj.kai_3, obj.kai_4] = ...
                        logCumulant@PolSAR_AnalyzeTool(obj);
        end

        function varargout = eigenDecomposition(obj, Calculate, Filename, varargin)
            % parse output arguments
            minArgs=0;
            maxArgs=2;
            nargoutchk(minArgs,maxArgs)
            % run eigen-decomposition
            [obj.H, obj.alpha_bar] = eigenDecomposition@PolSAR_AnalyzeTool(obj,...
                Calculate, Filename, varargin{:});
            % output
            if nargout>=1, varargout{1} = obj.H;
                if nargout>=2, varargout{2} = obj.alpha_bar;
                end
            end
        end
            
        function imageKCluster(obj, nColors)
            temp = imageKCluster@PolSAR_AnalyzeTool(obj, obj.im, nColors);
            creatLogCumulantRGB(obj);
            %obj.y_hat = padarray(temp, [1,1], -1,'both');
            obj.y_hat = temp;
        end

        function HMRF(obj, nColors)
            temp_label = obj.y_hat;
            temp_label_MRF = HMRF@PolSAR_AnalyzeTool(obj, obj.im, temp_label, nColors);
            obj.y_hat_MRF = padarray(temp_label_MRF, [1,1], -1,'both');
        end
        
        %%
        creatLogCumulantRGB(obj);
        varargout = myDataAugmentation(obj, im, label, varargin);
        varargout = generateImage4Classification(obj, input_vector, varargin)
    end
end
%{
classdef Sub < Super
   methods
      function foo(obj)
         % preprocessing steps
          ...
         foo@Super(obj);
         % postprocessing steps
          ...
      end
   end
end
%}