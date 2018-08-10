classdef PolSAR_AnalyzeTool < handle
   properties
        outputPath;
        inputPath;
        polsarFileName;
        image_size;
   end
   methods
        function obj = PolSAR_AnalyzeTool()
            obj.polsarFileName = '';
            
        end
        function set.outputPath(obj, value)
            if(ischar(value))
                obj.outputPath = value;
            else
                error('Output path must be a char array')
            end
        end
        function set.inputPath(obj, value)
            obj.inputPath = value;
        end
        % class function declearation
        output = eigenDecomp(varargin) 
        
   end
end