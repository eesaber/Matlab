function [p] = intensityMapping(A, varargin)
    % INTENSITYMAPPING function is used to scale the image.
    % USAGE:
    % 1. [p] = INTENSITYMAPPING(A), A is the image matrix and p is the 1-d array, which 
    % can be recovered to image size by using RESHAPE(p, SIZE(A)). The
    % return value is unsigned 8-bit integer.
    % 2. [p] = INTENSITYMAPPING(A, bit), bit can be 'uint8', 'uint16' or 'uint32'. The
    % defalut data type is uint8.
    % ALGORITHM:
    % 1. Turn matrix A to dB unit.
    % 2. Determines the mean and standard deviation of whole image.
    % 3. Map range [μ-2σ, μ+2σ] to [0 255], set the value which is smaller
    %    than μ-2σ to μ-2σ; and set the value which is larger than μ+2σ to
    %    μ+2σ.
    %....................................................................%
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'char'},{});
    validationFcn_2_ = @(x) validateattributes(x,{'numeric'},{'numel', 2});
    addParameter(parse_,'Bit','uint8',validationFcn_1_);
    addParameter(parse_,'Method','2sigma',validationFcn_1_);
    addParameter(parse_,'Range',[0 1],validationFcn_2_);
    parse(parse_,varargin{:})

    if strcmp(parse_.Results.Bit, 'uint8')
        scale = 2^8-1;
        f = @(x) uint8(x);
    elseif strcmp(parse_.Results.Bit, 'uint16')
        scale = 2^16-1;
        f = @(x) uint16(x);
    elseif strcmp(parse_.Results.Bit, 'uint32')
        scale = 2^32-1;
        f = @(x) uint32(x);
    else
        error('Available Bit option: uint8, uint16, uint32')
    end

    p = A(:);
    if strcmp(parse_.Results.Method, '2sigma')
        p_sigma = std(p(~isinf(p)),'omitnan');
        p_mean = mean(p(~isinf(p)),'omitnan');
        p(p > p_mean+2*p_sigma) = p_mean+2*p_sigma;
        p(p < p_mean-2*p_sigma) = p_mean-2*p_sigma;
        p = f((p - (p_mean-2*p_sigma))*scale/(4*p_sigma));
    elseif strcmp(parse_.Results.Method, 'truncate')
        p(p > parse_.Results.Range(2)) = parse_.Results.Range(2);
        p(p < parse_.Results.Range(1)) = parse_.Results.Range(1);
        p = f((p - parse_.Results.Range(1))/ ...
            (parse_.Results.Range(2) - parse_.Results.Range(1))*scale);
    end
end