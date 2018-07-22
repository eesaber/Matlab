function Pauli_decomp(R, G, B, varargin)
    parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'numeric'},{'nonempty'});
    validationFcn_2_ = @(x) validateattributes(x,{'logical'},{});
    validationFcn_3_ = @(x) validateattributes(x,{'char'},{});
	addParameter(parse_,'Contour',[],validationFcn_1_);
    addParameter(parse_,'Saibu',false,validationFcn_2_);
    addParameter(parse_,'Filename','',validationFcn_3_);
	parse(parse_,varargin{:})
    %% Pauli-decomposition	
    image(cat(3, reshape(stat(R), size(R)), reshape(stat(G), size(G)), ......
        reshape(stat(B), size(B)))) % Explicit explanation is in stat() 

    if numel(parse_.Results.Contour) ~= 0
        hold on 
        [rr1,rr2] = contour(parse_.Results.Contour,100:100:500,'LineColor','y','Linewidth',1,'ShowText','on');
        clabel(rr1,rr2,'Color','w','Fontsize',18)
        hold off
    end
    if numel(parse_.Results.Filename)
        plot_para('Maximize',true,'Filename', parse_.Results.Filename);
    end
    
    if parse_.Results.Saibu
        %% Plot the dominant channel
        dum = ones(size(R));
        figure
        imagesc(dum.*(R>G).*(R>B))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('Azimuth (pixel)', 'Fontsize', 40)
        ylabel('Range (pixel)', 'Fontsize', 40)
        if numel(parse_.Results.Filename)
            plot_para('Maximize',true,'Filename',[parse_.Results.Filename, '_r'])
            movefile([parse_.Results.Filename, '_r.jpg'],  'output/')
        end
        figure
        imagesc(dum.*(G>R).*(G>B))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('Azimuth (pixel)', 'Fontsize', 40)
        ylabel('Range (pixel)', 'Fontsize', 40)
        if numel(parse_.Results.Filename)
            plot_para('Maximize',true,'Filename',[parse_.Results.Filename, '_g'])
            movefile([parse_.Results.Filename, '_g.jpg'],  'output/')
        end
        figure
        imagesc(dum.*(B>R).*(B>G))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('Azimuth (pixel)', 'Fontsize', 40)
        ylabel('Range (pixel)', 'Fontsize', 40)
        if numel(parse_.Results.Filename)
            plot_para('Maximize',true,'Filename',[parse_.Results.Filename, '_b'])
            movefile([parse_.Results.Filename, '_b.jpg'],  'output/')
        end

    end
end
function [p] = stat(A)
    % STAT function is used to scale the image.
    % USAGE:
    % [p] = STAT(A), A is the image matrix and p is the 1-d array, which 
    % can be recovered to image size by using RESHAPE(p, SIZE(A)). The
    % return value is unsigned 8-bit integer.
    % ALGORITHM:
    % 1. Turn matrix A to dB unit.
    % 2. Determines the mean and standard deviation of whole image.
    % 3. Map range [μ-2σ, μ+2σ] to [0 255], set the value which is smaller
    %    than μ-2σ to μ-2σ; and set the value which is larger than μ+2σ to
    %    μ+2σ.
    %....................................................................%
    
    p = reshape(10*log10(A),[1,numel(A)]);
    p_sigma = std(p(~isinf(p)),'omitnan');
    p_mean = mean(p(~isinf(p)),'omitnan');
    p(p > p_mean+2*p_sigma) = p_mean+2*p_sigma;
    p(p < p_mean-2*p_sigma) = p_mean-2*p_sigma;
    p = uint8((p - (p_mean-2*p_sigma))*255/(4*p_sigma));
end