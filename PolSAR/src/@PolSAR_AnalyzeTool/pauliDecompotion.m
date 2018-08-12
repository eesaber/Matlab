function pauliDecompotion(obj, R, G, B, varargin)
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'char'},{});
    
    addParameter(parse_,'Filename','',validationFcn_1_);
    addParameter(parse_,'Method','2sigma',validationFcn_1_);
    parse(parse_,varargin{:})
    %% Pauli-decomposition
    method = '2sigma';
    figure
    if strcmp(method, '2sigma')
        % Show the Pauli vector with scaleing each component to [0, 255].
        image(cat(3, reshape(Stat(10*log10(R)), size(R)), ......
            reshape(Stat(10*log10(G)), size(G)), reshape(Stat(10*log10(B)), size(B))))
    else
        % Show the Pauli vector based on its power without scaleing each
        % component to [0, 255].
        up_ = 10; low_ = -35;
        Pauli = zeros([size(R), 3]);	
        % |S_vv - S_hh|^2 -> double bounce scattering 
        t_p = 10*log10(R);	
        t_p(t_p < low_) = low_;
        t_p(t_p > up_ ) = up_;
        Pauli(:,:,1) = (t_p-low_)/(up_-low_);	
        % |S_hv|^2 -> volume scattering
        t_p= 10*log10(G);
        t_p(t_p < low_) = low_;
        t_p(t_p > up_ ) = up_;
        Pauli(:,:,2) = (t_p-low_)/(up_-low_);
        % |S_vv + S_hh|^2 -> single scattering
        t_p = 10*log10(B);
        t_p(t_p < low_) = low_;
        t_p(t_p > up_ ) = up_;
        Pauli(:,:,3) = (t_p-low_)/(up_-low_);
        image(Pauli)
    end 
    set(gca,'Ydir','normal')   
    if numel(parse_.Results.Filename)
        plot_para('Maximize',true,'Filename', parse_.Results.Filename);
    end    
end
